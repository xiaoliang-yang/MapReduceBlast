/**
 * modified at 2011-7-2
 * added IS_LAST_SPLIT_BYTES field
 * 
 * MapReduce BLAST
 * @author Yang Xiao-liang
 * @email yangxiaoliang2006@gmail.com
 * @copyright June, 2011
 */
package edu.nju.cs.mapreduceblast;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.util.ReflectionUtils;

public class Fasta2SequenceFileConverter {
	// NCBI BLAST sequence database format, each sequence begins with ">"
	public static final char NEW_SEQUENCE_IDENTIFIER = '>';
	// the length of returned the sequence in bp (base pair)
	// must be even
	// default value: 10Mbp
	public static final int bpAmountPerLineDefault = 10 * 1024 * 1024;
	// the default length of overlapped flank
	// bust be even
	// default value: 1Mbp
	public static final int bpAmountOverlapDefault = 1024 * 1024;
	// starting offset field, int type, 4bytes
	public static final int OFFSET_FIELD_BYTES = 4;
	// sequcne length field, int type, 4bytes
	public static final int SEQ_LENGTH_FIELD_BYTES = 4;
	// overlap length field, int type, 4bytes
	public static final int OVERLAP_LENGTH_FIELD_BYTES = 4;
	
	public static final int IS_LAST_SPLIT_BYTES =1;
	
	
	//the offset of sequence split in a entire sequence	
	private int start_offset=0;
	//sequence ID, as the key of the <key, value> pair
	private String seqID = null;
	
	//buffered array for the value field of a key-value pair
	private byte[] valBuf;
	//the start position of sequence split in valBuf 
	private final int valBufSeqSplitStartIndex=OFFSET_FIELD_BYTES
		+SEQ_LENGTH_FIELD_BYTES+OVERLAP_LENGTH_FIELD_BYTES
		+IS_LAST_SPLIT_BYTES;
	//point to the byte in which is ready to put DNA base 
	private int valBufIndex=valBufSeqSplitStartIndex;
	//the current bases number in valBuf
	private int counterSeqSplit_bp=0;
	
	private SequenceFile.Writer writer;
	private BufferedReader reader;
	
	private Text key = new Text();
	private BytesWritable value = new BytesWritable();
	
	// buffer for key value pair
	private String outKeyBuf;
	private byte[] outValueBuf;
	private int outValueBufLen;
	private boolean isOutKeyValueBufFull=false;
	
	/**
	 * 
	 * @param fastaInFile input fasta file
	 * @param fs file system for the output file, e.g. hdfs or others
	 * @param conf Hadoop Configuration
	 * @param fileName output file name
	 */
	public Fasta2SequenceFileConverter(File fastaInFile, FileSystem fs,
			Configuration conf, Path fileName){
	
		try {
			reader = new BufferedReader(new FileReader(fastaInFile));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			writer = new SequenceFile.Writer(fs, conf, fileName, Text.class, BytesWritable.class);
		} catch (IOException e) {
			e.printStackTrace();
		}
		// value: |start_offset|seq.length|overlap.length|sequence split|overlap flank|
		// each byte contains two bp
		// assume this.bpAmountOverlapDefautl and this.bpAmountPerLineDefault are both even numbers
		if(valBuf==null){
			valBuf = new byte[OFFSET_FIELD_BYTES+SEQ_LENGTH_FIELD_BYTES
			                  +OVERLAP_LENGTH_FIELD_BYTES+IS_LAST_SPLIT_BYTES
			                  +(bpAmountOverlapDefault+bpAmountPerLineDefault)/2];
			outValueBuf=new byte[valBuf.length];
		}
	}
	
	/**
	 * write the content of valBuf to SequenceFile format file,
	 * after this method is called, valBuf will be cleared
	 * if valBuf is empty, no action taken
	 * @return true if succeed
	 * @throws IOException 
	 */
	public boolean writeKeyValue() throws IOException{
		//if nothing to write
		if(this.counterSeqSplit_bp==0){
			return false;
		}
		//key: sequence id (e.g. 224589800) 
		//value: |start_offset (int)|seq.length (int)|overlap.length (int)|sequence split|overlap|
		
		//start offset field bytes
		byte[] startOffsetField = BlastUtils.int2bytes(this.start_offset);
		for(int i=0; i<startOffsetField.length; i++){
			this.valBuf[0+i]=startOffsetField[i];
			
		}
//		//@debug
//		System.out.println("seqID: "+this.seqID);
//		System.out.println("start_offset: "+this.start_offset+" bytes: "+Arrays.toString(startOffsetField));
//		System.out.print("bytes: [");
//		for(int i =0; i<4; i++){
//			System.out.printf("%02X ", startOffsetField[i]);
//		}
//		System.out.println("]");
			
			
		
		//sequence split length field bytes
		int seqLen = this.counterSeqSplit_bp > bpAmountPerLineDefault? 
				bpAmountPerLineDefault : this.counterSeqSplit_bp;
		byte[] seqLenField = BlastUtils.int2bytes(seqLen);
		for(int i=0; i<seqLenField.length; i++){
			this.valBuf[OFFSET_FIELD_BYTES+i]=seqLenField[i];
		}
//		//@debug
//		System.out.println("seq split length: "+seqLen+" bytes: "+Arrays.toString(seqLenField));
//		//System.out.println("bytes in record: "+Arrays.toString(Arrays.copyOfRange(this.valBuf, 4, 8)));
//		System.out.print("bytes: [");
//		for(int i =0; i<4; i++){
//			System.out.printf("%02X ", seqLenField[i]);
//		}
//		System.out.println("]");
		
		
		//overlap length field
		int overlapLen = this.counterSeqSplit_bp - bpAmountPerLineDefault > 0 ?
				this.counterSeqSplit_bp- bpAmountPerLineDefault : 0;
		byte[] overlapField = BlastUtils.int2bytes(overlapLen);
		for(int i=0; i<overlapField.length; i++){
			this.valBuf[SEQ_LENGTH_FIELD_BYTES+OFFSET_FIELD_BYTES+i]=overlapField[i];
		}
//		//@debug
//		System.out.println("overlap length: "+overlapLen+" bytes: "+Arrays.toString(overlapField));
//		//System.out.println("bytes in record: "+Arrays.toString(Arrays.copyOfRange(this.valBuf, 8, 12)));
//		System.out.print("bytes: [");
//		for(int i =0; i<4; i++){
//			System.out.printf("%02X ", overlapField[i]);
//		}
//		System.out.println("]");
		
		
		// if a key-value pair is buffered already
		if(isOutKeyValueBufFull){
			// if the buffered key-value record is last split of the sequence
			if(!seqID.equals(outKeyBuf)){
				outValueBuf[0+OFFSET_FIELD_BYTES+SEQ_LENGTH_FIELD_BYTES
			                  +OVERLAP_LENGTH_FIELD_BYTES]=0x01;
			}else{
				outValueBuf[0+OFFSET_FIELD_BYTES+SEQ_LENGTH_FIELD_BYTES
			                  +OVERLAP_LENGTH_FIELD_BYTES]=0x00;
			}
			//write the buffered key-value pair to the HDFS SequenceFile file
			this.key.set(outKeyBuf);
			this.value.set(outValueBuf, 0, outValueBufLen);
			this.writer.append(this.key,this.value);
			
			// buffer the current key-value pair
			// copy key-value pair to out key value buffer
			this.outKeyBuf = this.seqID;
			for(int i=0; i<this.valBufIndex; i++){
				this.outValueBuf[i]=valBuf[i];
			}
			this.outValueBufLen=this.valBufIndex;
			this.isOutKeyValueBufFull=true;
			
		}else{
			// copy key-value pair to out key value buffer
			this.outKeyBuf = this.seqID;
			for(int i=0; i<this.valBufIndex; i++){
				this.outValueBuf[i]=valBuf[i];
			}
			this.outValueBufLen=this.valBufIndex;
			this.isOutKeyValueBufFull=true;
		}
		
		//clear valBuf
		this.valBufIndex = this.valBufSeqSplitStartIndex;
		
		//update start_offset
		this.start_offset += bpAmountPerLineDefault; //this.current_bp_offset;
		
//		//@debug
//		for(int i=0; i<this.valBuf.length; i++){
//			this.valBuf[i]=0;
//		}
		return true;
	}
	
	// called only when the last line of input multi-fasta file is processed  
	// suppose that the buffered key-value record is last split of the sequence
	private void flashOutKeyValueBuf() throws IOException{
		if(isOutKeyValueBufFull){
			// suppose that the buffered key-value record is last split of the sequence
			outValueBuf[0+OFFSET_FIELD_BYTES+SEQ_LENGTH_FIELD_BYTES
			                  +OVERLAP_LENGTH_FIELD_BYTES]=0x01;
			//write the buffered key-value pair to the HDFS SequenceFile file
			this.key.set(outKeyBuf);
			this.value.set(outValueBuf, 0, outValueBufLen);
			this.writer.append(this.key,this.value);
			isOutKeyValueBufFull=false;
		}
		
	}
	/**
	 * append sequence split line to seq split field of valBuf
	 * if valBuf is full, method writeKeyValue() will be called internally 
	 * @param seqSplitLine
	 * @return
	 * @throws Exception 
	 */
	public boolean appendToBuf(String seqSplitLine) throws Exception{
		byte baseByte;
		boolean leftHalf;
		for(int seqSplitLineIndex=0; seqSplitLineIndex<seqSplitLine.length(); seqSplitLineIndex++){
			if(leftHalf = this.counterSeqSplit_bp%2 ==0){
				// true : false;
				baseByte = BlastUtils.encode_bp(seqSplitLine.charAt(seqSplitLineIndex), true);
				this.valBuf[this.valBufIndex] &= (byte)0x0F;
				this.valBuf[this.valBufIndex] |= baseByte;
			}else{
				baseByte = BlastUtils.encode_bp(seqSplitLine.charAt(seqSplitLineIndex), false);
				this.valBuf[this.valBufIndex] &= (byte)0xF0;
				this.valBuf[this.valBufIndex] |= baseByte;
			}
			this.counterSeqSplit_bp++;
			//if right half byte, then index will proceed in the next byte
			if(!leftHalf){
				this.valBufIndex++;
			}
			
			//if valBuf is full, invoke the write operation
			if(this.counterSeqSplit_bp == 
				bpAmountPerLineDefault+bpAmountOverlapDefault){
				//write key-value
				this.writeKeyValue();
				//the valBuf is cleared here
				
				//copy the overlap to seqSplit bytes
				int overlapCount = counterSeqSplit_bp - bpAmountPerLineDefault;
				int overlapIndex = valBufSeqSplitStartIndex+bpAmountPerLineDefault/2;
				for(int bpCopyed=0; bpCopyed<=overlapCount; bpCopyed++){
					if(bpCopyed%2==0){
						continue;
					}
					this.valBuf[this.valBufIndex++]=this.valBuf[overlapIndex++];
				}
				//update the counter
				this.counterSeqSplit_bp=overlapCount;
			}
		}
		return true;
	}
	
	
	/**
	 * convert a DNA sequences file in multi-fasta format to HDFS SequenceFile
	 * format
	 * read the fasta file line by line, and convert the characters into 4bit-code,
	 * stored them in a byte array, then write to SequcenceFile Record when it reaches
	 * the default length of the recored (10Mbp sequence+1Mbp overlap) 
	 * 
	 * DNA fasta format exaple: file: ref_GRCH37_chr37_1.fa content as follows:
	 * >gi|224589800|ref|NC_000001.10| Homo sapiens chromosome 1, GRCh37 primary reference assembly
	 * NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	 * TATTAGTGATTTGGGCTGGGGCCTGGCCATGTGTATTTTTTTAAATTTCCACTGATGATTTTGCTGCATG
	 * GCCGGTGTTGAGAATGACTGCGCAAATTTGCCGGATTTCCTTTGCTGTTCCTGCATGTAGTTTAAACGAG
	 * 
	 * DNA SequenceFile format: 
	 * stored as <key, value> pair records,
	 * a long sequence is split into fixed-length record, 
	 * the default length is 10Mbp, with 1Mbp-flank overlapping 
	 * key: sequence id (e.g. 224589800) 
	 * value: |start_offset|seq.length|overlap.length|is last split flag|sequence split|overlap|
	 * {int start_offset in sequence}, 
	 * {int seq. length} default 10Mbp, 
	 * {int overlap length} 
	 * {byte[] seq split} 
	 * {byte[] overlap seq} 
	 * {byte is last split flag}
	 * DNA >>1 bp
	 * represented by 4 bits, 
	 * A: 0000 T: 0001 C: 0010 G: 0011 
	 * N: 0100 (unknown base)
	 * 
	 * @param fastaInFile
	 * @param fs
	 * @param conf
	 * @param fileName
	 * @return true if succeed, false else
	 * @throws Exception 
	 * @throws Exception
	 */
	public boolean convert2HadoopSequenceFile() throws Exception {
		// read sequence split lines, fill them in valBuf, 
		String seqSplit=null;
		boolean validSeq = false;
		while((seqSplit=reader.readLine()) != null){
			if(seqSplit.length()==0) {
				continue;
			}
			//the header line
			if(seqSplit.charAt(0) == NEW_SEQUENCE_IDENTIFIER){
				//a new sequence begins, write the last key-value of previous sequence if any 
				this.writeKeyValue();
				//get the sequence ID
				this.seqID = BlastUtils.getSeqID(seqSplit);
				this.start_offset=0;
				this.counterSeqSplit_bp=0;
				validSeq=true;
				continue;
			}
			//each sequence must begin with a header line, then followed with sequence split lines
			if(!validSeq){
				continue;
			}
			//here is the sequence split lines
			//fill the bases in the seqSplit in seqBuf
			//the appendToBuf() method will call the writeKeyValue() when buffer is full
			this.appendToBuf(seqSplit);
		}
		
		//here, read completed, write the last key-value if any
		this.writeKeyValue();
		flashOutKeyValueBuf();
		return true;
	}
	
	/**
	 * untested
	 * @param recordVal
	 * @return
	 */
	//TODO: untested
	public static int getSeqOffset(byte[] recordVal){
		return BlastUtils.bytes2int(recordVal, 0);
	}
	
	// get the sequence split length
	//TODO: untested
	public static int getSeqSplitLength(byte[] recordVal){
		return BlastUtils.bytes2int(recordVal, OFFSET_FIELD_BYTES);
	}
	//TODO: untested
	public static int getSeqOverlapLength(byte[] recordVal){
		return BlastUtils.bytes2int(recordVal, OFFSET_FIELD_BYTES+SEQ_LENGTH_FIELD_BYTES); 
	}
	//TODO: untested
	public static int getSeqBeginIndex(byte[] recordVal){
		return OFFSET_FIELD_BYTES+SEQ_LENGTH_FIELD_BYTES+OFFSET_FIELD_BYTES;
	}
	
	public void close() throws IOException{
		this.reader.close();
		this.writer.close();
	}

	private static void validate(String seqfile) throws Exception{
		// (File fastaInFile, FileSystem fs, Configuration conf, Path fileName)
		Configuration conf = new Configuration();
		FileSystem fs = FileSystem.get(conf);
		Path path = new Path(seqfile);
		SequenceFile.Reader reader = new SequenceFile.Reader(fs, path, conf);
		fs.open(path);
		PrintStream out = new PrintStream(new File(seqfile+".txt"));
		Writable key = (Writable)ReflectionUtils.newInstance(reader.getKeyClass(), conf);
		Writable value = (Writable)ReflectionUtils.newInstance(reader.getValueClass(), conf);
		int start_offset ;
		int seqSplitLen ;
		int overlapLen ;
		int counter_bp ;
	
		while(reader.next(key, value)){
			StringBuilder strBuilder = new StringBuilder();
			strBuilder.append(key.toString()+"@");
			byte[] bytes = ((BytesWritable)value).getBytes();
			start_offset = BlastUtils.bytes2int(bytes, 0);
			seqSplitLen = BlastUtils.bytes2int(bytes, 4);
			overlapLen = BlastUtils.bytes2int(bytes, 8);
			strBuilder.append(start_offset+" ");
			strBuilder.append(seqSplitLen+" ");
			strBuilder.append(overlapLen+" ");
			counter_bp = 0;
			for(int i=12; i< bytes.length; i++){
				strBuilder.append((char)BlastUtils.decode_bp(bytes[i], true));
				counter_bp++;
				if(counter_bp==seqSplitLen+overlapLen){
					break;
				}
				strBuilder.append((char)BlastUtils.decode_bp(bytes[i], false));
				counter_bp++;
				if(counter_bp==seqSplitLen+overlapLen){
					break;
				}
			}
			out.println(strBuilder.toString());
		}
		
	}
	public static void main(String[] args) throws Exception{
		// (File fastaInFile, FileSystem fs, Configuration conf, Path fileName)
		if(args.length < 2 || args.length > 3 ){
			System.out.println("usage: fasta2seq <input fasta file> <output file name> [-validate]");
			System.out.println("notes: \n" +
					"\tinput fasta file shuold be on local disk file system \n" +
					"\tsuch as Linux ext3, instead of HDFS\n" +
					"\toutput file will be written to the HDFS\n" +
					"\t-validate option generates a *.txt file \n" +
					"\twhich is recoverted from the seq file");
			System.exit(-1);
		}
		Configuration conf = new Configuration();
		FileSystem fs = FileSystem.get(conf);
		String file= args[0]; 
		
		File fastaFile = new File(file);
		System.out.println("reading fasta data from file: "+file);
		
		Path path = new Path(args[1]);
		System.out.println("writing SequenceFile format NDA sequences to file: "+path.toString());
		
		System.out.println("Please wait...");
		Fasta2SequenceFileConverter converter = new Fasta2SequenceFileConverter(fastaFile, fs, conf, path);
		
		converter.convert2HadoopSequenceFile();
		converter.close();
		System.out.println("convert completed!");
		if(args.length == 3 && args[2].equals("-validate")){
			System.out.println("generating validate file...");
			validate(file+".seq");
			System.out.println("validate file :"+file+".seq.txt generated.");
		}
	}

}
