/**
 * MapReduce BLAST
 * @author Yang Xiao-liang
 * @email yangxiaoliang2006@gmail.com
 * @copyright June, 2011
 */

package edu.nju.cs.mapreduceblast;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.mapreduce.Mapper;

import edu.nju.cs.mapreduceblast.automaton.Scanner;
import edu.nju.cs.mapreduceblast.automaton.TwoHit;

/**
 * last modified 2011-7-1
 * 
 * @author Yang Xiao-liang
 * created at June 2011
 * read sequence record from SequenceFile file
 * input: key: seqID(Text), value: metaInfo_seqSplit(BytesWritable)
 * 
 * for file formant details, see Fasta2SequenceFileConverter.java
 * 
 * scanning 'two-hit' in database sequences by using an AC automu
 * 
 * extending 'two-hit' to ungapped (then gapped) alignment
 */
public class BlastnMapper extends Mapper<Text, BytesWritable, 
		Text, BytesWritable>{ // Alignment> {
	
	enum GappedExtend { DP_TIMES}
	
	private boolean DEBUG=false;
	// key-value to write out
	private IntWritable outKeyAlignScore = new IntWritable();
	private Alignment outValueAlign ;

	//blastn arguments
	private int wordLen;
	
	private int match;
	private int mismatch;
	private int gap;
	
	private int twoHitDistanceA;
	
	private int ungappedScoreThreshold;
	private int gappedScoreThreshold;

	private int ungappedXDrop;
	private int gappedXDrop;

	//for query and db Sequences 
	private Scanner scanner;
	private byte[] queryBytes;
	//workspace for buffered db seq 
	private byte[] seqbytes = new 
			byte[Fasta2SequenceFileConverter.bpAmountPerLineDefault+
	        Fasta2SequenceFileConverter.bpAmountOverlapDefault];
	private int seqStartOffset ;
	private int seqSplitLen ;
	private int overlapLen ;
	private boolean isLastRecordOfSeq=false;
	
	//fields for the alignment
	private DPExtender dpExtender;
	
	//
	private HashSet<String> twoHitsHashSet = new HashSet<String>();
	private ArrayList<GapFreeExtension> filteredExt = 
			new ArrayList<GapFreeExtension>();
	
	private Text outKey = new Text();
	private BytesWritable outValue = new BytesWritable();
	/**
	 * load the scanner and query bytes file 
	 */
	@Override
	public void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		this.DEBUG = conf.getBoolean("DEBUG", false);
		
		//TODO:
		//get blastn arguments/////////////////
		this.match = conf.getInt("blastn.argument.wordLen", 11);
		this.match= conf.getInt("blastn.argument.matchScore", 2);
		this.mismatch= conf.getInt("blastn.argument.mistmatchScore",
				-3);
		this.gap = conf.getInt("blastn.argument.gapPenalty", -5);
		
		this.twoHitDistanceA = conf.getInt(
				"blastn.argument.twoHitDistanceA", 40);

		this.ungappedScoreThreshold = conf.getInt(
				"blastn.argument.ungappedScoreThreshold", 100*match);
		this.gappedScoreThreshold = conf.getInt(
				"blastn.argument.gappedScoreThreshold", 200*match);
		
		this.ungappedXDrop = conf.getInt("blastn.argument.ungappedX",
				10*match);
		this.gappedXDrop = conf.getInt("balstn.argument.gappedX",
				10*match);
		
		// load scanner object from cached file on local disk
		String scannerFile = conf.get("blastn.scanner.file");
		try {
			scanner = Scanner.readFromFile(scannerFile);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
		// load query bytes from cache file 
		String queryBytesFile = conf.get("blastn.query.byte.file");
		queryBytes = BlastUtils.readBytesFromLoalFile(queryBytesFile);
		
		// create the extender
		dpExtender = new DPExtender();
	}
	

	/**
	 * scanning 'two-hit' in database sequences by using an AC automu
	 * extending 'two-hit' to ungapped (then gapped) alignment
	 */
	@Override
	public void map(Text key, BytesWritable val, Context context) 
			throws IOException, InterruptedException{
		// extract bases in val, into a byte array,///////////// 
		// in which each byte contains one bp(encoded in 4 bits) on
		// the right half byte
		byte[] bytes = val.getBytes();
		int counter_bp = 0;
		seqStartOffset = BlastUtils.bytes2int(bytes, 0);
		seqSplitLen = BlastUtils.bytes2int(bytes, 4);
		overlapLen = BlastUtils.bytes2int(bytes, 8);
		for(int i=12; i< bytes.length; i++){
			// in val, one byte contains two bases,
			seqbytes[counter_bp] = (byte)(bytes[i] >>>4);
			counter_bp++;
			if(counter_bp==seqSplitLen+overlapLen){
				break;
			}
			seqbytes[counter_bp] = (byte)(bytes[i] & (byte)0x0F);
			counter_bp++;
			if(counter_bp==seqSplitLen+overlapLen){
				break;
			}
		}
		//scan for two-hits////////////////////
		// if is last record  of the sequence
		isLastRecordOfSeq = bytes[4+4+4]==(byte)0x01 ? true : false; 
		ArrayList<TwoHit> twoHitList=null;
		if(isLastRecordOfSeq){
			twoHitList = scanner.scan(seqbytes, 0, 
					seqSplitLen+overlapLen,
					twoHitDistanceA);
		}else{
			twoHitList = scanner.scan(seqbytes, 0, 
					seqSplitLen+wordLen, twoHitDistanceA);
		}
		//TODO:
		
		
//		System.out.println("seqID:newone");
		//extending////////////////////////////
		GapFreeExtension gapFreeExt;
		int preTwoHitQueryIndex=-10;
		int preTwoHitSubjectIndex=-10;
		twoHitsHashSet.clear();
		filteredExt.clear();
		for(TwoHit twoHit: twoHitList){
			//avoid repeated two-hit extending
			// if consecutive hits
			if((twoHit.getSecondIndexInDB_Sequence()
						==preTwoHitSubjectIndex
					&& twoHit.getSecondIndexInQuery()
						==preTwoHitQueryIndex)
					|| (twoHit.getSecondIndexInDB_Sequence()+1 
							== preTwoHitSubjectIndex 
						&& twoHit.getSecondIndexInQuery()+1 
							== preTwoHitQueryIndex)
					|| (twoHit.getSecondIndexInDB_Sequence()-1 
							== preTwoHitSubjectIndex 
						&& twoHit.getSecondIndexInQuery()-1 
							== preTwoHitQueryIndex)){
				preTwoHitSubjectIndex 
					= twoHit.getSecondIndexInDB_Sequence();
				preTwoHitQueryIndex = twoHit.getSecondIndexInQuery();
				continue;
			}
			preTwoHitQueryIndex = twoHit.getSecondIndexInQuery();
			preTwoHitSubjectIndex 
				= twoHit.getSecondIndexInDB_Sequence();
			
			//avoid gap free extension repeat
			// if repeative hits
			boolean isRepeatHit = false;
			for(GapFreeExtension ext: filteredExt){
				// if on the same diagonal and in the range of 
				// ext's extension
				int dist1 = twoHit.getSecondIndexInDB_Sequence()
						-ext.twoHit.getSecondIndexInDB_Sequence(); 
				int dist2 = twoHit.getSecondIndexInQuery()
						-ext.twoHit.getSecondIndexInQuery(); 
				if(dist1 == dist2 
						&& dist1 >= ext.leftEnd 
						&& dist1 <= ext.rightEnd){
					isRepeatHit = true;
					break;
				}
			}
			if(isRepeatHit){
				continue;
			}
			
			//gap free extending
			gapFreeExt = ungappedExtend(
					queryBytes, 0, queryBytes.length,
					seqbytes, 0, seqbytes.length, twoHit, 
					this.match, this.mismatch, this.ungappedXDrop);
			
			if(gapFreeExt.score < ungappedScoreThreshold){
				continue;
			}
			
			filteredExt.add(gapFreeExt);
		}	
		
		Collections.sort(filteredExt);
		
		// transmit seq bytes data to reducer, and do DP-alignment 
		// at reducer side
		for(int i = 0; i<10 && i<filteredExt.size(); i++){
			// emit key=(SeqID,startOffset,XBaseIndex,YBaseIndex)
			// value=(seqbytes)
			int XBaseIndex = 
					filteredExt.get(i).twoHit.getSecondIndexInQuery();
			int YBaseIndex = filteredExt.get(i)
					.twoHit.getSecondIndexInDB_Sequence();
			outKey.set(key+","+seqStartOffset+","+XBaseIndex
					+","+YBaseIndex);
			outValue.set(seqbytes, 0, seqSplitLen+overlapLen);
			
			context.write(outKey, outValue);
			
			context.getCounter(
					GappedExtend.DP_TIMES).increment(1);
			
//			// TODO: move gapped extending to reduce side
//			//DP gapped extending
//			dpExtender.setParameters(
//					queryBytes, 0, queryBytes.length,
//					seqbytes, 0, 0+seqSplitLen+overlapLen,
//					filteredExt.get(i).twoHit.getSecondIndexInQuery(),
//					filteredExt.get(i).twoHit
//						.getSecondIndexInDB_Sequence(),
//					gappedXDrop, match, mismatch, gap);
//			try {
//				if(dpExtender.extend(gappedScoreThreshold)){
					//outValueAlign = dpExtender.getAlignment();
//					outValueAlign.setSubjectSeqID(key.toString());
//					outValueAlign.setSubjectSeqStartOffset(
//							seqStartOffset+
//							outValueAlign.getSubjectSeqStartOffset());
//					outKeyAlignScore.set(outValueAlign.getScore());
//					//context.write(outKeyAlignScore, outValueAlign);
//					Text align=new Text();
//					byte[] bytes1 = outValueAlign.getSubjectSeq();
//					align.set(Arrays.toString(bytes1));
//					context.write(key, align);
//					//context.write(outKeyAlignScore, align);
//					context.getCounter(
//							GappedExtend.DP_TIMES).increment(1);
//				}
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
		
		}//~DP extending
	}//~map()
		
	
	/**
	 * 
	 * 
	 * @param querySeq
	 * @param queryBegin
	 * @param queryEnd
	 * @param dbSeq
	 * @param dbBeging
	 * @param dbEnd
	 * @param queryHitIndex
	 * @param dbHitIndex
	 * @param match
	 * @param mistmatch
	 * @param X the allowed distance below the 
	 * best score found for shorter extensions 
	 * @return
	 */
	static GapFreeExtension ungappedExtend(byte[] querySeq, 
			int queryBegin, int queryEnd, 
			byte[] dbSeq, int dbBegin, 
			int dbEnd, TwoHit twoHit,
			int match, int mismatch, int X){
		int queryHitIndex = twoHit.getSecondIndexInQuery();
		int dbHitIndex = twoHit.getSecondIndexInDB_Sequence();
		int forwardBestScore =0;
		int forwardCurrentScore=0;
		//extend forward
		int forwardIndex=0;
		for(; ;forwardIndex++){
			if(forwardIndex+queryHitIndex>=queryEnd 
					|| forwardIndex+dbHitIndex >= dbEnd){
				break;
			}
			//match
			if(querySeq[queryHitIndex+forwardIndex]
					==dbSeq[dbHitIndex+forwardIndex]){
				forwardCurrentScore+=match;
				
			}else{//mismatch
				forwardCurrentScore+=mismatch;
			}
			
			if(forwardCurrentScore>forwardBestScore){
				forwardBestScore = forwardCurrentScore;
			}else{//terminate
				if(forwardBestScore - forwardCurrentScore > X){
					break;
				}
			}
		}
		
		int backwardBestScore=0;
		int backwardCurrentScore=0;
		//extend backward
		int backwardIndex=0;
		for(; ;backwardIndex++){
			if(queryHitIndex-1-backwardIndex<queryBegin
					|| dbHitIndex-1-backwardIndex < dbBegin){
				break;
			}
			
			if(querySeq[queryHitIndex-1-backwardIndex]
					==dbSeq[dbHitIndex-1-backwardIndex]){
				backwardCurrentScore+=match;
			}else{
				backwardCurrentScore+=mismatch;
			}
			
			if(backwardCurrentScore>backwardBestScore){
				backwardBestScore = backwardCurrentScore;
			}else{
				if(backwardBestScore - backwardCurrentScore > X){
					break;
				}
			}
		}
		return new GapFreeExtension(
				forwardBestScore+backwardBestScore,
				twoHit, backwardIndex, forwardIndex);
	}
	


	public static class GapFreeExtension implements 
			Comparable<GapFreeExtension>, Writable{
		public int score;
		public TwoHit twoHit;
		public int leftEnd;
		public int rightEnd;
		public GapFreeExtension(int score, TwoHit hit, 
				int leftEnd, int rightEnd ){
			this.score = score;
			this.twoHit = hit;
			this.leftEnd = leftEnd;
			this.rightEnd = rightEnd;
		}
		@Override
		public int compareTo(GapFreeExtension other) {
			
			return -(this.score - other.score);
		}
		@Override
		public void readFields(DataInput arg0) throws IOException {
			// TODO Auto-generated method stub
			
		}
		@Override
		public void write(DataOutput arg0) throws IOException {
			// TODO Auto-generated method stub
			
		}
	}

}