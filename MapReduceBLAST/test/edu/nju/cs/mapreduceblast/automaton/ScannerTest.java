package edu.nju.cs.mapreduceblast.automaton;

import java.util.ArrayList;
import java.util.Iterator;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.util.ReflectionUtils;

import edu.nju.cs.mapreduceblast.BlastUtils;
import edu.nju.cs.mapreduceblast.Fasta2SequenceFileConverter;
import edu.nju.cs.mapreduceblast.automaton.AutomataMachine;
import edu.nju.cs.mapreduceblast.automaton.Scanner;
import edu.nju.cs.mapreduceblast.automaton.TwoHit;


public class ScannerTest { 
	
	public static void main(String [] args) throws Exception{
		//construct the scanner
		String queryFile = "query.fa";
		String query=BlastUtils.readQueryFromLocalFile(queryFile);
		int wordLength=11;
		AutomataMachine automaton=new AutomataMachine(query, wordLength);
		Scanner scanner=new Scanner(automaton);

		//read db sequences and scan two-hits
		Path dbSeqFile = new Path("nt-small.fa.seq");
		Configuration conf = new Configuration();
		FileSystem fs = FileSystem.get(conf);
		SequenceFile.Reader reader = new SequenceFile.Reader(fs, dbSeqFile, conf);
		fs.open(dbSeqFile);
		
		Writable key = (Writable)ReflectionUtils.newInstance(reader.getKeyClass(), conf);
		Writable value = (Writable)ReflectionUtils.newInstance(reader.getValueClass(), conf);
		
		int start_offset ;
		int seqSplitLen ;
		int overlapLen ;
		int counter_bp ;
		byte[] seqbytes = new byte[Fasta2SequenceFileConverter.bpAmountPerLineDefault+
		                           Fasta2SequenceFileConverter.bpAmountOverlapDefault];
		while(reader.next(key, value)){
			byte[] bytes = ((BytesWritable)value).getBytes();
			start_offset = BlastUtils.bytes2int(bytes, 0);
			seqSplitLen = BlastUtils.bytes2int(bytes, 4);
			overlapLen = BlastUtils.bytes2int(bytes, 8);
			
			counter_bp = 0;
			for(int i=12; i< bytes.length; i++){
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
			ArrayList<TwoHit> twoHitList = scanner.scan(seqbytes, 0, seqSplitLen+overlapLen, 10);
			Iterator<TwoHit> itr=twoHitList.iterator();
			while(itr.hasNext()){
				System.out.println(itr.next().toString());
			}
		}
		
		
		
	}
}  


