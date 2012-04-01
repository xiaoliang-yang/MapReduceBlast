package edu.nju.cs.mapreduceblast;


import java.io.IOException;
import java.util.Arrays;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import edu.nju.cs.mapreduceblast.automaton.Scanner;

public class BlastnReducer extends Reducer<Text, BytesWritable,
		Text,BytesWritable>{
	
	//blastn arguments
	private int match;
	private int mismatch;
	private int gap;
	
	
	private int gappedScoreThreshold;

	private int gappedXDrop;
	
	//fields for the alignment
	private DPExtender dpExtender;
	
	private byte[] queryBytes;
	/**
	 * load the scanner and query bytes file 
	 */
	@Override
	public void setup(Context context) throws IOException {
		
		//get blastn arguments/////////////////
		Configuration conf = context.getConfiguration();
		this.match = conf.getInt("blastn.argument.wordLen", 11);
		this.match= conf.getInt("blastn.argument.matchScore", 2);
		this.mismatch= conf.getInt("blastn.argument.mistmatchScore",
				-3);
		this.gap = conf.getInt("blastn.argument.gapPenalty", -5);
		
		this.gappedScoreThreshold = conf.getInt(
				"blastn.argument.gappedScoreThreshold", 200*match);
		
		this.gappedXDrop = conf.getInt("balstn.argument.gappedX",
				10*match);
		
		// load query bytes from cache file 
		String queryBytesFile = conf.get("blastn.query.byte.file");
		queryBytes = BlastUtils.readBytesFromLoalFile(queryBytesFile);
		
		dpExtender = new DPExtender();
				
	}
	
	
	// emit key=(SeqID,startOffset,XBaseIndex,YBaseIndex)
	// value=(seqbytes)
	@Override
	public void reduce(Text key, Iterable<BytesWritable> values, 
			Context context)throws IOException, InterruptedException{
		
		byte[] seqBytes = values.iterator().next().getBytes();
		String[] params = key.toString().split(",");
		int seqID = Integer.valueOf(params[0]);
		int seqStartOffset = Integer.valueOf(params[1]);
		int XBaseIndex = Integer.valueOf(params[2]);
		int YBaseIndex = Integer.valueOf(params[3]);
		
		Alignment outValueAlign ;
		
		//DP gapped extending
		try {
			
			dpExtender.setParameters(
					queryBytes, 0, queryBytes.length,
					seqBytes, 0, seqBytes.length, 
					XBaseIndex, YBaseIndex, gappedXDrop, 
					match, mismatch, gap);
			
			if(dpExtender.extend(gappedScoreThreshold)){
				outValueAlign = dpExtender.getAlignment();
				
				//TODO: output format
				BytesWritable bytesSegment = new BytesWritable();
				byte[] dbSegment = outValueAlign.getSubjectSeq();
				bytesSegment.set(dbSegment, 0, dbSegment.length);
				
				context.write(key, bytesSegment);
				//context.write(outKeyAlignScore, outValueAlign);
				//context.write(outKeyAlignScore, align);
				context.getCounter(
						BlastnMapper.GappedExtend.DP_TIMES).increment(1);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
