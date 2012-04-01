package edu.nju.cs.mapreduceblast;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import org.apache.hadoop.io.Writable;

public class Alignment implements Writable{
	
	private byte[] querySeq;
	private byte[] subjectSeq;
	private int lenOfAlignment;
	
	private int score;
	
	private String subjectSeqID;

	// start offset in query sequence
	private int queryStartOffset;
	private int queryEndOffset;
	//
	private int subjectSeqStartOffset;
	private int subjectSeqEndOffset;
	
	
	
	public String getSubjectSeqID() {
		return subjectSeqID;
	}
	public void setSubjectSeqID(String subjectSeqID) {
		this.subjectSeqID = subjectSeqID;
	}
	public int getQueryStartOffset() {
		return queryStartOffset;
	}
	public void setQueryStartOffset(int queryStartOffset) {
		this.queryStartOffset = queryStartOffset;
	}
	public int getQueryEndOffset() {
		return queryEndOffset;
	}
	public void setQueryEndOffset(int queryEndOffset) {
		this.queryEndOffset = queryEndOffset;
	}
	public int getSubjectSeqStartOffset() {
		return subjectSeqStartOffset;
	}
	public void setSubjectSeqStartOffset(int subjectSeqStartOffset) {
		this.subjectSeqStartOffset = subjectSeqStartOffset;
	}
	public int getSubjectSeqEndOffset() {
		return subjectSeqEndOffset;
	}
	public void setSubjectSeqEndOffset(int subjectSeqEndOffset) {
		this.subjectSeqEndOffset = subjectSeqEndOffset;
	}
	public byte[] getQuerySeq() {
		return querySeq;
	}
	public void setQuerySeq(byte[] querySeq) {
		this.querySeq = querySeq;
	}
	public byte[] getSubjectSeq() {
		return subjectSeq;
	}
	public void setSubjectSeq(byte[] subjectSeq) {
		this.subjectSeq = subjectSeq;
	}
	public int getScore() {
		return score;
	}
	public void setScore(int score) {
		this.score = score;
	}
	public void setFields(byte[] query, byte[] subject, int score, 
			int queryStartOffset, int queryEndOffset,
			int subjectSeqStartOffset, int subjectSeqEndOffset) throws Exception{
		if(query.length != subject.length){
			throw new Exception("query length and subject lenght are not equal!");
		}
		this.querySeq=query;
		this.subjectSeq=subject;
		this.score=score;
		this.queryEndOffset = queryEndOffset;
		this.queryStartOffset = queryStartOffset;
		this.subjectSeqStartOffset = subjectSeqStartOffset;
		this.subjectSeqEndOffset = subjectSeqEndOffset;
		this.lenOfAlignment = query.length;
	}
	
	public Alignment(byte[] query, byte[] subject, int score, 
			int queryStartOffset, int queryEndOffset,
			int subjectSeqStartOffset, int subjectSeqEndOffset) throws Exception{
		
		setFields(query, subject, score, queryStartOffset, queryEndOffset,
				subjectSeqStartOffset, subjectSeqEndOffset);
	}
	
	public Alignment(){
		
	}
	
	@Override
	public void readFields(DataInput in) throws IOException {
		score = in.readInt();
		queryStartOffset = in.readInt();
		queryEndOffset = in.readInt();
		subjectSeqStartOffset = in.readInt();
		subjectSeqEndOffset = in.readInt();
		
		lenOfAlignment = in.readInt();
		querySeq = new byte[lenOfAlignment];
		in.readFully(querySeq, 0, lenOfAlignment);
		subjectSeq = new byte[lenOfAlignment];
		in.readFully(subjectSeq, 0, lenOfAlignment);
	}
	
	@Override
	public void write(DataOutput out) throws IOException {
		out.writeInt(score);
		out.writeInt(queryStartOffset);
		out.writeInt(queryEndOffset);
		out.writeInt(subjectSeqStartOffset);
		out.writeInt(subjectSeqEndOffset);
		
		out.writeInt(lenOfAlignment);
		out.write(querySeq);
		out.write(subjectSeq);
		
	}
	
}
