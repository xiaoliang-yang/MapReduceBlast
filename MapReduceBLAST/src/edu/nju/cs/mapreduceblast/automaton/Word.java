/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

import java.io.Serializable;

public class Word implements Serializable{
	private int index;
	private byte[] content;
	
	public Word(byte[] content,int index){
		this.content=content;
		this.index=index;
	}
	
	public byte[] getContent(){
		return content;
	}
	
	public int getIndex(){
		return index;
	}
}
