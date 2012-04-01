/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

import java.io.Serializable;
import java.util.ArrayList;

public class OutPut implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -7934438831311225279L;
	private ArrayList<Word> output=new ArrayList<Word>();
	
	public OutPut(ArrayList<Word> output){
		this.output=output;
	}
	
	public ArrayList<Word> getOutput(){
		return output;
	}
}
