/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

import java.util.ArrayList;
import java.util.Iterator;

public class State {
	private int identity;
	private State[] Goto=new State[5];
	private State[] Move=new State[5];
	private ArrayList<Word> output=new ArrayList<Word>();
	private State fail=null;
	
	public State(int identity){
		this.identity=identity;
		Goto[Constant.A]=null;
		Goto[Constant.T]=null;
		Goto[Constant.C]=null;
		Goto[Constant.G]=null;
		Goto[Constant.N]=null;
		output.clear();
	}
	
	public int getIdentity(){
		return identity;
	}
	
	public void setGoto(byte character,State targetState){
		Goto[character]=targetState;
	}
	
	public State getGoto(byte character){
		return Goto[character];
	}
	
	public void setMove(byte character,State targetState){
		Move[character]=targetState;
	}
	
	public State getMove(byte character){
		return Move[character];
	}
	
	public void addOutput(Word word){
		output.add(word);
	}
	
	public void addOutput(ArrayList<Word> output){
		Iterator<Word> itr=output.iterator();
		while(itr.hasNext())
			this.output.add(itr.next());
	}
	
	public ArrayList<Word> getOutput(){
		return output;
	}
	
	public ArrayList<Word> getWords(){
		return output;
	} 
	
	public void setFail(State failState){
		fail=failState;
	}
	
	public State getFail(){
		return fail;
	}
	
}
