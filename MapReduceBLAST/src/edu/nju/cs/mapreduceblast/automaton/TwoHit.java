/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

public class TwoHit {

	int firstIndexInDB_Sequence;
	int secondIndexInDB_Sequence;
	int firstIndexInQuery;
	int secondIndexInQuery;
	
	public TwoHit(int firstIndexInDB_Sequence,int secondIndexInDB_Sequence,
			int firstIndexInQuery,int secondIndexInQuery){
		this.firstIndexInDB_Sequence=firstIndexInDB_Sequence;
		this.secondIndexInDB_Sequence=secondIndexInDB_Sequence;
		this.firstIndexInQuery=firstIndexInQuery;
		this.secondIndexInQuery=secondIndexInQuery;
	}
	
	public int getFistIndexInDB_Sequence(){
		return firstIndexInDB_Sequence;
	}
	
	public int getSecondIndexInDB_Sequence(){
		return secondIndexInDB_Sequence;
	}
	
	public int getFistIndexInQuery(){
		return firstIndexInQuery;
	}
	
	public int getSecondIndexInQuery(){
		return secondIndexInQuery;
	}
	
	public String toString(){
		return firstIndexInDB_Sequence+" "+secondIndexInDB_Sequence+";"+
		firstIndexInQuery+" "+secondIndexInQuery;
	}
}
