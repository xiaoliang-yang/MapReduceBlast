/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

/**
 * @author Administrator
 *
 */
public class Hit {

	private int indexInDB_Sequence;
	private int[] indexInQuery;
	
	public Hit(int indexInDB_Sequence,int[] indexInQuery){
		this.indexInDB_Sequence=indexInDB_Sequence;
		this.indexInQuery=indexInQuery;
	}
	
	public int getIndexInDB_Sequence(){
		return indexInDB_Sequence;
	}
	
	public int[] getIndexInQuery(){
		return indexInQuery;
	}
}
