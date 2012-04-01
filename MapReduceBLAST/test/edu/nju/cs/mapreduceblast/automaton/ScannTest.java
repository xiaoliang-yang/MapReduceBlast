package edu.nju.cs.mapreduceblast.automaton;

import java.util.ArrayList;

import edu.nju.cs.mapreduceblast.BlastUtils;

public class ScannTest {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		//construct the scanner
		String queryFile = "testquery.fa";
		String query=BlastUtils.readQueryFromLocalFile(queryFile);
		int wordLength=11;
		AutomataMachine automaton=new AutomataMachine(query, wordLength);
		Scanner scanner=new Scanner(automaton);
//		byte[] queryseq = BlastUtils.DNA_String2bytes(queryFile);
		
		String dbFile = "testdb.fa";
		String dbstring = BlastUtils.readQueryFromLocalFile(dbFile);
		byte[] dbseq = BlastUtils.DNA_String2bytes(dbstring);
		
		ArrayList<TwoHit> twoHitList = scanner.scan(dbseq, 0, dbseq.length, 40);
		for(TwoHit twoHit: twoHitList){
			System.out.print("q1:"+twoHit.getFistIndexInQuery()+" d1:"+twoHit.getFistIndexInDB_Sequence()+" ");
			System.out.println("q2:"+twoHit.getSecondIndexInQuery()+" d2:"+twoHit.getSecondIndexInDB_Sequence()+
					" q2-q1 == d2-d1:" + (
					twoHit.getSecondIndexInQuery()-twoHit.getFistIndexInQuery() ==
					twoHit.getSecondIndexInDB_Sequence()-twoHit.getFistIndexInDB_Sequence()) );
		}

	}

}
