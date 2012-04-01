/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

/**
 * @author Administrator
 *
 */
public class Scanner implements Serializable{

	private int[][] moveFunction;
	private OutPut[] outputList;
	private int lengthOfWord;
	
	/**
	 * 
	 * @param automataMachine
	 * @param DB_Suquence
	 */
	public Scanner(AutomataMachine automataMachine){
		moveFunction=automataMachine.getMoveFunction();
		outputList=automataMachine.getOutputList();
		lengthOfWord=automataMachine.getLengthOfWord();
	}
	/**
	 * 
	 * @param firstHit
	 * @param secondHit
	 * @return ArrayList<TwoHit>
	 */
	private ArrayList<TwoHit> outputTwoHit(Hit firstHit,Hit secondHit){
		ArrayList<TwoHit> listOfTwoHit=new 	ArrayList<TwoHit>();;
		int[] indexInquery_A=firstHit.getIndexInQuery();
		int[] indexInquery_B=secondHit.getIndexInQuery();
		
		int indexInDB_Sequence_A=firstHit.getIndexInDB_Sequence();
		int indexInDB_Sequence_B=secondHit.getIndexInDB_Sequence();
		int distance=indexInDB_Sequence_B-indexInDB_Sequence_A;
		if(distance>=lengthOfWord)
		for(int i=0;i<indexInquery_A.length;i++)
			for(int j=0;j<indexInquery_B.length;j++){
				if(indexInquery_B[j]-indexInquery_A[i]==distance){
					listOfTwoHit.add( new TwoHit(indexInDB_Sequence_A,
							indexInDB_Sequence_B,indexInquery_A[i],indexInquery_B[j]));
				}
			}
		
		return listOfTwoHit;
	}
	
	
	private ArrayList<TwoHit> addHit(ArrayList<Hit> listOfHit,Hit hit,int distance_A){
		ArrayList<TwoHit> listOfTwoHit=new 	ArrayList<TwoHit>();
		TwoHit twoHit=null;
		if(listOfHit.isEmpty())
			listOfHit.add(hit);
		else{
			int indexOfCurrentHit=hit.getIndexInDB_Sequence();
			int indexOfPreviousHit;
			Iterator<Hit> itr=listOfHit.iterator();
			Hit previousHit=null;
			while(itr.hasNext()){
				previousHit=itr.next();
				indexOfPreviousHit=previousHit.getIndexInDB_Sequence();
				if(indexOfCurrentHit-indexOfPreviousHit<distance_A){
					merge(listOfTwoHit,outputTwoHit(previousHit,hit));
				}else 
					if(indexOfCurrentHit-indexOfPreviousHit==distance_A){
						merge(listOfTwoHit,outputTwoHit(previousHit,hit));
						itr.remove();
					}else{
							itr.remove();
					}
			}
			listOfHit.add(hit);
				
		}
		return listOfTwoHit;
	}
	
	private void merge(ArrayList<TwoHit> list1,ArrayList<TwoHit> list2){
		Iterator<TwoHit> itr=list2.iterator();
		while(itr.hasNext()){
			list1.add(itr.next());
		}
	}
	
	/**
	 * 
	 * @param DB_Sequence
	 * @param beginIndex
	 * @param endIndex
	 * @param twoHitDistanceA the allowed range for two hits 
	 * (the distance of the indexes of the first bases of each hit ) 
	 * @return
	 */
	public ArrayList<TwoHit> scan(byte[]DB_Sequence,int beginIndex,int endIndex,int twoHitDistanceA){
		ArrayList<TwoHit> listOfTwoHit=new ArrayList<TwoHit>();
		ArrayList<Hit> listOfHit=new ArrayList<Hit>();
		int currentState=0;
		for(int i=beginIndex;i<endIndex;i++){
			currentState=moveFunction[DB_Sequence[i]][currentState];
			ArrayList<Word> words=outputList[currentState].getOutput();
			if(!words.isEmpty()){
				int[]indexInQuery=new int[words.size()];
				for(int j=0;j<words.size();j++)
					indexInQuery[j]=(words.get(j)).getIndex();
				Hit hit=new Hit(i-lengthOfWord+1,indexInQuery);
				merge(listOfTwoHit,addHit(listOfHit,hit,twoHitDistanceA));
			}
		}
		return listOfTwoHit;
	}

	/**
	 * write the Scanner object to a HDFS file
	 * @param path
	 * @throws IOException 
	 */
	public void writeToHDFS_File(FileSystem fs, Path path) throws IOException{
		ObjectOutputStream ostream = new ObjectOutputStream(fs.create(path));
		ostream.writeObject(this);
		ostream.close();
		//TODO:
	}
	
	/**
	 * write Scanner object to file on local disk 
	 * @param file
	 * @throws IOException
	 */
	public void writeToFile(String file) throws IOException{
		FileOutputStream fileOut = new FileOutputStream(file);
		ObjectOutputStream ostream = new ObjectOutputStream(fileOut);
		ostream.writeObject(this);
		ostream.flush();
		ostream.close();
		fileOut.flush();
		fileOut.close();
	}
	
	/**
	 * read a Scanner object from a file on local disk
	 * @param file
	 * @return scanner object
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static Scanner readFromFile(String file) throws IOException, ClassNotFoundException{
		Scanner scanner= null;
		ObjectInputStream ostream = new ObjectInputStream(new FileInputStream(file));
		scanner = (Scanner)ostream.readObject();
		return scanner;
	}
}
