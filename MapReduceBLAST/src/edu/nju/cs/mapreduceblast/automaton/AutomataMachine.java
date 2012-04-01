/**
 * @author Liu Yu-long
 * @date	June, 2011
 */
package edu.nju.cs.mapreduceblast.automaton;

import java.util.LinkedList;
import java.util.Queue;

import edu.nju.cs.mapreduceblast.BlastUtils;

public class AutomataMachine {
	
	private State startStateOfAutomataMachine;
	private int lengthOfWord;
	private int stateCount=0;
	
	public AutomataMachine(String query,int lengthOfWord) throws Exception{
		this.lengthOfWord=lengthOfWord;
		startStateOfAutomataMachine=constructionOfAutomataMachine(query,lengthOfWord);
	}
	
	public State getStartStateOfAutomataMachine(){
		return startStateOfAutomataMachine;
	}
	
	public int getStateCount(){
		return stateCount;
	}
	
	public int getLengthOfWord(){
		return lengthOfWord;
	}
	
	public int[][] getMoveFunction(){
		int [][] moveFunction=new int[5][stateCount+1];
		Queue<State> queue = new LinkedList<State>();
		State currentState,targetState,tempState;
		int currentIdentity,targetIdentity;
		
		currentIdentity=startStateOfAutomataMachine.getIdentity();
		for (byte character : Constant.alphabet) { 
			targetState=startStateOfAutomataMachine.getGoto(character);
			targetIdentity=targetState.getIdentity();
			moveFunction[character][currentIdentity]=targetIdentity;
			if(targetState!=startStateOfAutomataMachine)
				queue.offer(targetState);
         } 
		
		while((currentState=queue.poll())!=null){
			currentIdentity=currentState.getIdentity();
			for (byte character : Constant.alphabet) { 
				targetState=currentState.getGoto(character);
				if(targetState!=null){
					queue.offer(targetState);
					targetIdentity=targetState.getIdentity();
					moveFunction[character][currentIdentity]=targetIdentity;
				}else{
					tempState=currentState.getFail();
					targetState=tempState.getMove(character);
					targetIdentity=targetState.getIdentity();
					moveFunction[character][currentIdentity]=targetIdentity;
				}
	         } 
		}
		return moveFunction;
	}
	
	public OutPut[] getOutputList(){
		OutPut[] outputList=new OutPut[stateCount+1];
		Queue<State> queue = new LinkedList<State>();
		State currentState,targetState;
		int currentIdentity;
		
		currentIdentity=startStateOfAutomataMachine.getIdentity();
		outputList[currentIdentity]=new OutPut(startStateOfAutomataMachine.getOutput());
		for (byte character : Constant.alphabet) { 
			targetState=startStateOfAutomataMachine.getGoto(character);
			if(targetState!=startStateOfAutomataMachine)
				queue.offer(targetState);
         } 
		
		while((currentState=queue.poll())!=null){
			currentIdentity=currentState.getIdentity();
			outputList[currentIdentity]=new OutPut(currentState.getOutput());
			for (byte character : Constant.alphabet) { 
				targetState=currentState.getGoto(character);
				if(targetState!=null){
					queue.offer(targetState);
				}
	         } 
		}
		
		return outputList;
	}
	
	
	private byte[] string2bytes(String s) throws Exception{
		char[] temp=s.toCharArray();
		byte[] bytes=new byte[temp.length];
		
		for(int i=0;i<temp.length;i++){
			bytes[i]=BlastUtils.encode_bp(temp[i],false);
		}
		return bytes;
	}
	
	private Word[] cutSequenceIntoWords(String sequence,int wordLength) throws Exception{
		int sequenceLength=sequence.length();
		
		if(sequenceLength<wordLength){
			System.out.println("The sequence can't be cut into words!");
			return null;
		}
		
		int numOfWords=sequenceLength-wordLength+1;
		Word[] wordList=new Word[numOfWords];
		int beginIndex=0,endIndex=wordLength;
		
		for(int i=0;i<numOfWords;i++){
			String segment=sequence.substring(beginIndex, endIndex);
			byte[] wordContent=string2bytes(segment);
			wordList[i]=new Word(wordContent,beginIndex);
			beginIndex++;
			endIndex++;
		}
		
		return wordList;
	}
	
	private  void enter(State startState,Word word){
		byte[] content=word.getContent();
		int lengthOfWord=content.length;
		
		State currentState=startState;
		int index=0;
//		���������޸ı�֤���鲻Խ�磬������ͬword��content��ͬ��ʱ����������Խ������
		while(index<lengthOfWord&&currentState.getGoto(content[index])!=null){
			currentState=currentState.getGoto(content[index]);
			index++;
		}
		
		for(;index<lengthOfWord;index++){
			stateCount++;
			State newState=new State(stateCount);
			currentState.setGoto(content[index], newState);
			currentState=newState;
		}
		
		currentState.addOutput(word);
		
	}
	
	private  State constructionOfAutomataMachine(String query,int lengthOfWord) throws Exception{
		
		Word[] wordList=cutSequenceIntoWords(query,lengthOfWord);
		/*
		 * ����Goto��output/
		 */
		State startState=new State(0);
		
		for(int i=0;i<wordList.length;i++)
			enter(startState,wordList[i]);
		
		for (byte character : Constant.alphabet) { 
			if(startState.getGoto(character)==null)
				startState.setGoto(character, startState);
         } 
		
		/*
		 *����fail��output/
		 */
		Queue<State> queue = new LinkedList<State>();  
		State currentState,targetState,tempState;
		
		for (byte character : Constant.alphabet) { 
			targetState=startState.getGoto(character);
			if(targetState!=startState){
				queue.offer(targetState);
				targetState.setFail(startState);
			}
         } 

			
			while((currentState=queue.poll())!=null){  
				for (byte character : Constant.alphabet) { 
					targetState=currentState.getGoto(character);
					if(targetState!=null){
						queue.offer(targetState);
						tempState=currentState.getFail();
						while(tempState.getGoto(character)==null){
							tempState=tempState.getFail();
						}
						State failState=tempState.getGoto(character);
						targetState.setFail(failState);
						targetState.addOutput(failState.getOutput());
					}
		         } 
			}
			
		/*
		 *�������յ�ȷ������״̬�Զ���ȥ����failure transitions
		 *�����move function/
		 */	
		queue.clear();
		
		for (byte character : Constant.alphabet) { 
			targetState=startState.getGoto(character);
			startState.setMove(character, targetState);
			if(targetState!=startState)
				queue.offer(targetState);
         } 
		
		while((currentState=queue.poll())!=null){
			for (byte character : Constant.alphabet) { 
				targetState=currentState.getGoto(character);
				if(targetState!=null){
					queue.offer(targetState);
					currentState.setMove(character, targetState);
				}else{
					tempState=currentState.getFail();
					targetState=tempState.getMove(character);
					currentState.setMove(character, targetState);
				}
	         } 
		}
		return startState;
	}
}
