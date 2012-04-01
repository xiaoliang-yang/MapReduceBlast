package edu.nju.cs.mapreduceblast.automaton;

import java.util.ArrayList;
import java.util.Iterator;

public class TestByte {

	public static void main(String [] args){
		
		ArrayList<String> list=new ArrayList<String>();
		list.add("liu");
		list.add("yu");
		list.add("long");
		
		Iterator<String> itr=list.iterator();
		
		while(itr.hasNext()){
			System.out.println(itr.next());
			itr.remove();
		}
		
		for(int i=0;i<list.size();i++)
			System.out.println(list.get(i));
	}
}
