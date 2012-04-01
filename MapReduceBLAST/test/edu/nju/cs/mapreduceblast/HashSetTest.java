package edu.nju.cs.mapreduceblast;

import java.util.HashSet;

public class HashSetTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		HashSet<String> set = new HashSet<String>();
		String str = 13+","+884;
		set.add(str);
		String str2 = 13+","+884;
		System.out.println("set contains"+str2+":"+set.contains(str2));

	}

}
