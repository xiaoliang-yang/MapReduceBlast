/**
 * MapReduce BLAST
 * @author Yang Xiao-liang
 * @email yangxiaoliang2006@gmail.com
 * @copyright June, 2011
 */
package edu.nju.cs.mapreduceblast;

import org.apache.hadoop.util.ProgramDriver;

public class MapReduceBLASTDriver {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int exitCode = -1;
		ProgramDriver pgd = new ProgramDriver();
		try {
			pgd.addClass("convert2seq", Fasta2SequenceFileConverter.class,
							"convert a fasta DNA sequcne file to a Hadoop SequenceFile file");
			pgd.addClass("blastn", BlastnDriver.class, 
					"blastn search on db sequence");
			
			pgd.driver(args);
			exitCode = 0;
		} catch (Throwable e) {
			e.printStackTrace();
		}
		System.exit(exitCode);
	}

}
