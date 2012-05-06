/**
 * MapReduce BLAST
 * @author Yang Xiao-liang
 * @email yangxiaoliang2006@gmail.com
 * @copyright June, 2011,
 */
package edu.nju.cs.mapreduceblast;

import java.io.File;
import java.net.URI;
import java.util.Date;
import java.util.HashMap;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.filecache.DistributedCache;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.SequenceFileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

import edu.nju.cs.mapreduceblast.automaton.AutomataMachine;
import edu.nju.cs.mapreduceblast.automaton.Scanner;

public class BlastnDriver {
	static class Argument {
		public Argument(String name, String value, boolean isNecessary, boolean isSet, String description){
			this.name = name;
			this.value = value;
			this.isNecessary = isNecessary;
			this.isSet = isSet;
			this.description = description;
		}
		public String name;
		public String value;
		public boolean isNecessary = false;
		public boolean isSet = false;
		public String description;
	}
	public static void printUsage(){
		System.out.println("<-query <local query_file_name>> <-db <database file>> <-r reduceNum> <-out <out path>>" +
		" [-UT <ungapped threshold (int)>] [-w <word length (int)>] [- d <distance (int)>]");
	}
	public static void main(String[] args) throws Exception {
		if(args.length==0){
			printUsage();
			System.exit(-1);
		}
		HashMap<String, Argument> arguments = new HashMap<String, Argument>();
		arguments.put("-query", new Argument("-query", null, true, false, " query file path"));
		arguments.put("-db", new Argument("-db", null, true, false, "sequence database path"));
		arguments.put("-out", new Argument("-out", null, true, false, "out put file path"));
		arguments.put("-r", new Argument("-r", null, false, false, "reduce num"));
		arguments.put("-w", new Argument("-w", null, false, false, "word length"));
		arguments.put("-d", new Argument("-d", null, false, false, "two hits distance A"));
		arguments.put("-UT", new Argument("-UT", null, false, false, "ungapped extend threshold score"));
		//TODO: match, mismatch, gap costs
		
		//set arguments values
		Argument arg;
		for(int i=0; i<args.length; i++){
			arg = arguments.get(args[i]);
			if(arg != null){
				if(i+1<args.length){
					arg.value = args[i+1];
					arg.isSet = true;
				}else{
					System.out.println(arg.name+" incorrect");
					printUsage();
					System.exit(-1);
				}
			}
		}
		
		//get values
		String queryFile=null;
		arg = arguments.get("-query");
		if(arg.isSet){
			queryFile = arg.value;
		}else{
			System.out.println("input query is not set");
			printUsage();
			System.exit(-1);
		}

		String dbFile=null;
		arg = arguments.get("-db");
		if(arg.isSet){
			dbFile= arg.value;
		}else{
			System.out.println("db is not set");
			printUsage();
			System.exit(-1);
		}
		
		String outFile=null;
		arg = arguments.get("-out");
		if(arg.isSet){
			outFile= arg.value;
		}else{
			System.out.println("out file is not set");
			printUsage();
			System.exit(-1);
		}
		
		//TODO: parse the args, and set balstn arguments in configuration
		Configuration conf = new Configuration();
		FileSystem fs = FileSystem.get(conf);
		
		arg = arguments.get("-UT");
		if(arg.isSet){
			int ungappedThreshold = Integer.parseInt(arg.value);
			conf.setInt("blastn.argument.ungappedScoreThreshold", ungappedThreshold);
		}
		arg = arguments.get("-d");
		if(arg.isSet){
			int twoHitDistance = Integer.parseInt(arg.value);
			conf.setInt("blastn.argument.twoHitDistanceA",twoHitDistance);
		}
		
		int wordLength=11;
		arg = arguments.get("-w");
		if(arg.isSet){
			wordLength = Integer.parseInt(arg.value);
			conf.setInt("blastn.argument.wordLen", wordLength);
		}
		
		int reduceNum = 1;
		arg = arguments.get("-r");
		if(arg.isSet){
			reduceNum = Integer.parseInt(arg.value);
		}
		////////////////////////////////////////////////////////////////////////////////
		System.out.println(new Date());
		//construct the scanner
		String query = BlastUtils.readQueryFromLocalFile(queryFile);
		AutomataMachine automaton=new AutomataMachine(query, wordLength);
		Scanner scanner=new Scanner(automaton);
		System.out.println("Scanner constructed.");
		
		//and write scanner object  to file
		String scannerObjFileName = "_"+(new File(queryFile)).getName()+"-Scanner.obj";
		Path scannerObjPath = new Path(scannerObjFileName); 
		scanner.writeToHDFS_File(fs, scannerObjPath);
		System.out.println("Scanner object wrriten to file: "+scannerObjFileName);
		
		//distribute the scanner object file to slave nodes by using DistributedCache
		DistributedCache.addCacheFile(new URI(scannerObjFileName+"#"+scannerObjFileName), conf);
		conf.set("blastn.scanner.file", scannerObjFileName);
		System.out.println("scanner object distributed to local disks of slave nodes");
		
		//distribute query bytes to slave nodes
		byte[] queryBytes = BlastUtils.DNA_String2bytes(query);
		String queryBytesFile = "_query.byte";
		System.out.print("writing query bytes to HDFS...");
		Path queryBytesFilePath = new Path(queryBytesFile);
		BlastUtils.writeBytesToHDFS(fs, conf, queryBytesFilePath, queryBytes);
		System.out.println("...Done!");
		DistributedCache.addCacheFile(new URI(queryBytesFile+"#"+queryBytesFile), conf);
		conf.set("blastn.query.byte.file",queryBytesFile);
		
		//create symbolic link to cache files
		DistributedCache.createSymlink(conf);
		
		//balstn mapreduce job
		Job job = new Job(conf);
		job.setJarByClass(BlastnDriver.class);
		job.setJobName("blastn");
		job.setMapperClass(BlastnMapper.class);
		job.setNumReduceTasks(reduceNum);
		
		job.setOutputKeyClass(Text.class);
		job.setOutputValueClass(BytesWritable.class);
		
		job.setInputFormatClass(SequenceFileInputFormat.class);
//		job.setOutputFormatClass(SequenceFileOutputFormat.class);
		
		FileInputFormat.addInputPath(job, new Path(dbFile));
		FileOutputFormat.setOutputPath(job, new Path(outFile));
		System.out.println(new Date().toString()+"job prepared!");
		job.waitForCompletion(true);

		//clear
		fs.deleteOnExit(scannerObjPath);
		fs.deleteOnExit(queryBytesFilePath);
	}

}
