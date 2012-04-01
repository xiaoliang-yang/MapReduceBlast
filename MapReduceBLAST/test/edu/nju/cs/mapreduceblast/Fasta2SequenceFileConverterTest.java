package edu.nju.cs.mapreduceblast;


import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.util.ReflectionUtils;

public class Fasta2SequenceFileConverterTest {

	public static void main(String[] args) throws Exception{
		// (File fastaInFile, FileSystem fs, Configuration conf, Path fileName)
		Configuration conf = new Configuration();
		FileSystem fs = FileSystem.get(conf);
		String file="nt-small.fa"; //"test.fa";
		File fastaFile = new File(file);
		Path path = new Path(file+".seq");
		Fasta2SequenceFileConverter converter = new Fasta2SequenceFileConverter(fastaFile, fs, conf, path);
		converter.convert2HadoopSequenceFile();
		converter.close();
	}

}
