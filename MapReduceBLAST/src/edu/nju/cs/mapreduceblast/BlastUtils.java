/**
 * MapReduce BLAST
 * @author Yang Xiao-liang
 * @email yangxiaoliang2006@gmail.com
 * @copyright June, 2011
 */
package edu.nju.cs.mapreduceblast;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

public class BlastUtils {
	/**
	 * encode one bp in 4 bits 
	 * DNA >>
	 * 1 bp represented by 4 bits,
	 * A: 0000 
	 * T: 0001 (U in RNA equals to T in DNA)
	 * C: 0010 
	 * G: 0011 
	 * N: 0100 (unknown base, include N, I, V, D, B, H, Y, M, K, R, W, S)
	 * @param bp
	 * @param leftmost4bits
	 *            if true, the leftmost 4 bits represents the bp (ATCGN) else,
	 *            the rightmost 4bits
	 * @return
	 * @throws Exception
	 */
	public static byte encode_bp(char bp, boolean leftHalf)
			throws Exception {
		byte code;
		if (leftHalf) {
			switch (bp) {
			case 'A':	case 'a':
				code = 0x00;
				break;
			case 'T':	case 't':
			case 'U':	case 'u':
				code = 0x10;
				break;
			case 'C':	case 'c':
				code = 0x20;
				break;
			case 'G':	case 'g':
				code = 0x30;
				break;
			case 'N': 	case 'n':	
			case 'I':	case 'i':
			case 'V':	case 'v':	
			case 'D':	case 'd':
			case 'B':	case 'b':
			case 'H':	case 'h':
			case 'W':	case 'w':
			case 'S':	case 's':
			case 'K':	case 'k':
			case 'M':	case 'm':
			case 'Y':	case 'y':
			case 'R':	case 'r':
				code = 0x40;
				break;
			default:
				throw new Exception("Invalid input DNA base: " + bp);
			}
		} else {
			switch (bp) {
			case 'A':	case 'a':
				code = 0x00;
				break;
			case 'T':	case 't':
			case 'U':	case 'u':
				code = 0x01;
				break;
			case 'C':	case 'c':
				code = 0x02;
				break;
			case 'G':	case 'g':
				code = 0x03;
				break;
			case 'N': 	case 'n':
			case 'I':	case 'i':
			case 'V':	case 'v':
			case 'D':	case 'd':
			case 'B':	case 'b':
			case 'H':	case 'h':
			case 'W':	case 'w':
			case 'S':	case 's':
			case 'K':	case 'k':
			case 'M':	case 'm':
			case 'Y':	case 'y':
			case 'R':	case 'r':
				code = 0x04;
				break;
			default:
				throw new Exception("Invalid input DNA base: " + bp);
			}
		}
		return code;
	}

	/**
	 * decode half byte to DNA base
	 * * DNA >>
	 * 1 bp represented by 4 bits,
	 * A: 0000 T: 0001
	 * C: 0010 G: 0011 
	 * N: 0100 (unknown base)
	 * @param code
	 * @param leftHalf if true, real code is on the left half byte, false else 
	 * @return DNA base (byte)'A'
	 * @throws Exception 
	 */
	public static byte decode_bp(byte code, boolean leftHalf) throws Exception{
		byte bp;
		if(leftHalf){
			code >>>= 4;
		}else{
			code &= 0x0F;
		}
		switch(code){
		case (byte)0x00:
			bp = (byte)'A'; break;
		case (byte)0x01:
			bp = (byte)'T'; break;
		case (byte)0x02:
			bp = (byte)'C'; break;
		case (byte)0x03:
			bp = (byte)'G'; break;
		case (byte)0x04:
			bp = (byte)'N'; break;
		default:
			throw new Exception("Error: bad DNA base code :"+code);
		}
		
		return bp;
	}
	
	/**
	 * big-endian alignment
	 * @param value
	 * @return
	 */
	public static byte[] int2bytes(int value){
		byte[] bytes = new byte[4];
		
		for(int i=0; i<4; i++){
			bytes[3-i] = (byte)(value>>>i*8);
		}
		return bytes;
	}
	

	/**
	 * convert a 4byte-array to an integer
	 * big-endian alignment
	 * @param bytes 
	 * @param begin the starting index of the integer
	 */
	public static int bytes2int(byte[] bytes, int begin){
		int value = 0;
		int mask = 0x000000FF;
		for(int i= 0; i< 4; i++){
			value <<= 8;
			value |= mask&bytes[begin+i];
			//@debug
			//System.out.printf("i=%d: val: %08X\n",i,value);
		}
		return value;
	}
	/**
	 * get the sequence ID from the sequence header line 
	 * the blast sequence header is like this:
	 * >gi|224589800|ref|NC_000001.10| Homo sapiens chromosome 1, GRCh37 primary reference assembly
	 * @param seqHeader the header line
	 * the default pattern is ">gi|.+|" 
	 * @return the sequence id
	 */
	public static String getSeqID(String seqHeader) throws Exception {
		String reg=">gi\\|[^\\|]+?\\|";
		Pattern pattern = Pattern.compile(reg);
		Matcher matcher = pattern.matcher(seqHeader);
		if(matcher.find()){
			int start=matcher.start();
			int end = matcher.end();
			if (end-start < ">gi||".length()){
				throw new Exception("bad sequence header line");
			}
			return seqHeader.substring(start+4, end-1);
		}else{
			throw new Exception("regex matches failed, no \""+reg+"\"found");
		}
		
	}
	
	/**
	 * read one DNA sequecne from local file,
	 * a sequecnce header line is begin with indentifier '>',
	 * if multiple sequences are in the file, only the first is read
	 * @param file
	 * @return 
	 * @throws IOException 
	 */
	public static String readQueryFromLocalFile(String file) throws IOException{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		StringBuilder builder = new StringBuilder();
		String line ;
		boolean isFirstSeq = true;
		// only the lines of the first sequence is read
		while((line=reader.readLine())!=null){
			if(line.length()==0){ continue;}
			//skip the header line
			if(line.charAt(0)=='>'){
				if(isFirstSeq){
					continue;
				}else{
					break;
				}
			}
			isFirstSeq = false;
			builder.append(line);
		}
		reader.close();
		return builder.toString();
	}
	/**
	 * convert DNA string to byte array, each byte represent a base, 
	 * encoded by BlastUtils.encode_bp()
	 * @param DNA_Str
	 * @return byte array
	 * @throws Exception 
	 */
	public static byte[] DNA_String2bytes(String DNA_Str) throws Exception{
		byte[] bytes = new byte[DNA_Str.length()];
		for(int i=0; i< DNA_Str.length(); i++){
			bytes[i] = BlastUtils.encode_bp(DNA_Str.charAt(i), false);
		}
		return bytes;
	}
	
	public static void writeBytesToHDFS(FileSystem fs, Configuration conf, Path filePath, byte[] bytes) throws IOException{
		FSDataOutputStream out = fs.create(filePath);
		out.write(bytes);
		out.close();
	}
	
	/**
	 * caution: this method can only deal with small files
	 * @throws IOException 
	 */
	public static byte[] readBytesFromLoalFile(String file) throws IOException{
		File bytesFile = new File(file);
		byte[] bytes = new byte[(int)bytesFile.length()];
		FileInputStream in = new FileInputStream(bytesFile);
		in.read(bytes);
		in.close();
		return bytes;
	}
	

}
