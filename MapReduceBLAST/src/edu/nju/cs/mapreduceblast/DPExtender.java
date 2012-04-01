/**
 * @author Yang Xiao-liang
 * @email yangxiaoliang2006@gmail.com
 * @copyright June, 2011
 */
package edu.nju.cs.mapreduceblast;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

class Area {
	//coordinates for a rectangle area
	public int X1;
	public int X2;
	public int Y1;
	public int Y2;
	//the optimal score
	public int score;

	public Area(int x1, int x2, int y1, int y2, int score) {
		this.X1 = x1;
		this.X2 = x2;
		this.Y1 = y1;
		this.Y2 = y2;
		this.score = score;
	}
	public String toSting(){
		return "<"+X1+","+X2+","+Y1+","+Y2+"> "+score;
	}
}

/**
 * modified 2011-6-28 
 * by yang xiaoliang
 * extend in two directions
 * trace back path in two directions
 * split the big method dbExtend() into two (the other is divideIntoTwoAreas())
 * 
 * modified 2011-6-27
 * by yang xiaoliang
 * bounded trace back path
 * write points to file for debugging
 * 
 * @author yang xiaoliang
 *
 */

public class DPExtender {
	
	// for recording the X-drop explored strip's boundary (for the Y-coordinate)
	private HashMap<Integer, Integer> boundaryYUpper = new HashMap<Integer, Integer>();
	private HashMap<Integer, Integer> boundaryYLower = new HashMap<Integer, Integer>();
	// parameters
	private byte[] XStr;
	private byte[] YStr;
	private int XBegin;
	private int XEnd;
	private int YBegin;
	private int YEnd;
	private int XBaseIndex;
	private int YBaseIndex;
	private int match;
	private int mismatch;
	private int gap;
	private int drop;
	
	// score
	private int leftScore;
	private int rightScore;
	// relative coordinates
	private int leftScoreXInMatrix;
	private int leftScoreYInMatrix;
	private int rightScoreXInMatrix;
	private int rightScoreYInMatrix;
	
	// record critical path points in the matrix
	// the coordinates may be positive or negative values 
	private HashSet<String> criticalPathRelative = new HashSet<String>();
	
	//@debug
	private boolean DEBUG=false;
	private static boolean DEBUG_WRITE_PATH = false;
	private boolean DEBUG_VERBOSE=false;
	
	private PrintWriter tracebackWriter=null;
	private PrintWriter criticalPathWriter = null;
	
	/**
	 * note: it is supposed that the follow condition must be satisfied: 
	 * XStr[XBaseIndex]==YStr[YBaseIndex]
	 * 
	 * @param XStr
	 * @param XBegin
	 * @param XEnd
	 * @param YStr
	 * @param YBegin
	 * @param YEnd
	 * @param XBaseIndex
	 * @param YBaseIndex
	 * @param drop
	 * @param match
	 * @param mismatch
	 * @param gap
	 * @throws IOException
	 */
	public DPExtender(
			byte[] XStr, int XBegin, int XEnd, 
			byte[] YStr, int YBegin, int YEnd,
			int XBaseIndex, int YBaseIndex,
			int drop, int match, int mismatch,
			int gap) throws IOException{
		this.XStr = XStr;
		this.XBegin = XBegin;
		this.XEnd = XEnd;
		this.YStr = YStr;
		this.YBegin = YBegin;
		this.YEnd = YEnd;
		this.XBaseIndex = XBaseIndex;
		this.YBaseIndex = YBaseIndex;
		this.match = match;
		this.mismatch = mismatch;
		this.gap = gap;
		this.drop = drop;
		
		// workspace for store of optimal values
		if(DEBUG_WRITE_PATH){
			criticalPathWriter = new PrintWriter("pointsCriticalPath.txt");
			tracebackWriter = new PrintWriter("pointsBandTracePath.txt");
		}
	}
	
	public DPExtender() throws IOException{
				// workspace for store of optimal values
		if(DEBUG_WRITE_PATH){
			criticalPathWriter = new PrintWriter("pointsCriticalPath.txt");
			tracebackWriter = new PrintWriter("pointsBandTracePath.txt");
		}
	}
	
	public void clear(){
		leftScore=0;
		rightScore=0;
		leftScoreXInMatrix=0;
		leftScoreYInMatrix=0;
		rightScoreXInMatrix=0;
		rightScoreYInMatrix=0;
		
		criticalPathRelative.clear();
		boundaryYLower.clear();
		boundaryYUpper.clear();
	}
	
	/**
	 * @param XStr the query sequence
	 * @param XBegin 
	 * @param XEnd
	 * @param YStr the subject sequence
	 * @param YBegin
	 * @param YEnd
	 * @param XBaseIndex the hit index in query
	 * @param YBaseIndex the hit index in subject
	 * @param drop X-drop Threshold
	 * @param match
	 * @param mismatch
	 * @param gap
	 * @throws IOException
	 */
	public void setParameters(
			byte[] XStr, int XBegin, int XEnd, 
			byte[] YStr, int YBegin, int YEnd,
			int XBaseIndex, int YBaseIndex,
			int drop, int match, int mismatch,
			int gap) throws IOException{
		// set the parameters
		this.XStr = XStr;
		this.XBegin = XBegin;
		this.XEnd = XEnd;
		this.YStr = YStr;
		this.YBegin = YBegin;
		this.YEnd = YEnd;
		this.XBaseIndex = XBaseIndex;
		this.YBaseIndex = YBaseIndex;
		this.match = match;
		this.mismatch = mismatch;
		this.gap = gap;
		this.drop = drop;
		
		// initialize the state of the extender
		this.clear();
	}
	/**
	 * extend the two-hit to get the alignment 
	 * to get the alignment result call the .getAlignment() method
	 * @param ThresholdScore
	 * @return true if the attained score is equivalent to 
	 * or greater than the thresholdScore, false else.  
	 * @throws Exception
	 */
	public boolean extend(int thresholdScore) throws Exception{
		//////////////extend////////////////////////////////////////////////
		// @debug
		long timeExtensionStart = System.currentTimeMillis();
		if(DEBUG){
			System.out.println("  >>>X-drop scoring......");
			//TODO:
			System.out.println("   YBaseIndex:"+this.YBaseIndex+"XBaseIndex:"+XBaseIndex);
			System.err.println("   YBaseIndex:"+this.YBaseIndex+"XBaseIndex:"+XBaseIndex);
		}
		// ~debug
		
		//X-drop extend rightwards, gain scoreRight
		//X-drop extend leftwards, gain scoreLeft
		xDropScore();
		if(leftScore+rightScore-match < thresholdScore){
			return false;
		}
		
		//@debug
		long timeExtensionFinished = System.currentTimeMillis();
		if(DEBUG){
			System.out.println("   left score:["+leftScoreXInMatrix+","+leftScoreYInMatrix+"]="+leftScore
					+"\n   right score:["+rightScoreXInMatrix+","+rightScoreYInMatrix+"]:"+rightScore);
			System.out.println("   total score(left+right-match)):"+
					(leftScore+rightScore-match));
			System.out.println("   boundary array size, Upper:"+
					boundaryYUpper.size()+" Lower:"+boundaryYLower.size());
			System.out.println("   [extension time costs: "
					+(timeExtensionFinished-timeExtensionStart)+" milliseconds]");
		}
		if(DEBUG_WRITE_PATH){
			PrintWriter boundaryWriter=null;
			try {
				boundaryWriter = new PrintWriter("pointsBoundary.txt");
				for(int i=leftScoreXInMatrix; i<=rightScoreXInMatrix; i++){
					boundaryWriter.println(i+" "+boundaryYUpper.get(i));
					boundaryWriter.println(i+" "+boundaryYLower.get(i));
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}finally{
				boundaryWriter.close();
			}
		}
		//~debug
		/////////////////////////////////////////////////////////////
		
		////////////////trace back//////////////////////////////////
		long timeTracBackStart = System.currentTimeMillis();
		if(DEBUG){
			System.out.println("  >>>trace back path......");
		}
		
		//Hirsberg's linear space trace back path of right part
		//Hirsberg's linear space trace back path of left part
		traceBackPath();
		
		if(DEBUG){
			long timeTraceBackFinished = System.currentTimeMillis();
			System.out.println("   [trace back time costs: "
				+(timeTraceBackFinished-timeTracBackStart)+" milliseconds]");
			System.out.println("[total time costs:"+
				(timeTraceBackFinished-timeExtensionStart)+"]");
		}
		///////////////////////////////////////////////////////////
		return true;
	}
	
	/**
	 * the extension score results are saved in class fields
	 * 
	 */
	private void xDropScore(){
		
		/***** Zheng Zhang's X-Drop bounded dynamic programming gapped extension ***/
		// "Zheng Zhang, el. A Greedy Algorithm for Aligning DNA Sequences, 2000"
		// to read the codes, read the paper first
		
		// workspace
		// <(i,j), optimal score of (i,j)>
		HashMap<String, Integer> antidiagonalK0 = new HashMap<String, Integer>();
		HashMap<String, Integer> antidiagonalK1 = new HashMap<String, Integer>();
		HashMap<String, Integer> antidiagonalK2 = new HashMap<String, Integer>();
		// the initial point
		boundaryYUpper.put(0, 0);
		boundaryYLower.put(0, 0);
		
		//////////the base-right-hand part//////////
		//in the first quadrant
		//M, N correspond to the M, N in the paper
		//string A a1a2...ai, string B b1b2...bj, i<=M, j<=N
		//i.e. the x, y limit of the matrix
		int M_Right = XEnd - XBaseIndex;
		int N_Right = YEnd - YBaseIndex;
		// best score so far
		int T_RightBestScore = 0;
		int T_RightScoreXCoordinate = 0;
		int T_RightScoreYCoordinate = 0;
		// current Score
		int T_RightCurrentScore = 0;
		// bound for X-coordinate
		int L_Right = 0;
		int U_Right = 0;
		// anti-diagonal k
		int k = 0;
		// the initial anti-diagonal
		antidiagonalK0.put("0,0", 0);
		while (L_Right <= U_Right + 1) {
			k += 1;
			int L_RightNext = Integer.MAX_VALUE;
			int U_RightNext = Integer.MIN_VALUE;
			
			// swap the workspaces
			HashMap<String, Integer> tempK;
			tempK = antidiagonalK0;
			antidiagonalK0 = antidiagonalK2;
			antidiagonalK0.clear();
			antidiagonalK2 = antidiagonalK1;
			antidiagonalK1 = tempK;

			int xCoordinate = 0;
			int yCoordinate = 0;

			for (xCoordinate = L_Right; xCoordinate <= U_Right + 1; xCoordinate++) {
				yCoordinate = k - xCoordinate;
				int score = Integer.MIN_VALUE;
				// three cases
				// case1: come from diagonal
				Integer diagnal = antidiagonalK2
						.get((xCoordinate - 1 + "," + (yCoordinate - 1)));
				if (diagnal != null) {
					int diagnoalScore = diagnal;
					// note: XStr, YStr index from 0, so :XBegin+xCoordinate-1
					int misOrMatchScore = XStr[XBaseIndex + xCoordinate - 1] == YStr[YBaseIndex
							+ yCoordinate - 1] ? match : mismatch;
					score = diagnoalScore + misOrMatchScore > score ? 
							diagnoalScore + misOrMatchScore	: score;
				}
				// case2: come from below
				Integer below = antidiagonalK1
						.get((xCoordinate + "," + (yCoordinate - 1)));
				if (below != null) {
					int yScore = below;
					score = score < yScore + gap ? yScore + gap : score;
				}
				// case3: come from left
				Integer left = antidiagonalK1.get((xCoordinate - 1) + ","
						+ (yCoordinate));
				if (left != null) {
					int xScore = left;
					score = score < xScore + gap ? xScore + gap : score;
				}
				if (score > T_RightCurrentScore) {
					T_RightCurrentScore= score;
					T_RightScoreXCoordinate= xCoordinate;
					T_RightScoreYCoordinate= yCoordinate;
				}
				T_RightCurrentScore = T_RightCurrentScore > score ? T_RightCurrentScore : score;

				if (score >= T_RightBestScore - drop) {
					antidiagonalK0.put(xCoordinate + "," + yCoordinate, score);
					L_RightNext = L_RightNext < xCoordinate ? L_RightNext : xCoordinate;
					U_RightNext = U_RightNext > xCoordinate ? U_RightNext : xCoordinate;
					//record boundary
					if(boundaryYLower.get(xCoordinate)!=null){
						int YLower = boundaryYLower.get(xCoordinate);
						YLower = YLower > yCoordinate? yCoordinate : YLower;
						boundaryYLower.put(xCoordinate, YLower);
					}else{
						boundaryYLower.put(xCoordinate, yCoordinate);
					}
					if(boundaryYUpper.get(xCoordinate)!=null){
						int YUpper= boundaryYUpper.get(xCoordinate);
						YUpper= YUpper < yCoordinate? yCoordinate : YUpper;
						boundaryYUpper.put(xCoordinate, YUpper);
					}else{
						boundaryYUpper.put(xCoordinate, yCoordinate);
					}
				}
			}
			// band closed, extension terminated
			if (antidiagonalK0.size() == 0) {
				break;
			}
			// update L, U
			L_Right = L_RightNext;
			U_Right = U_RightNext;
			// keep to legitimate point
			L_Right = L_Right > k + 1 - N_Right ? L_Right : k + 1 - N_Right;
			U_Right = U_Right < M_Right - 1 ? U_Right : M_Right - 1;
			// update TBestScore
			T_RightBestScore = T_RightCurrentScore;
		}//while (L <= U + 1)
		//save the results
		this.rightScoreXInMatrix = T_RightScoreXCoordinate;
		this.rightScoreYInMatrix = T_RightScoreYCoordinate;
		this.rightScore = T_RightBestScore;
		///////////////////////////////////////////////////////////////////////////
		
		////////////the base-left-hand part////////////////////////
		//in the forth quadrant
		//the left most coordinate (usually be negative)
		int M_Left= XBegin- XBaseIndex-1;
		//the bottom most coordinate (usually be negative)
		int N_Left = YBegin - YBaseIndex-1;
		//best score so far
		int T_LeftBestScore = 0;
		int T_LeftScoreXCoordinate = 0;
		int T_LeftScoreYCoordinate = 0;
		//current Score
		int T_LeftCurrentScore = 0;
		//bound for X-coordinate
		int L_Left= 0;
		int U_Left = 0;
		//anti-diagonal k
		int kLeft = 0;
		antidiagonalK0.clear();
		antidiagonalK1.clear();
		antidiagonalK2.clear();
		
		//the initial anti-diagonal
		antidiagonalK0.put("0,0", 0);
		while (L_Left >= U_Left-1) {
			//next anti-diagonal
			kLeft -= 1;
			int L_LeftNext = Integer.MIN_VALUE;
			int U_LeftNext = Integer.MAX_VALUE;
			
			// swap the workspaces
			HashMap<String, Integer> tempK;
			tempK = antidiagonalK0;
			antidiagonalK0 = antidiagonalK2;
			antidiagonalK0.clear();
			antidiagonalK2 = antidiagonalK1;
			antidiagonalK1 = tempK;

			int xCoordinate = 0;
			int yCoordinate = 0;
			for (xCoordinate = L_Left; xCoordinate >= U_Left - 1; xCoordinate--) {
				yCoordinate = kLeft - xCoordinate;
				int score = Integer.MIN_VALUE;
				// three cases
				// case1: come from diagonal
				Integer diagnal = antidiagonalK2.get(((xCoordinate+1)+","+(yCoordinate+1)));
				if (diagnal != null) {
					int diagnoalScore = diagnal;
					// note: XStr, YStr index from 0,
					// xCoordinate: 0 ,-1, -2, ....
					// so :XBegin+xCoordinate+1
					int misOrMatchScore = XStr[XBaseIndex + xCoordinate+1] == 
						YStr[YBaseIndex	+ yCoordinate+1] ? match : mismatch;
					score = diagnoalScore + misOrMatchScore > score ? 
							diagnoalScore + misOrMatchScore : score;
				}
				// case2: come from below
				Integer below = antidiagonalK1
						.get((xCoordinate + "," + (yCoordinate + 1)));
				if (below != null) {
					int yScore = below;
					score = score < yScore + gap ? yScore + gap : score;
				}
				// case3: come from left
				Integer left = antidiagonalK1.get((xCoordinate + 1) + ","
						+ (yCoordinate));
				if (left != null) {
					int xScore = left;
					score = score < xScore + gap ? xScore + gap : score;
				}
				if (score > T_LeftCurrentScore) {
					T_LeftCurrentScore= score;
					T_LeftScoreXCoordinate= xCoordinate;
					T_LeftScoreYCoordinate= yCoordinate;
				}
				T_LeftCurrentScore = T_LeftCurrentScore > score ? T_LeftCurrentScore : score;
				if (score >= T_LeftBestScore - drop) {
					antidiagonalK0.put(xCoordinate + "," + yCoordinate, score);
					//record boundary
					if(boundaryYLower.get(xCoordinate)!=null){
						int YLower = boundaryYLower.get(xCoordinate);
						YLower = YLower > yCoordinate? yCoordinate : YLower;
						boundaryYLower.put(xCoordinate, YLower);
					}else{
						boundaryYLower.put(xCoordinate, yCoordinate);
					}
					if(boundaryYUpper.get(xCoordinate)!=null){
						int YUpper= boundaryYUpper.get(xCoordinate);
						YUpper= YUpper < yCoordinate? yCoordinate : YUpper;
						boundaryYUpper.put(xCoordinate, YUpper);
					}else{
						boundaryYUpper.put(xCoordinate, yCoordinate);
					}
					L_LeftNext = L_LeftNext > xCoordinate ? L_LeftNext : xCoordinate;
					U_LeftNext = U_LeftNext < xCoordinate ? U_LeftNext : xCoordinate;
				}
			}
			// band closed, extension terminated
			if (antidiagonalK0.size() == 0) {
				break;
			}
			// update L, U
			L_Left = L_LeftNext;
			U_Left= U_LeftNext;
			// keep to legitimate point
			L_Left = L_Left < kLeft-1-N_Left ? L_Left: kLeft-1-N_Left;
			U_Left = U_Left > M_Left+1 ? U_Left: M_Left+1;
			// update TBestScore
			T_LeftBestScore = T_LeftCurrentScore;
		}//while (L <= U + 1)
		// save the results
		this.leftScoreXInMatrix = T_LeftScoreXCoordinate;
		this.leftScoreYInMatrix = T_LeftScoreYCoordinate;
		this.leftScore = T_LeftBestScore;
		////////////////////////////////////////////////////////////////
	}
	
	/**
	 * the path points are saved in class field
	 * @throws Exception 
	 */
	private void traceBackPath() throws Exception{
		
		/******* Hirsberg's divide and conquer linear space sequence alignment ****/
		// trace back the alignment path
		// in a non-iterated method with a queue of area list
		// workspace
		int WORK_SPACE_SIZE=rightScoreYInMatrix > -leftScoreYInMatrix ? 
				rightScoreYInMatrix+1 : (-leftScoreYInMatrix)+1;
		int[] Acolumn = new int[WORK_SPACE_SIZE];
		int[] Bcolumn = new int[WORK_SPACE_SIZE];
		int[] Ccolumn = new int[WORK_SPACE_SIZE];
		
		// a queue for areas to be dealt with
		LinkedList<Area> queueOfArea = new LinkedList<Area>();
		
		//put the initial area in the queue 
		queueOfArea.addLast(new Area(0, rightScoreXInMatrix, 0, rightScoreYInMatrix, rightScore));
		queueOfArea.addLast(new Area(leftScoreXInMatrix, 0, leftScoreYInMatrix, 0, leftScore));
		
		//pop the queue, process the areas
		while(queueOfArea.size()>0){
			Area area = (Area)(queueOfArea.removeFirst());
			if(DEBUG_VERBOSE){
				System.out.println("   while: get area:"+area.toSting());
			}
			if(area.Y2==area.Y1){
				for(int index1=area.X1; index1<=area.X2; index1++){
					criticalPathRelative.add(index1+","+area.Y1);
					//@debug
					if(DEBUG_WRITE_PATH){
						tracebackWriter.println(index1+" "+area.Y1);
					}
				}
			}else if(area.X2==area.X1){
				for(int index2=area.Y1; index2 <= area.X2; index2++){
					criticalPathRelative.add(area.X1+","+index2);
					//@debug
					if(DEBUG_WRITE_PATH){
						tracebackWriter.println(area.X1+" "+index2);
					}
				}
				// when area is small enough
			}else if((area.X2-area.X1+1)<=300 ||(area.Y2-area.Y1+1)<=300){
				int[][] matrix = new int[area.X2-area.X1+1][area.Y2-area.Y1+1];
				for(int indexY=0; indexY<matrix[0].length; indexY++){
					matrix[0][indexY] = indexY*gap;
				}
				// fill the matrix
				// in the first quadrant
				for(int xIndex=1; xIndex<matrix.length; xIndex++){
					matrix[xIndex][0]=xIndex*gap;
					for(int yIndex=1;yIndex<matrix[0].length; yIndex++){
						//three cases:
						//from diagonal
						int alphaXY ;
						if(area.X2>0){
							alphaXY = XStr[XBaseIndex+area.X1+xIndex-1]==
								YStr[YBaseIndex+area.Y1+yIndex-1]?match:mismatch;
						}else{
							alphaXY = XStr[XBaseIndex+area.X2-xIndex+1]==
								YStr[YBaseIndex+area.Y2-yIndex+1]?match:mismatch;
						}
						int optimalScore = matrix[xIndex-1][yIndex-1]+alphaXY;
						//from left
						optimalScore = optimalScore < matrix[xIndex-1][yIndex]+gap?
								matrix[xIndex-1][yIndex]+gap : optimalScore;
						//from below
						optimalScore = optimalScore < matrix[xIndex][yIndex-1]+gap?
								matrix[xIndex][yIndex-1]+gap : optimalScore;
						matrix[xIndex][yIndex] = optimalScore;
					}
				}
				//trace back the path
				//last point<x,y>
				int x=matrix.length-1;
				int y=matrix[0].length-1;
				criticalPathRelative.add((x+area.X1)+","+(y+area.Y1));
				while(x!=0 || y!=0){
					//three cases:
					//from diagonal
					if(x-1>=0 && y-1>=0){
						int alphaXY2;
						if(area.X2>0){
							alphaXY2 = XStr[XBaseIndex+area.X1+x-1]==
								YStr[YBaseIndex+area.Y1+y-1]?match:mismatch;
						}else{
							alphaXY2 = XStr[XBaseIndex+area.X2-x+1]==
								YStr[YBaseIndex+area.Y2-y+1]?match:mismatch;
						}
						if(matrix[x][y]==alphaXY2+matrix[x-1][y-1]){
							x=x-1;
							y=y-1;
							criticalPathRelative.add((x+area.X1)+","+(y+area.Y1));
							continue;
						}
					}
					//from left
					if(x-1 >= 0){
						if(matrix[x][y]==gap+matrix[x-1][y]){
							x=x-1;
							criticalPathRelative.add((x+area.X1)+","+(y+area.Y1));
							continue;
						}
					}
					//from below
					if(y-1 >= 0){
						if(matrix[x][y]==gap+matrix[x][y-1]){
							y=y-1;
							criticalPathRelative.add((x+area.X1)+","+(y+area.Y1));
							continue;
						}
					}
				}//while(x!=0 || y!=0)
				
				// divide the area into two sub-areas
			}else{
				Area[] twoAreas = divideIntoTwoAreas(area, Acolumn, Bcolumn, Ccolumn);
				if(DEBUG_VERBOSE){
					System.out.println("divide into:\n"+twoAreas[0].toSting()+"\n"+
							twoAreas[1].toSting());
				}
				queueOfArea.addLast(twoAreas[0]);
				queueOfArea.addLast(twoAreas[1]);
				//criticalPathRelative.add(twoAreas[0].X2+" "+twoAreas[0].Y2);
				
			}//while(queueOfArea.size()>0)
		}
		//@debug
		if(DEBUG_WRITE_PATH){
			for(String point: criticalPathRelative){
				criticalPathWriter.println(point);
			}
			criticalPathWriter.close();
			tracebackWriter.close();
		}
		//~debug

	}

		/**
		 * @param area must be at least of three columns
		 * @param T_RightScoreYCoordinate
		 * @param T_LeftScoreYCoordinate
		 * @return two sub-areas 
		 * notes: the left sub-area must be on the [0]
		 * @throws Exception 
		 * 
		 */
	private Area[] divideIntoTwoAreas(Area area, int[] Acolumn, int[] Bcolumn, int[] Ccolumn ) throws Exception{
		Area[] twoAreas = new Area[2];
		// XStr middle index
		int Xmid = (area.X1+area.X2)/2;
		if(DEBUG_VERBOSE){
			System.out.println("   X-mid:"+Xmid);
		}
		/////////////////////fill the left half matrix///////////////////////
		int YLowerBound = area.Y1 > boundaryYLower.get(area.X1) ? area.Y1 : boundaryYLower.get(area.X1);
		int YUpperBound = area.Y2 < boundaryYUpper.get(area.X1)	? area.Y2 : boundaryYUpper.get(area.X1);
		int preYLowerBound;
		int preYUpperBound;
		if(DEBUG_VERBOSE){
			System.out.println("   linear trace back:");
		}
		// fill the initial column: area.X1
		// in the first quadrant
		if(area.X2>0){
			for(int indexB=YLowerBound; indexB<=YUpperBound; indexB++ ){
				Bcolumn[indexB] = (indexB-area.Y1)*gap;
				//@debug
				if(DEBUG_VERBOSE){
					System.out.println("["+area.X1+","+indexB+"]:"+Bcolumn[indexB]);
				}
				if(DEBUG_WRITE_PATH){
					tracebackWriter.println(area.X1+" "+indexB);
				}
				//~debug
			}
			// in the third quadrant
		}else{
			for(int indexB=YLowerBound; indexB <=YUpperBound; indexB++ ){
				Bcolumn[-indexB]= (indexB-area.Y1)*gap;
				//@debug
				if(DEBUG_VERBOSE){
					System.out.println("["+area.X1+","+indexB+"]:"+Bcolumn[-indexB]);
				}
				if(DEBUG_WRITE_PATH){
					tracebackWriter.println(area.X1+" "+indexB);
				}
				//~debug
			}
		}
		
		preYLowerBound=YLowerBound;
		preYUpperBound=YUpperBound;
		// fill the other columns: <X1+1 -- Xmid>
		for(int indexX = area.X1+1; indexX<=Xmid; indexX++){
			//swap A column and B column
			int[] tempCol = Bcolumn; Bcolumn=Acolumn; Acolumn=tempCol;
			//update boundary
			YLowerBound = area.Y1 > boundaryYLower.get(indexX) ? area.Y1 : boundaryYLower.get(indexX);
			YUpperBound = area.Y2 < boundaryYUpper.get(indexX) ? area.Y2 : boundaryYUpper.get(indexX);

			// fill column: indexX
			for(int indexB2=YLowerBound; indexB2<=YUpperBound; indexB2++){
				int alphaXY;
				//////////// in the first quadrant
				if(area.X2>0){
//					System.out.println("<"+indexX+","+indexB2+">: Upper:"+YUpperBound
//							+" Lower:"+YLowerBound
//							+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
//							+" indexB2-1:"+(indexB2-1)
//							+" indexB2"+indexB2);
					boolean taken = false;
					Bcolumn[indexB2]=Integer.MIN_VALUE;
					//three cases: value from left, below, diagonal
//					System.out.print("case 1:");
					//case 1: from diagonal
					if(indexB2-1>=preYLowerBound && indexB2-1 <=preYUpperBound){
//						System.out.println("["+indexX+","+indexB2+"]=>"+
//						"["+(XBaseIndex+indexX-1)+","+(YBaseIndex+indexB2-1)+"]"
//						+XStr[XBaseIndex+indexX-1]+"=="+YStr[YBaseIndex+indexB2-1]);

						alphaXY=XStr[XBaseIndex+indexX-1]==
							YStr[YBaseIndex+indexB2-1]? match : mismatch;
						Bcolumn[indexB2] = alphaXY+Acolumn[indexB2-1];
//						System.out.print(alphaXY+Acolumn[indexB2-1]);
						taken = true;
					}
//					System.out.println();
					//case 2: from left
//					System.out.print("case 2: ");
					if(indexB2<=preYUpperBound && indexB2>=preYLowerBound){
						Bcolumn[indexB2] = Bcolumn[indexB2] > gap+Acolumn[indexB2]? 
									Bcolumn[indexB2] : gap+Acolumn[indexB2];
//						System.out.print(gap+Acolumn[indexB2]);
						taken = true;
					}
//					System.out.println();
					//case3: from below
//					System.out.print("case 3: ");
					if(indexB2-1<=YUpperBound && indexB2-1 >= YLowerBound){
						Bcolumn[indexB2]=Bcolumn[indexB2] > gap+Bcolumn[indexB2-1]?
									Bcolumn[indexB2] : gap+Bcolumn[indexB2-1];
//						System.out.print(gap+Bcolumn[indexB2-1]);
						taken = true;
					}
					if(DEBUG_VERBOSE){
						System.out.println("["+indexX+","+indexB2+"]:"+Bcolumn[indexB2]);
					}
					if(!taken){
//						YLowerBound++;
						YLowerBound=indexB2+1;
						
						//TODO:
						throw new Exception("non of three conditions is taken!\n"+
								"<"+indexX+","+indexB2+">: Upper:"+YUpperBound
								+" Lower:"+YLowerBound
								+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
								+" indexB2:"+indexB2);
						
					}
					if(DEBUG_WRITE_PATH&&taken){
						tracebackWriter.println(indexX+" "+indexB2);
					}
//					System.out.println("=======================");
					
//					System.out.println("["+indexX+","+indexB2+"]"+"=>["
//							+(XBaseIndex+indexX+1)+","+(YBaseIndex+indexB2+1)+
//							"]: "+XStr[XBaseIndex+indexX]+"=="+YStr[YBaseIndex+indexB2]);
					////////////in the third quadrant////////////
				}else{
					boolean taken = false;
					Bcolumn[-indexB2]=Integer.MIN_VALUE;
					//three cases: value from left, below, diagonal
//					System.out.print("case 1:");
					//case 1: from diagonal
					if(indexB2-1>=preYLowerBound && indexB2-1 <=preYUpperBound){
//					System.out.println("["+indexX+","+indexB2+"]"+"=>["
//							+(XBaseIndex+indexX)+","+(YBaseIndex+indexB2)+
//							"]: "+XStr[XBaseIndex+indexX]+"=="+YStr[YBaseIndex+indexB2]);
						//be caution here,"XBaseIndex+indexX" only, "+1" is needn't 
						alphaXY=XStr[XBaseIndex+indexX]==
							YStr[YBaseIndex+indexB2]? match : mismatch;
						Bcolumn[-indexB2] = alphaXY+Acolumn[-(indexB2-1)];
						taken = true;
//						System.out.println("alpha:"+alphaXY);
//						System.out.print(alphaXY+Acolumn[-indexB2-1]);
					}
					//case 2: from left
//					System.out.print("case 2: ");
					if(indexB2<=preYUpperBound && indexB2>=preYLowerBound){
						Bcolumn[-indexB2] = Bcolumn[-indexB2] > gap+Acolumn[-indexB2]? 
							Bcolumn[-indexB2] : gap+Acolumn[-indexB2];
						taken = true;
//						System.out.print(gap+Acolumn[-indexB2]);
					}
//					System.out.println();
					//case3: from below
//					System.out.print("case 3: ");
					if(indexB2-1<=YUpperBound && indexB2-1 >= YLowerBound){
						Bcolumn[-indexB2]=Bcolumn[-indexB2] > gap+Bcolumn[-(indexB2-1)]?
							Bcolumn[-indexB2] : gap+Bcolumn[-(indexB2-1)];
						taken = true;
//						System.out.print(gap+Bcolumn[indexB2-1]);
					}
					if(DEBUG_VERBOSE){
						System.out.println("["+indexX+","+indexB2+"]:"+Bcolumn[-indexB2]);
					}
					if(!taken){
						//TODO:
//						YLowerBound++;
//						System.err.println("non of three conditions is taken!"+
//								"<"+indexX+","+indexB2+">: Upper:"+YUpperBound
//								+" Lower:"+YLowerBound
//								+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
//								+" indexB2:"+indexB2);
						YLowerBound=indexB2+1;
						continue;
						//re-limit the upper bound
//						throw new Exception("non of three conditions is taken!"+
//								"<"+indexX+","+indexB2+">: Upper:"+YUpperBound
//								+" Lower:"+YLowerBound
//								+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
//								+" indexB2:"+indexB2);
					}
					if(DEBUG_WRITE_PATH&&taken){
						tracebackWriter.println(indexX+" "+indexB2);
					}

				}
			}
			preYLowerBound = YLowerBound;
			preYUpperBound = YUpperBound;
		}//for(int indexX = area.X1+1; indexX<=Xmid; indexX++)
		
		
		//////////////////////// fill the right half matrix ////////////////////////////
		
//		System.out.println("right half");
		// update boundaries
		YLowerBound = area.Y1 > boundaryYLower.get(area.X2) ? area.Y1 : boundaryYLower.get(area.X2);
		YUpperBound = area.Y2 < boundaryYUpper.get(area.X2) ? area.Y2 : boundaryYUpper.get(area.X2); 
		// the initial column area.X2
		// in the first quadrant
		if(area.X2>0){
			for(int indexC=YUpperBound; indexC>=YLowerBound; indexC--){
				Ccolumn[indexC] = (-(indexC-area.Y2))*gap;
				if(DEBUG_VERBOSE){
					System.out.println("["+area.X2+","+indexC+"]:"+Ccolumn[indexC]);
				}
				if(DEBUG_WRITE_PATH){
					tracebackWriter.println(area.X2+" "+indexC);
				}
			}
			// in the third quadrant
		}else{
			for(int indexC=YUpperBound; indexC>=YLowerBound; indexC--){
				Ccolumn[-indexC] = (-(indexC-area.Y2))*gap;
				if(DEBUG_VERBOSE){
					System.out.println("["+area.X2+","+indexC+"]:"+Ccolumn[-indexC]);
				}
				if(DEBUG_WRITE_PATH){
					tracebackWriter.println(area.X2+" "+indexC);
				}
			}
		}
		preYLowerBound = YLowerBound;
		preYUpperBound = YUpperBound;
		// fill columns: <Xmid -- X2-1>
		for(int indexX2=area.X2-1; indexX2>=Xmid; indexX2--){
			// swap A and C
			int[] tempCol2=Ccolumn; Ccolumn = Acolumn; Acolumn=tempCol2;
			// update boundaries
			YLowerBound = area.Y1 > boundaryYLower.get(indexX2) ? area.Y1 : boundaryYLower.get(indexX2);
			YUpperBound = area.Y2 < boundaryYUpper.get(indexX2) ? area.Y2 : boundaryYUpper.get(indexX2);
			for(int indexC2=YUpperBound; indexC2>=YLowerBound; indexC2--){
				int alphaXY;
				//three cases: value from right, upper, diagonal
				//////////// in the first quadrant
				if(area.X2>0){
//					System.out.println("<"+indexX2+","+indexC2+">: Upper:"+YUpperBound
//							+" Lower:"+YLowerBound
//							+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
//							+" indexC2-1:"+(indexC2-1)
//							+" indexC2:"+indexC2);
					boolean taken = false;
					Ccolumn[indexC2]=Integer.MIN_VALUE;
//					System.out.print("case 1: ");
					//case 1: from diagonal
					if(indexC2+1<=preYUpperBound && indexC2+1 >= preYLowerBound){
//						System.out.println(indexX2+","+indexC2+"=>"+
//							"["+(XBaseIndex+indexX2)+","+(YBaseIndex+indexC2)+"]:"
//							+XStr[XBaseIndex+indexX2]+"=="+YStr[YBaseIndex+indexC2]);
						// be caution here, "[XBaseIndex+indexX2" only, "-1" is not needed 
						alphaXY=XStr[XBaseIndex+indexX2]==
							YStr[YBaseIndex+indexC2]? match : mismatch;
						Ccolumn[indexC2] = alphaXY+Acolumn[indexC2+1];
						taken = true;
//						System.out.print(alphaXY+Acolumn[indexC2+1]);
//						System.out.println("form diagonal:"+Ccolumn[indexC2]);
					}
//					System.out.println();
//					System.out.print("case 2:");
					//case 2: from right
					if(indexC2>=preYLowerBound && indexC2 <= preYUpperBound){
//						System.out.println("pre L~U:"+preYLowerBound+"~"+preYUpperBound);
//						System.out.print(gap+Acolumn[indexC2]);
						Ccolumn[indexC2] = Ccolumn[indexC2] > gap+Acolumn[indexC2] ?
								Ccolumn[indexC2] : gap+Acolumn[indexC2];
								taken = true;
//						System.out.println("form right:"+Ccolumn[indexC2]);
					}
//					System.out.println();
					//case3: from upper
//					System.out.print("case 3:");
					if(indexC2+1 <=YUpperBound && indexC2+1 >= YLowerBound){
//						System.out.print(gap+Ccolumn[indexC2+1]);
						
						Ccolumn[indexC2] = Ccolumn[indexC2] > gap+Ccolumn[indexC2+1] ?
								Ccolumn[indexC2] : gap+Ccolumn[indexC2+1];
						taken = true;
//						System.out.println("form upper:"+Ccolumn[indexC2]
//						    +"= gap:"+gap+"+upper["+(indexC2+1)");
					}
//					System.out.println();
					if(DEBUG_VERBOSE){
						System.out.println("["+indexX2+","+indexC2+"]:"+Ccolumn[indexC2]);
					}
					if(!taken){
						//re-limit the upper bound
//						YUpperBound--;
						//TODO:
						YUpperBound=indexC2-1;
//						System.out.println("<"+indexX2+","+indexC2+">: Upper:"+YUpperBound
//								+" Lower:"+YLowerBound
//								+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
//								+" indexC2-1:"+(indexC2-1)
//								+" indexC2:"+indexC2);
						if(DEBUG_VERBOSE){
							System.out.println("new UpperBound:"+YUpperBound);
						}
						continue;
					}
					if(DEBUG_WRITE_PATH&&taken){
						tracebackWriter.println(indexX2+" "+indexC2);
					}
					////////////in the third quadrant////////////
				}else{
//					System.out.println("<"+indexX2+","+indexC2+">: Upper:"+YUpperBound
//					+" Lower:"+YLowerBound
//					+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
//					+" indexB2-1:"+(indexC2-1)
//					+" indexB2"+indexC2);
					boolean taken = false;
					Ccolumn[-indexC2]=Integer.MIN_VALUE;
//					System.out.print("case 1: ");
					//case 1: from diagonal
					if(indexC2+1<=preYUpperBound && indexC2+1 >= preYLowerBound){
//						System.out.println(indexX2+","+indexC2+"=>"+
//							"["+(XBaseIndex+indexX2+1)+","+(YBaseIndex+indexC2+1)+"]");
						alphaXY=XStr[XBaseIndex+indexX2+1]==
							YStr[YBaseIndex+indexC2+1]? match : mismatch;
						Ccolumn[-indexC2] = alphaXY+Acolumn[-(indexC2+1)];
						taken = true;
//						System.out.print(alphaXY+Acolumn[-indexC2+1]);
					}
//					System.out.println();
//					System.out.print("case 2:");
					//case 2: from right
					if(indexC2>=preYLowerBound && indexC2 <= preYUpperBound){
//						System.out.print(gap+Acolumn[indexC2]);
						Ccolumn[-indexC2] = Ccolumn[-indexC2] > gap+Acolumn[-indexC2] ?
								Ccolumn[-indexC2] : gap+Acolumn[-indexC2];
						taken = true;
					}
//					System.out.println();
					//case3: from upper
//					System.out.print("case 3:");
					if(indexC2+1 <=YUpperBound && indexC2+1 >= YLowerBound){
//						System.out.print(gap+Ccolumn[indexC2+1]);
						Ccolumn[-indexC2] = Ccolumn[-indexC2] > gap+Ccolumn[-(indexC2+1)] ?
								Ccolumn[-indexC2] : gap+Ccolumn[-(indexC2+1)];
						taken = true;
					}
					if(DEBUG_VERBOSE){
						System.out.println("["+indexX2+","+indexC2+"]:"+Ccolumn[-indexC2]);
					}
					if(!taken){
						//re-limit the upper bound
//						YUpperBound--;
						//TODO:
						YUpperBound=indexC2-1;
						if(DEBUG_VERBOSE){
							System.out.println("<"+indexX2+","+indexC2+">: Upper:"+YUpperBound
									+" Lower:"+YLowerBound
									+" PreUpper,PreLower:"+preYUpperBound+","+preYLowerBound
									+" indexC2-1:"+(indexC2-1)
									+" indexC2:"+indexC2);
							System.out.println("new UpperBound:"+YUpperBound);
						}
					}
					if(DEBUG_WRITE_PATH&&taken){
						tracebackWriter.println(indexX2+" "+indexC2);
					}

				}
			}//for(int indexC2=YUpperBound; indexC2>=YLowerBound; indexC2--)
			preYLowerBound = YLowerBound;
			preYUpperBound = YUpperBound;
		}//for(int indexX2=area.X2-1; indexX2>=Xmid; indexX2--)
		///////////////////////////////////////////////////////
		
		///////// find the critical point//////////////////////
		if(area.X2>0){
			boolean foundCriticalPoint = false;
			for(int criticalIndex = YLowerBound; criticalIndex<=YUpperBound; criticalIndex++){
				if(DEBUG_VERBOSE){
					System.out.println("score:"+area.score+" left["+criticalIndex+"]:"+Bcolumn[criticalIndex]
					    +" + right["+criticalIndex+"]:"+Ccolumn[criticalIndex]+"="
					    +(Bcolumn[criticalIndex]+Ccolumn[criticalIndex]));
				}
				if(Bcolumn[criticalIndex]
				           +Ccolumn[criticalIndex]>=area.score){
					if(DEBUG_VERBOSE){
						System.out.println("split point:<"+Xmid+","+criticalIndex+">");
					}
					twoAreas[0]=new Area(area.X1,Xmid,area.Y1,criticalIndex,
							Bcolumn[criticalIndex]);
					twoAreas[1]=new Area(Xmid,area.X2,criticalIndex,area.Y2,
							Ccolumn[criticalIndex]);
					foundCriticalPoint=true;
					break;
				}
			}
			
			if(!foundCriticalPoint){
				throw new Exception("critical point not found");
			}
			
		}else{
			boolean foundCriticalPoint = false;
			for(int criticalIndex = YLowerBound; criticalIndex<=YUpperBound; criticalIndex++){
				if(DEBUG_VERBOSE){
					System.out.println("left["+criticalIndex+"]:"+Bcolumn[-criticalIndex]
					       +" + right["+criticalIndex+"]:"+Ccolumn[-criticalIndex]+"="
					       +(Bcolumn[-criticalIndex]+Ccolumn[-criticalIndex]));
				}
				if(Bcolumn[-criticalIndex]
				           +Ccolumn[-criticalIndex]>=area.score){
					if(DEBUG_VERBOSE){
						System.out.println("split point:<"+Xmid+","+criticalIndex+">");
					}
					twoAreas[0]=new Area(area.X1,Xmid,area.Y1,criticalIndex,
							Bcolumn[-criticalIndex]);
					twoAreas[1]=new Area(Xmid,area.X2,criticalIndex,area.Y2,
							Ccolumn[-criticalIndex]);
					foundCriticalPoint=true;
					break;
				}
			}
			if(!foundCriticalPoint){
				throw new Exception("critical point not found"+
						area.toSting());
			}
		}
		if(DEBUG_WRITE_PATH){
			tracebackWriter.flush();
		}
		return twoAreas;
	}
		
	/**
	 * this method must be called after dpExtend() is called
	 * @return
	 * @throws Exception 
	 */
	public Alignment getAlignment() throws Exception{
		Alignment result=null;
		if(DEBUG){
			System.out.println("path set size:"+criticalPathRelative.size());
		}
		if(criticalPathRelative.size()==0){
			return null;
		}
		byte[] XAlignResult = new byte[this.criticalPathRelative.size()-1];
		byte[] YAlignResult = new byte[this.criticalPathRelative.size()-1];
		int Xindex=0;
		int Yindex=0;
		int x=leftScoreXInMatrix;
		int y=leftScoreYInMatrix;
		//<xLeft,yleft> to <0,0>
		while(x<0 || y<0){
			///match or mismatch
			if(this.criticalPathRelative.contains((x+1)+","+(y+1))){
				XAlignResult[Xindex++] = BlastUtils.decode_bp(this.XStr[XBaseIndex+x+1],false);
				YAlignResult[Yindex++] = BlastUtils.decode_bp(this.YStr[YBaseIndex+y+1],false);
				x++;
				y++;
				// y gap
			}else if(criticalPathRelative.contains((x+1)+","+(y))){
				XAlignResult[Xindex++] = BlastUtils.decode_bp(this.XStr[XBaseIndex+x+1], false);
				YAlignResult[Yindex++] = (byte)('-');
				x++;
				// x gap
			}else if(criticalPathRelative.contains((x)+","+(y+1))){
				XAlignResult[Xindex++] = (byte)('-');
				YAlignResult[Yindex++] = BlastUtils.decode_bp(this.YStr[YBaseIndex+y+1], false);
				y++;
			}else{
				throw new Exception("critical path incorrect!");
			}
		}
		//because the base index overlapped
		if(Xindex>0){
			Xindex--;
			Yindex--;
		}
		// <0,0> to <xRight, yRight>
		x=0; y=0;
		while(x<rightScoreXInMatrix || y<rightScoreYInMatrix){
			// match or mismatch
			if(this.criticalPathRelative.contains((x+1)+","+(y+1))){
				x++;
				y++;
				XAlignResult[Xindex++] = BlastUtils.decode_bp(this.XStr[XBaseIndex+x-1],false);
				YAlignResult[Yindex++] = BlastUtils.decode_bp(this.YStr[YBaseIndex+y-1],false);
				// y gap
			}else if(criticalPathRelative.contains((x+1)+","+(y))){
				x++;
				XAlignResult[Xindex++] = BlastUtils.decode_bp(this.XStr[XBaseIndex+x-1], false);
				YAlignResult[Yindex++] = (byte)('-');
				// x gap
			}else if(criticalPathRelative.contains((x)+","+(y+1))){
				y++;
				XAlignResult[Xindex++] = (byte)('-');
				YAlignResult[Yindex++] = BlastUtils.decode_bp(this.YStr[YBaseIndex+y-1], false);
			}else{
				throw new Exception("critical path incorrect!");
			}
		}
		result = new Alignment(XAlignResult, YAlignResult,leftScore+rightScore-match,
				leftScoreXInMatrix+XBaseIndex+1,rightScoreXInMatrix+XBaseIndex-1,
				leftScoreYInMatrix+YBaseIndex+1,rightScoreYInMatrix+YBaseIndex-1);
		if(DEBUG){
			System.out.println("alingment size:"+(Xindex+1));
		}
		return result;
	}
	
	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
//		System.out.println("Integer.MIN_MAX:"+Integer.MIN_VALUE);
		
		String dbFile =  "querydata/95kb.txt";
		String queryFile = "querydata/95kb.txt";
		int XBaseIndex = 1000;
		int YBaseIndex = 1000;
		
//		String queryFile = "DB.txt";
//		String dbFile =  "query.txt";
//		int XBaseIndex = 79;//80;//800-1;//957;
//		int YBaseIndex = 79;//80;//800-1;//957;
		String query = BlastUtils.readQueryFromLocalFile(queryFile);
		String dbstring = BlastUtils.readQueryFromLocalFile(dbFile);
		
//		String query =    "caGTCGCTCTGCGGAGAGGCTGGCAGATTGAGCCCTGGGAGGTTCTCTCCAGCACTAGCA";
//		String dbstring = "GGTAGAGCCTGGGTGTTCCCTGGCCTGCTAGACTCTCACCAGCACTTGGCCGGTGCTGGG";
//		int XBaseIndex = 2;
//		int YBaseIndex = 1;
//		
		byte[] queryseq = BlastUtils.DNA_String2bytes(query);
		byte[] dbseq = BlastUtils.DNA_String2bytes(dbstring);
		System.out.println("queryseq.length:"+queryseq.length);
		System.out.println("base bp:["+XBaseIndex+"]:"+query.charAt(XBaseIndex)+
				"== ["+YBaseIndex+"]:"+dbstring.charAt(YBaseIndex));
		DPExtender extender = new DPExtender();
		extender.setParameters(queryseq, 0, queryseq.length,
				dbseq, 0, dbseq.length, 
				XBaseIndex, YBaseIndex,
				50, 2,-3, -5);
				//Integer.MAX_VALUE, 2,-2, -2);
//	
//		byte[] X = new byte[]{0x01,0x02,0x03,0x04, 0x01, 0x03,0x02,0x02,0x03};
//		byte[] Y = new byte[]{0x01,0x02,0x03,0x03, 0x01, 0x03,0x02,0x02,0x03};
//			int XBaseIndex = 0;//80;//800-1;//957;
//			int YBaseIndex = 0;//80;//800-1;//957;
//			DPExtender extender = new DPExtender(X, 0, X.length,
//					Y, 0, Y.length, 
//					XBaseIndex, YBaseIndex,
//					100, 2,-3, -5);
		
		System.out.println("extending...");
		extender.extend(100);
		
		System.out.println("getAlignment...");
		Alignment result = extender.getAlignment();
		
		/*if(result!=null){
			System.out.println("the result:");
			for(byte bp: result.getQuerySeq()){
				System.out.print((char)bp);
			}
			System.out.println();
			for(byte bp: result.getSubjectSeq()){
				System.out.print((char)bp);
			}
			System.out.println();
			
			System.out.println("write alignment to file...");
			FileOutputStream out = new FileOutputStream("alignment.txt");
			DataOutput dataOut = new DataOutputStream(out);
			result.write(dataOut);
			out.close();
			
			System.out.println("read alignment from file...");
			Alignment alignment = new Alignment();
			FileInputStream in = new FileInputStream("alignment.txt");
			DataInput dataInput = new DataInputStream(in); 
			alignment.readFields(dataInput);
			in.close();
			
			System.out.println("the result:");
			for(byte bp: alignment.getQuerySeq()){
				System.out.print((char)bp);
			}
			System.out.println();
			for(byte bp: alignment.getSubjectSeq()){
				System.out.print((char)bp);
			}
		}
		*/
	}

}
