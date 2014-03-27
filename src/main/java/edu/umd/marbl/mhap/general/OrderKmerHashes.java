/* 
 * MHAP package
 * 
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2014 by Konstantin Berlin and Sergey Koren
 * University Of Maryland
 * 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */
package edu.umd.marbl.mhap.general;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

import edu.umd.marbl.mhap.utils.FastAlignRuntimeException;
import edu.umd.marbl.mhap.utils.Pair;
import edu.umd.marbl.mhap.utils.Utils;

public class OrderKmerHashes 
{
	private static final class SortableIntPair implements Comparable<SortableIntPair>, Serializable
	{
		public final int x;
		public final int y;
		/**
		 * 
		 */
		private static final long serialVersionUID = 2525278831423582446L;

		public SortableIntPair(int x, int y)
		{
			this.x = x;
			this.y = y;
		}

		@Override
		public int compareTo(SortableIntPair p)
		{
			return Integer.compare(this.x, p.x);

			// int result = Integer.compare(this.x, p.x);

			// if (result!=0)
			// return result;

			// return Integer.compare(this.y, p.y);
		}
	}
	
	private final int[][][] orderedHashes;

	public final static int MAX_ARRAY_SIZE = 1000;
	
	private final static int[][][] allocateMemory(int size)
	{
		int subDivisions = size/MAX_ARRAY_SIZE+1;
		if (size%MAX_ARRAY_SIZE==0)
			subDivisions--;		

		//allocate the memory
		int[][][] completeHash = new int[subDivisions][][];
		
		int division = 0;
		for (int iter = 0; iter < size; iter+=MAX_ARRAY_SIZE)
		{
			completeHash[division] = new int[Math.min(MAX_ARRAY_SIZE, size-iter)][2];
			
			division++;
		}
		
		return completeHash;
	}

	public final static OrderKmerHashes fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			//dos.writeInt(size());
			int hashLength = input.readInt();			
			
			int[][][] orderedHashes = allocateMemory(hashLength);

			int division = 0;
			int i = 0;
			for (int iter=0; iter<hashLength; iter++)
			{
				if (orderedHashes[division].length<=i)
				{
					division++;
					i = 0;
				}

				//dos.writeInt(this.completeHash[iter][iter2]);
				orderedHashes[division][i][0] = input.readInt();
				orderedHashes[division][i][1] = input.readInt();
				
				i++;
			}
			
			return new OrderKmerHashes(orderedHashes);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}

	private OrderKmerHashes(int[][][] orderedHashes)
	{
		this.orderedHashes = orderedHashes;
	}
	
	public OrderKmerHashes(Sequence seq, int kmerSize)
	{
		this.orderedHashes = getFullHashes(seq, kmerSize);
	}

	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream(size()*2);
    DataOutputStream dos = new DataOutputStream(bos);
    
    try
		{
			dos.writeInt(size());
	    for (int iter=0; iter<this.orderedHashes.length; iter++)
	    	for (int iter2=0; iter2<this.orderedHashes[iter].length; iter2++)
	    	{
	    		dos.writeInt(this.orderedHashes[iter][iter2][0]);
	    		dos.writeInt(this.orderedHashes[iter][iter2][1]);
	    	}

			
	    dos.flush();
	    return bos.toByteArray();
		}
    catch (IOException e)
    {
    	throw new FastAlignRuntimeException("Unexpected IO error.");
    }
	}
	
	public int getHash(int index)
	{
		int i1 = index/this.orderedHashes[0].length;
		int i2 = index%this.orderedHashes[0].length;
		
		return this.orderedHashes[i1][i2][0];
	}
	
	private int[][][] storeAsArray(SortableIntPair[] completeHashAsPair)
	{
		//allocate the memory
		int[][][] completeHash = allocateMemory(completeHashAsPair.length);

		int division = 0;
		int i = 0;
		for (int iter = 0; iter < completeHashAsPair.length; iter++)
		{
			if (completeHash[division].length<=i)
			{
				division++;
				i = 0;
			}
			
			completeHash[division][i][0] = completeHashAsPair[iter].x;
			completeHash[division][i][1] = completeHashAsPair[iter].y;
			
			i++;
		}

		return completeHash;
	}
	
	private int[][][] getFullHashes(Sequence seq, int subKmerSize)
	{
		// compute just direct hash of sequence
		int[][] hashes = Utils.computeKmerHashesInt(seq, subKmerSize, 1);

		SortableIntPair[] completeHashAsPair = new SortableIntPair[hashes.length];
		for (int iter = 0; iter < hashes.length; iter++)
			completeHashAsPair[iter] = new SortableIntPair(hashes[iter][0], iter);

		// sort the results, sort is in place so no need to look at second
		Arrays.sort(completeHashAsPair);
		
		return storeAsArray(completeHashAsPair);
	}
	
	public OverlapInfo getFullScoreExperimental(OrderKmerHashes s, double maxShiftPercent)
	{
		int[][][] allKmerHashes = this.orderedHashes;

		// get the kmers of the second sequence
		int[][][] sAllKmerHashes = s.orderedHashes;
		
		//get sizes
		int size1 = this.size();
		int size2 = s.size();
		
		int count = 0;
		//int[] posShift = new int[Math.min(size1, size2)/8+1];
		int[] pos1Index = new int[Math.min(size1, size2)/8+1];
		int[] pos2Index = new int[pos1Index.length];
		
		//init counters
		count = 0;
		int ii1 = 0;
		int ii2 = 0;
		int i1 = 0;
		int i2 = 0;			
		
		//init the loop storage
		int hash1 = 0;
		int hash2 = 0;
		int pos1;
		int pos2;

		// perform merge operation to get the shift and the kmer count
		while (true)
		{
			//store previous value
			//int prevHash1 = hash1;
			//int prevHash2 = hash2;
			//int prevIndex1 = i1;
			//int prevIndex2 = i2;

			if (i1>=allKmerHashes[ii1].length)
			{
				ii1++;
				i1 = 0;
				
				//break if reached end
				if (ii1>=allKmerHashes.length)
					break;
			}
			if (i2>=sAllKmerHashes[ii2].length)
			{
				ii2++;
				i2 = 0;

				//break if reached end
				if (ii2>=sAllKmerHashes.length)
					break;
			}				
			
			//get the values in the array
			hash1 = allKmerHashes[ii1][i1][0];
			pos1 = allKmerHashes[ii1][i1][1];

			hash2 = sAllKmerHashes[ii2][i2][0];
			pos2 = sAllKmerHashes[ii2][i2][1];

			if (hash1 < hash2)
				i1++;
			else if (hash2 < hash1)
				i2++;
			else
			{
				//check if current shift makes sense positionally
				//int currShift = pos2-pos1;
				
				//adjust array size if needed
				if (pos1Index.length<=count)
				{
					//posShift = Arrays.copyOf(posShift, posShift.length*2);
					pos1Index = Arrays.copyOf(pos1Index, pos1Index.length*2);
					pos2Index = Arrays.copyOf(pos2Index, pos2Index.length*2);
				}
				
				// compute the shift
				//posShift[count] = currShift;								
				pos1Index[count] = pos1;
				pos2Index[count] = pos2;
				
				count++;
				i1++;
				i2++;
			}
		}

		//fit a regression line
		double eps = 1.0e-1;
		double outlier = 3.0;
		double[] r = new double[pos1Index.length];
		
		//create the copies of data
		int[] pos1GoodIndex = Arrays.copyOf(pos1Index, pos1Index.length);
		int[] pos2GoodIndex = Arrays.copyOf(pos1Index, pos1Index.length);
		
		Pair<Double,Double> linePrev = new Pair<Double,Double>(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		Pair<Double,Double> line = Utils.linearRegression(pos1GoodIndex, pos2GoodIndex, count);
		
		int validCount = 0;
  	int steps = 0;
		while (steps<10 && count>10 && (Math.abs(linePrev.x-line.x)>eps || Math.abs(linePrev.y-line.y)>eps))
		{
			//compute residuals
			for (int iter=0; iter<count; iter++)
			{
				r[iter] = (double)pos2Index[iter]-(line.x+line.y*(double)pos1Index[iter]);
				//System.err.print(" "+r[iter]);
			}
			
			double rmean = Utils.mean(r, count);
			double std = Utils.std(r, count, rmean);
			
			//go through the list and put valid regions at start
			validCount = 0;
			for (int iter=0; iter<count; iter++)
			{
				if (Math.abs(r[iter])<Math.max(200.0,outlier*std))
				{
					pos1GoodIndex[validCount] = pos1Index[iter];
					pos2GoodIndex[validCount] = pos2Index[iter];
					validCount++;
				}
			}
						
			//perform another round of regression
			linePrev = line;
			line = Utils.linearRegression(pos1GoodIndex, pos2GoodIndex, validCount);
			
			steps++;
		}
				
		//extrapolate to 0 and end point
		int a1 = Math.max(0,(int)Math.round(-line.x/line.y));
		int a2 = Math.max(0,(int)line.x.doubleValue());
		int b1 = Math.min(size(), (int)Math.round(((double)s.size()-line.x)/line.y));
		int b2 = Math.min(s.size(),(int)(line.x+line.y*(double)size()));
		
		int ahang = a1-a2;
		int bhang = (this.size()-b1>s.size()-b2) ? b1-this.size() : s.size() - b2;
		
		//compute correlation only in the proper region
		validCount = 0;
		for (int iter=0; iter<count; iter++)
		{
			if (pos2Index[iter]<a1 && pos2Index[iter]>b1)
				continue;
			
			double residual = (double)pos2Index[iter]-(line.x+line.y*(double)pos1Index[iter]);
			if (Math.abs(residual)<=400.0)
			{
				pos1GoodIndex[validCount] = pos1Index[iter];
				pos2GoodIndex[validCount] = pos2Index[iter];
				validCount++;
			}
		}

		//compute the correlation over the valid region
		double corr = 0.0;
		if (validCount>10)
			corr = Utils.pearsonCorr(pos1GoodIndex, pos2GoodIndex, validCount);

		
		//int shiftb = -medianShift - this.size() + s.size();
		//if (score>0.04)
			//System.out.format("A=%d %d, B=%d %d\n", -medianShift, a1-a2, shiftb, -medianShift+(-medianShift-(a1-a2))-this.size()+s.size());
		//	System.out.format("A=%d %d, B=%d %d\n", -medianShift, ahang, shiftb, bhang);
		//return new OverlapInfo(score, -medianShift, shiftb);
					
		//if (corr>0.06)
		//{
		//	int[] test = Arrays.copyOf(pos1Index, count);
		//	int[] test2 = Arrays.copyOf(pos2Index, count);

		//	System.err.println("corr="+corr+": ["+Arrays.toString(test)+"; "+Arrays.toString(test2)+"];");
		//}

		//the hangs are adjusted by the rate of slide*distance traveled relative to median, -medianShift-(a1-a2)
		return new OverlapInfo(corr, ahang, bhang);
	}

	
	public OverlapInfo getFullScore(OrderKmerHashes s, double maxShiftPercent)
	{
		int[][][] allKmerHashes = this.orderedHashes;

		// get the kmers of the second sequence
		int[][][] sAllKmerHashes = s.orderedHashes;
		
		//get sizes
		int size1 = this.size();
		int size2 = s.size();
		
		// init the ok regions
		int valid1Lower = 0;
		int valid1Upper = size1;
		int valid2Lower = 0;
		int valid2Upper = size2;

		int medianShift = 0;
		int overlapSize = Math.min(size1, size2);
		int absMaxShiftInOverlap = Math.max(size1, size2);

		int count = 0;
		int[] posShift = new int[Math.min(size1, size2)/8+1];
		int[] pos1Index = new int[posShift.length];
		int[] pos2Index = new int[posShift.length];
		
		//check the repeat flag
		int numScoringRepeats = 2;
		if (maxShiftPercent <= 0)
		{
			numScoringRepeats = 1;
			maxShiftPercent = Math.abs(maxShiftPercent);
		}

		// refine multiple times to get better interval estimate
		for (int repeat = 0; repeat < numScoringRepeats; repeat++)
		{
			//init counters
			count = 0;
			int ii1 = 0;
			int ii2 = 0;
			int i1 = 0;
			int i2 = 0;			
			
			//init the loop storage
			int hash1 = 0;
			int hash2 = 0;
			int pos1;
			int pos2;

			// perform merge operation to get the shift and the kmer count
			while (true)
			{
				//store previous value
				//int prevHash1 = hash1;
				//int prevHash2 = hash2;
				//int prevIndex1 = i1;
				//int prevIndex2 = i2;

				if (i1>=allKmerHashes[ii1].length)
				{
					ii1++;
					i1 = 0;
					
					//break if reached end
					if (ii1>=allKmerHashes.length)
						break;
				}
				if (i2>=sAllKmerHashes[ii2].length)
				{
					ii2++;
					i2 = 0;

					//break if reached end
					if (ii2>=sAllKmerHashes.length)
						break;
				}				
				
				//get the values in the array
				hash1 = allKmerHashes[ii1][i1][0];
				pos1 = allKmerHashes[ii1][i1][1];

				hash2 = sAllKmerHashes[ii2][i2][0];
				pos2 = sAllKmerHashes[ii2][i2][1];

				if (hash1 < hash2 || pos1 < valid1Lower || pos1 >= valid1Upper)
					i1++;
				else if (hash2 < hash1 || pos2 < valid2Lower || pos2 >= valid2Upper)
					i2++;
				else
				{
					//check if current shift makes sense positionally
					int currShift = pos2-pos1;
					if (Math.abs(currShift-medianShift) > absMaxShiftInOverlap)
					{
						//do not record this shift and increase counter
						i2++;
						continue;
					}
					
					//if (prevHash1==hash1 && i1!=prevIndex1)
					//{
					//	i1++;
					//	continue;
					//}
					//if (prevHash2==hash2 && i2!=prevIndex2)
					//{
					//	i2++;
					//	continue;
					//}


					//adjust array size if needed
					if (posShift.length<=count)
					{
						posShift = Arrays.copyOf(posShift, posShift.length*2);
						pos1Index = Arrays.copyOf(pos1Index, pos1Index.length*2);
						pos2Index = Arrays.copyOf(pos2Index, pos2Index.length*2);
					}
					
					// compute the shift
					posShift[count] = currShift;								
					pos1Index[count] = pos1;
					pos2Index[count] = pos2;
					
					count++;
					//i1++;
					i2++;
				}
			}

			// get the median
			if (count > 0)
			{
				//pick out only the best values
				if (repeat>0)
				{
					int reducedCount = -1;

					//copy over only the best values
					for (int iter=0; iter<count; iter++)
					{
						if (reducedCount>=0 && pos1Index[reducedCount]==pos1Index[iter])
						{
							//if better, record it
							if (Math.abs(posShift[reducedCount]-medianShift)>Math.abs(posShift[iter]-medianShift))
							{
								pos1Index[reducedCount] = pos1Index[iter];
								pos2Index[reducedCount] = pos2Index[iter];
								posShift[reducedCount] = posShift[iter];
							}
						}
						else
						{
							//add the new data
							reducedCount++;
							pos1Index[reducedCount] = pos1Index[iter];
							pos2Index[reducedCount] = pos2Index[iter];
							posShift[reducedCount] = posShift[iter];
						}
					}
					
					count = reducedCount+1;
				}

				medianShift = Utils.quickSelect(Arrays.copyOf(posShift, count), count / 2, count);
			}
			else
				medianShift = 0;
			
			// get the actual overlap size
			int leftPosition = Math.max(0, medianShift);
			int rightPosition = Math.min(size2, size1 + medianShift);
			overlapSize = Math.max(50,rightPosition - leftPosition);

			//compute the max possible allowed shift in kmers
			absMaxShiftInOverlap = Math.min(Math.max(size1, size2), (int)((double)overlapSize*maxShiftPercent));

			// get the updated borders
			valid1Lower = Math.max(0, -medianShift - absMaxShiftInOverlap);
			valid1Upper = Math.min(size1, size2 - medianShift + absMaxShiftInOverlap);
			valid2Lower = Math.max(0, medianShift - absMaxShiftInOverlap);
			valid2Upper = Math.min(size2, size1 + medianShift + absMaxShiftInOverlap);

			/*
			System.err.println(overlapSize);
			System.err.println("Size1= "+size1+" Lower:"+
			valid1Lower+" Upper:"+valid1Upper+" Shift="+shift);
			System.err.println("Size2= "+size2+" Lower:"+
			valid2Lower+" Upper:"+valid2Upper);
			*/
		}
			
		//storage for edge computation
		int leftEdge1 = Integer.MAX_VALUE;
		int leftEdge2 = Integer.MAX_VALUE;
		int rightEdge1 = Integer.MIN_VALUE;
		int rightEdge2 = Integer.MIN_VALUE;

		int validCount = 0;
		for (int iter=0; iter<count; iter++)
		{
			int pos1 = pos1Index[iter];
			int pos2 = pos2Index[iter];
			
			//take only valid values
			if (Math.abs(posShift[iter] - medianShift) > absMaxShiftInOverlap)
				continue;

			//get the edges
			if (pos1<leftEdge1)
				leftEdge1 = pos1;
			if (pos2<leftEdge2)
				leftEdge2 = pos2;
			if (pos1>rightEdge1)
				rightEdge1 = pos1;
			if (pos2>rightEdge2)
				rightEdge2 = pos2;

			validCount++;
		}

		//compute the score
		double score = (double) validCount / (double) (overlapSize);

		//get edge info based on the german tank problem
		int gap1 = (int)Math.round((rightEdge1-validCount)/(double)validCount);
		int gap2 = (int)Math.round((rightEdge2-validCount)/(double)validCount);
		int a1 = Math.max(0,leftEdge1-gap1);
		int a2 = Math.max(0,leftEdge2-gap2);
		int b1 = Math.min(this.size(),rightEdge1+gap1);
		int b2 = Math.min(s.size(),rightEdge2+gap2);
		
		int ahang = a1-a2;
		int bhang = (this.size()-b1>s.size()-b2) ? b1-this.size() : s.size() - b2;
		
		//int shiftb = -medianShift - this.size() + s.size();
		//if (score>0.04)
			//System.out.format("A=%d %d, B=%d %d\n", -medianShift, a1-a2, shiftb, -medianShift+(-medianShift-(a1-a2))-this.size()+s.size());
		//	System.out.format("A=%d %d, B=%d %d\n", -medianShift, ahang, shiftb, bhang);
		//return new OverlapInfo(score, -medianShift, shiftb);
					
		//if (score>0.06)
		//{
		//	int[] test = Arrays.copyOf(posShift, count);
		//	int[] test2 = Arrays.copyOf(pos1Index, count);

		//	System.err.println("Start = "+Math.max(0, -medianShift)+", Overlap="+overlapSize+" Maxshift="+absMaxShiftInOverlap+": ["+Arrays.toString(test)+"; "+Arrays.toString(test2)+"];");
		//	System.err.println("Overlap="+overlapSize+", Shift/overlap="+(double)(test[test.length-10]-test[10])/(double)overlapSize);
		//}

		//the hangs are adjusted by the rate of slide*distance traveled relative to median, -medianShift-(a1-a2)
		return new OverlapInfo(score, ahang, bhang);
	}

	public int size()
	{
		//make sure you get the last element even if its different size
		return (this.orderedHashes.length-1)*this.orderedHashes[0].length+this.orderedHashes[this.orderedHashes.length-1].length;
	}
}
