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
package edu.umd.marbl.mhap.sketch;

import it.unimi.dsi.fastutil.ints.IntArrays;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.Arrays;

import edu.umd.marbl.mhap.impl.OverlapInfo;
import edu.umd.marbl.mhap.utils.Utils;

public final class OrderedNGramHashes
{
	private final static class EdgeData
	{
		public final int a1;
		public final int a2;
		public final int b1;
		public final int b2;
		public final int count;
		
		public EdgeData(int a1, int a2, int b1, int b2, int count)
		{
			this.a1 = a1;
			this.a2 = a2;
			this.b1 = b1;
			this.b2 = b2;
			this.count = count;
		}
	}
	
	private final static class MatchData
	{
		private int absMaxShiftInOverlap;
		private int count; 
		private final double maxShiftPercent;
		private int medianShift;
		private boolean needRecompute;
		public int[] pos1Index;
		public int[] pos2Index;
		public int[] posShift;
		private final int seqLength1;
		private final int seqLength2;

		public MatchData(OrderedNGramHashes o1, OrderedNGramHashes o2, double maxShiftPercent)
		{
			this.seqLength1 = o1.getSequenceLength();
			this.seqLength2 = o2.getSequenceLength();
			
			this.posShift = new int[Math.max(o1.size(), o2.size())/4+1];
			this.pos1Index = new int[posShift.length];
			this.pos2Index = new int[posShift.length];
			
			this.maxShiftPercent = maxShiftPercent;
			reset();
		}
		
		public EdgeData computeEdges()
		{
			// storage for edge computation
			int leftEdge1 = Integer.MAX_VALUE;
			int leftEdge2 = Integer.MAX_VALUE;
			int rightEdge1 = Integer.MIN_VALUE;
			int rightEdge2 = Integer.MIN_VALUE;

			// count only the shifts in the correct place
			int validCount = 0;
			int medianShift = getMedianShift();
			int absMaxShiftInOverlap = getAbsMaxShift();
			int count = size();
			for (int iter = 0; iter < count; iter++)
			{
				int pos1 = this.pos1Index[iter];
				int pos2 = this.pos2Index[iter];

				// take only valid values
				if (Math.abs(this.posShift[iter] - medianShift) > absMaxShiftInOverlap)
					continue;

				// get the edges
				if (pos1 < leftEdge1)
					leftEdge1 = pos1;
				if (pos2 < leftEdge2)
					leftEdge2 = pos2;
				if (pos1 > rightEdge1)
					rightEdge1 = pos1;
				if (pos2 > rightEdge2)
					rightEdge2 = pos2;

				validCount++;
			}

			if (validCount < 3)
				return null;

			// get edge info uniformly minimum variance unbiased (UMVU) estimators
			// a = (n*a-b)/(n-1)
			// b = (n*b-a)/(n-1)
			int a1 = Math.max(0, (int) Math.round((double)(validCount * leftEdge1 - rightEdge1) / (double) (validCount - 1)));
			int a2 = Math.min(this.seqLength1, (int) Math.round((double)(validCount * rightEdge1 - leftEdge1) / (double) (validCount - 1)));
			int b1 = Math.max(0, (int) Math.round((double)(validCount * leftEdge2 - rightEdge2) / (double) (validCount - 1)));
			int b2 = Math.min(this.seqLength2, (int) Math.round((double)(validCount * rightEdge2 - leftEdge2) / (double) (validCount - 1)));
			
			return new EdgeData(a1, a2, b1, b2, validCount);
		}
		
		public int getAbsMaxShift()
		{
			performUpdate();
			return this.absMaxShiftInOverlap;
		}
		
		public int getMedianShift()
		{
			performUpdate();
			return this.medianShift;
		}
		
		public boolean isEmpty()
		{
			return this.count<=0;
		}
		
		public void optimizeShifts()
		{
			if (isEmpty())
				return;
			
			int reducedCount = -1;

			// copy over only the best values
			int medianShift = getMedianShift();
			for (int iter = 0; iter < this.count; iter++)
			{
				if (reducedCount >= 0 && pos1Index[reducedCount] == pos1Index[iter])
				{
					// if better, record it
					if (Math.abs(posShift[reducedCount] - medianShift) > Math.abs(posShift[iter] - medianShift))
					{
						pos1Index[reducedCount] = pos1Index[iter];
						pos2Index[reducedCount] = pos2Index[iter];
						posShift[reducedCount] = posShift[iter];
					}
				}
				else
				{
					// add the new data
					reducedCount++;
					pos1Index[reducedCount] = pos1Index[iter];
					pos2Index[reducedCount] = pos2Index[iter];
					posShift[reducedCount] = posShift[iter];
				}
			}

			this.count = reducedCount + 1;			
			this.needRecompute = true;
		}
		
		private void performUpdate()
		{
			if (this.needRecompute)
			{
				if (this.count>0)
				{
					this.medianShift = Utils.quickSelect(Arrays.copyOf(this.posShift, this.count), this.count / 2, this.count);
					
					// get the actual overlap size
					int leftPosition = Math.max(0, -this.medianShift);
					int rightPosition = Math.min(this.seqLength1, this.seqLength2 - this.medianShift);
					int overlapSize = Math.max(10, rightPosition - leftPosition);

					// compute the max possible allowed shift in kmers
					this.absMaxShiftInOverlap = Math.min(Math.max(this.seqLength1, this.seqLength2), (int) ((double) overlapSize * maxShiftPercent));
				}
				else
				{
					this.medianShift = 0;
					this.absMaxShiftInOverlap = Math.max(this.seqLength1, this.seqLength2);
				}
			}
			
			this.needRecompute = false;
		}
		
		public void recordMatch(int pos1, int pos2, int shift)
		{
			// adjust array size if needed
			if (posShift.length <= this.count)
			{
				this.posShift = Arrays.copyOf(this.posShift, this.posShift.length * 2);
				this.pos1Index = Arrays.copyOf(this.pos1Index, this.pos1Index.length * 2);
				this.pos2Index = Arrays.copyOf(this.pos2Index, this.pos2Index.length * 2);
			}
			
			posShift[this.count] = shift;
			pos1Index[this.count] = pos1;
			pos2Index[this.count] = pos2;
			
			this.count++;
			this.needRecompute = true;
		}
		
		public void reset()
		{
			this.count = 0;
			this.needRecompute = true;
		}

		public int size()
		{
			return this.count;
		}

		public int valid1Lower()
		{
			performUpdate();
			int valid = Math.max(0, -getMedianShift() - getAbsMaxShift());
			
			return valid;
		}
		
		public int valid1Upper()
		{
			performUpdate();
			int valid = Math.min(this.seqLength1, this.seqLength2 - getMedianShift() + getAbsMaxShift());
			
			return valid;
		}

		public int valid2Lower()
		{
			performUpdate();
			int valid = Math.max(0, getMedianShift() - getAbsMaxShift());
			
			return valid;
		}
		
		public int valid2Upper()
		{
			performUpdate();
			int valid = Math.min(this.seqLength2, this.seqLength1 + getMedianShift() + getAbsMaxShift());
			
			return valid;
		}
	}
	
	private final int[][] orderedHashes;
	private final int seqLength;

	private static double computeKBottomSketchJaccard(int[][] seq1Hashes, int[][] seq2Hashes, int medianShift, int absMaxShiftInOverlap, int a1, int a2, int b1, int b2)
	{
		//get k for first string
		int s1 = 0;
		int[][] array1 = new int[seq1Hashes.length][];
		for (int i=0; i<seq1Hashes.length; i++)
		{
			int pos = seq1Hashes[i][1];
			if (pos >= a1 && pos <= a2)
			{
				array1[s1] = seq1Hashes[i];
				s1++;
			}
		}
		
		//get k for second string
		int s2 = 0;
		int[][] array2 = new int[seq2Hashes.length][];
		for (int j=0; j<seq2Hashes.length; j++)
		{
			int pos = seq2Hashes[j][1];
			if (pos >= b1 && pos <= b2)
			{
				array2[s2] = seq2Hashes[j];
				s2++;
			}
		}
		
		//compute k
		int k = Math.min(s1,s2);
		
		//empty has jaccard of 1
		if (k==0)
			return 0;
			
		//perform the k-bottom count
		int i = 0;
		int j = 0;
		int intersectCount = 0;
		int unionCount = 0;
		while (unionCount<k)
		{
			if (array1[i][0]<array2[j][0])
				i++;
			else
			if (array1[i][0]>array2[j][0])
				j++;
			else
			{
				intersectCount++;
				i++;
				j++;
			}
			
			unionCount++;
		}
		
		double score = ((double)intersectCount)/(double)k;
		
		return score;
	}

	public final static OrderedNGramHashes fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			// dos.writeInt(this.seqLength);
			// dos.writeInt(size());
			int seqLength = input.readInt();
			int hashLength = input.readInt();

			int[][] orderedHashes = new int[hashLength][2];

			for (int iter = 0; iter < hashLength; iter++)
			{
				// dos.writeInt(this.completeHash[iter][iter2]);
				orderedHashes[iter][0] = input.readInt();
				orderedHashes[iter][1] = input.readInt();
			}

			return new OrderedNGramHashes(seqLength, orderedHashes);

		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	public static double jaccardToIdentity(double score, int gramSize)
	{
		return -1.0/(double)gramSize*Math.log(2.0*score/(1.0+score));
	}

	private static void recordMatchingKmers(
			MatchData matchData, 
			int[][] seq1KmerHashes, 
			int[][] seq2KmerHashes,
			int repeat)
	{
		// init the loop storage
		int hash1;
		int hash2;
		int pos1;
		int pos2;
		
		// init the borders
		int medianShift = matchData.getMedianShift();
		int absMaxShift = matchData.getAbsMaxShift();
		int valid1Lower = matchData.valid1Lower();
		int valid2Lower = matchData.valid2Lower();
		int valid1Upper = matchData.valid1Upper();
		int valid2Upper = matchData.valid2Upper();
		
		// init counters
		int i1 = 0;
		int i2 = 0;
		
		//reset the data, redo the shifts
		matchData.reset();
		
		// perform merge operation to get the shift and the kmer count
		while (true)
		{
			if (i1>=seq1KmerHashes.length)
				break;
			if (i2>=seq2KmerHashes.length)
				break;
			
			// get the values in the array
			hash1 = seq1KmerHashes[i1][0];
			pos1 = seq1KmerHashes[i1][1];
			hash2 = seq2KmerHashes[i2][0];
			pos2 = seq2KmerHashes[i2][1];

			if (hash1 < hash2 || pos1 < valid1Lower || pos1 >= valid1Upper)
				i1++;
			else if (hash2 < hash1 || pos2 < valid2Lower || pos2 >= valid2Upper)
				i2++;
			else
			{
				// check if current shift makes sense positionally
				int currShift = pos2 - pos1;
				int diffFromExpected = currShift - medianShift;
				if (diffFromExpected > absMaxShift)
					i1++;
				else
				if (diffFromExpected < -absMaxShift)
					i2++;
				else
				{				
					//record match
					matchData.recordMatch(pos1, pos2, currShift);
	
					// don't rely on repeats in the first iteration
					if (repeat == 0)
						i1++;
					i2++;
				}
			}
		}
	}

	private OrderedNGramHashes(int seqLength, int[][] orderedHashes)
	{
		this.seqLength = seqLength;
		this.orderedHashes = orderedHashes;
	}

	public OrderedNGramHashes(String seq, int kmerSize, int sketchSize)
	{
		this.seqLength = seq.length() - kmerSize + 1;
		
		if (this.seqLength<=0)
			throw new SketchRuntimeException("Sequence length must be greater or equal to kmerSize.");
		
		// compute just direct hash of sequence
		int[] hashes = HashUtils.computeSequenceHashes(seq, kmerSize);

		int[] perm = new int[hashes.length];

		//init the array
		for (int iter = 0; iter < hashes.length; iter++)
			perm[iter] = iter;
		
		//sort the array
		IntArrays.radixSortIndirect(perm, hashes, true);
		
		//sketchSize = (int)Math.round(0.25*(double)this.seqLength);

		//find the largest storage value
		int k = Math.min(sketchSize, hashes.length);
		
		//allocate the memory
		this.orderedHashes = new int[k][2];

		for (int iter = 0; iter < this.orderedHashes.length; iter++)
		{
			int index = perm[iter];
			this.orderedHashes[iter][0] = hashes[index];
			this.orderedHashes[iter][1] = index;
		}
	}

	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream(size() * 2);
		DataOutputStream dos = new DataOutputStream(bos);

		try
		{
			dos.writeInt(this.seqLength);
			dos.writeInt(size());
			for (int iter = 0; iter < this.orderedHashes.length; iter++)
			{
				dos.writeInt(this.orderedHashes[iter][0]);
				dos.writeInt(this.orderedHashes[iter][1]);
			}

			dos.flush();
			return bos.toByteArray();
		}
		catch (IOException e)
		{
			throw new SketchRuntimeException("Unexpected IO error.", e);
		}
	}
	
	public int getHash(int index)
	{
		return this.orderedHashes[index][0];
	}
	
	public OverlapInfo getOverlapInfo(OrderedNGramHashes toSequence, double maxShiftPercent)
	{
		//allocate the memory for the search
		MatchData matchData = new MatchData(this, toSequence, maxShiftPercent);

		//get the initial matches
		recordMatchingKmers(matchData, this.orderedHashes, toSequence.orderedHashes, 0);			
		if (matchData.isEmpty())
			return OverlapInfo.EMPTY;

		//get matches again, but now in a better region
		recordMatchingKmers(matchData, this.orderedHashes, toSequence.orderedHashes, 1);

		if (matchData.isEmpty())
			return OverlapInfo.EMPTY;

		matchData.optimizeShifts();
			
		if (matchData.isEmpty())
			return OverlapInfo.EMPTY;

		//get the edge data
		EdgeData edgeData = matchData.computeEdges();
		
		if (edgeData==null)
			return OverlapInfo.EMPTY;
		
		//compute the jaccard score using bottom-k sketching
		double score = computeKBottomSketchJaccard(this.orderedHashes, toSequence.orderedHashes, matchData.getMedianShift(), matchData.getAbsMaxShift(), edgeData.a1, edgeData.a2, edgeData.b1, edgeData.b2);
		
		double rawScore = (double)edgeData.count;
		//double rawScore = jaccardToIdentity(score, kmerSize);

		return new OverlapInfo(score, rawScore, edgeData.a1, edgeData.a2, edgeData.b1, edgeData.b2);
	}
	
	public int getSequenceLength()
	{
		return this.seqLength;
	}

	public int size()
	{
		return this.orderedHashes.length;
	}
}
