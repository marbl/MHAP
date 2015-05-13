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
import java.io.Serializable;
import java.util.Arrays;

import edu.umd.marbl.mhap.impl.OverlapInfo;
import edu.umd.marbl.mhap.utils.Utils;

public final class OrderedNGramHashes
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
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#toString()
		 */
		@Override
		public String toString()
		{
			return "["+x + ", " + y + "]";
		}
	
	}

	private final int[][] orderedHashes;
	private final int seqLength;

	public final static int REDUCTION = 4;

	private final static int[][] allocateMemory(int size)
	{
		// allocate the memory
		int[][] completeHash = new int[size][2];

		return completeHash;
	}

	public final static OrderedNGramHashes fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			// dos.writeInt(this.seqLength);
			// dos.writeInt(size());
			int seqLength = input.readInt();
			int hashLength = input.readInt();

			int[][] orderedHashes = allocateMemory(hashLength);

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

	private OrderedNGramHashes(int seqLength, int[][] orderedHashes)
	{
		this.seqLength = seqLength;
		this.orderedHashes = orderedHashes;
	}

	public OrderedNGramHashes(String seq, int kmerSize)
	{
		this.seqLength = seq.length() - kmerSize + 1;
		
		if (this.seqLength<=0)
			throw new SketchRuntimeException("Sequence length must be greater or equal to kmerSize.");
		
		this.orderedHashes = getFullHashes(seq, kmerSize);
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
			throw new SketchRuntimeException("Unexpected IO error.");
		}
	}

	public int getHash(int index)
	{
		return this.orderedHashes[index][0];
	}

	private int[][] storeAsArray(SortableIntPair[] completeHashAsPair)
	{
		// allocate the memory
		int[][] completeHash = allocateMemory(completeHashAsPair.length);

		for (int iter = 0; iter < completeHashAsPair.length; iter++)
		{
			completeHash[iter][0] = completeHashAsPair[iter].x;
			completeHash[iter][1] = completeHashAsPair[iter].y;
		}

		return completeHash;
	}

	private int[][] getFullHashes(String seq, int subKmerSize)
	{
		int cutoff = (int) ((long) Integer.MIN_VALUE + ((long) Integer.MAX_VALUE - (long) Integer.MIN_VALUE)
				/ (long) REDUCTION);

		// compute just direct hash of sequence
		int[] hashes = HashUtils.computeSequenceHashes(seq, subKmerSize);

		int count = 0;
		for (int val : hashes)
			if (val <= cutoff)
				count++;

		int[] cutHashes = new int[count];
		int[] perm = new int[count];
		int[] pos = new int[count];
		
		count = 0;
		for (int iter = 0; iter < hashes.length; iter++)
			if (hashes[iter] <= cutoff)
			{
				cutHashes[count] = hashes[iter];
				perm[count] = count;
				pos[count] = iter;

				count++;
			}
		
		//sort the array
		IntArrays.radixSortIndirect(perm, cutHashes, true);

		SortableIntPair[] completeHashAsPair = new SortableIntPair[count];		
		for (int iter=0; iter<count; iter++)
		{
			int index = perm[iter];
			completeHashAsPair[iter] = new SortableIntPair(cutHashes[index], pos[index]);
		}
		
		//System.err.println(Arrays.toString(completeHashAsPair));
		// sort the results, sort in place so no need to look at second
		//Arrays.sort(completeHashAsPair);

		return storeAsArray(completeHashAsPair);
	}

	public OverlapInfo getOverlapInfo(OrderedNGramHashes s, double maxShiftPercent)
	{
		int[][] allKmerHashes = this.orderedHashes;

		// get the kmers of the second sequence
		int[][] sAllKmerHashes = s.orderedHashes;

		// get sizes
		int size1 = this.size();
		int size2 = s.size();

		int kmerSize1 = this.seqLength;
		int kmerSize2 = s.seqLength;

		// init the ok regions
		int valid1Lower = 0;
		int valid1Upper = kmerSize1;
		int valid2Lower = 0;
		int valid2Upper = kmerSize2;

		int medianShift = 0;
		int overlapSize = Math.min(kmerSize1, kmerSize2);
		int absMaxShiftInOverlap = Math.max(kmerSize1, kmerSize2);

		int count = 0;
		int[] posShift = new int[Math.min(size1, size2) / 8 + 1];
		int[] pos1Index = new int[posShift.length];
		int[] pos2Index = new int[posShift.length];

		// check the repeat flag
		int numScoringRepeats = 2;
		if (maxShiftPercent <= 0)
		{
			numScoringRepeats = 1;
			maxShiftPercent = Math.abs(maxShiftPercent);
		}

		// refine multiple times to get better interval estimate
		for (int repeat = 0; repeat < numScoringRepeats; repeat++)
		{
			// init counters
			count = 0;
			int i1 = 0;
			int i2 = 0;

			// init the loop storage
			int hash1 = 0;
			int hash2 = 0;
			int pos1;
			int pos2;

			// perform merge operation to get the shift and the kmer count
			while (true)
			{
				if (i1>=allKmerHashes.length)
					break;
				if (i2>=sAllKmerHashes.length)
					break;
				
				// get the values in the array
				hash1 = allKmerHashes[i1][0];
				pos1 = allKmerHashes[i1][1];

				hash2 = sAllKmerHashes[i2][0];
				pos2 = sAllKmerHashes[i2][1];

				if (hash1 < hash2 || pos1 < valid1Lower || pos1 >= valid1Upper)
					i1++;
				else if (hash2 < hash1 || pos2 < valid2Lower || pos2 >= valid2Upper)
					i2++;
				else
				{
					// check if current shift makes sense positionally
					int currShift = pos2 - pos1;
					if (Math.abs(currShift - medianShift) > absMaxShiftInOverlap)
					{
						// do not record this shift and increase counter
						i2++;
						continue;
					}

					// adjust array size if needed
					if (posShift.length <= count)
					{
						posShift = Arrays.copyOf(posShift, posShift.length * 2);
						pos1Index = Arrays.copyOf(pos1Index, pos1Index.length * 2);
						pos2Index = Arrays.copyOf(pos2Index, pos2Index.length * 2);
					}

					// compute the shift
					posShift[count] = currShift;
					pos1Index[count] = pos1;
					pos2Index[count] = pos2;

					// if first round, store only first hit
					if (repeat == 0)
						i1++;
					i2++;

					count++;
				}
			}

			if (count <= 0)
				return new OverlapInfo(0.0, 0, 0, 0, 0, 0);

			// pick out only the matches that are best
			if (repeat > 0)
			{
				int reducedCount = -1;

				// copy over only the best values
				for (int iter = 0; iter < count; iter++)
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

				count = reducedCount + 1;
			}

			if (count <= 0)
				medianShift = 0;
			else
				medianShift = Utils.quickSelect(Arrays.copyOf(posShift, count), count / 2, count);

			// get the actual overlap size
			int leftPosition = Math.max(0, -medianShift);
			int rightPosition = Math.min(kmerSize1, kmerSize2 - medianShift);
			overlapSize = Math.max(this.seqLength - kmerSize1, rightPosition - leftPosition);

			// compute the max possible allowed shift in kmers
			absMaxShiftInOverlap = Math.min(Math.max(kmerSize1, kmerSize2),
					(int) ((double) overlapSize * maxShiftPercent));

			// get the updated borders
			valid1Lower = Math.max(0, -medianShift - absMaxShiftInOverlap);
			valid1Upper = Math.min(kmerSize1, kmerSize2 - medianShift + absMaxShiftInOverlap);
			valid2Lower = Math.max(0, medianShift - absMaxShiftInOverlap);
			valid2Upper = Math.min(kmerSize2, kmerSize1 + medianShift + absMaxShiftInOverlap);

			/*
			 * System.err.println(overlapSize);
			 * System.err.println("Size1= "+size1+" Lower:"+
			 * valid1Lower+" Upper:"+valid1Upper+" Shift="+shift);
			 * System.err.println("Size2= "+size2+" Lower:"+
			 * valid2Lower+" Upper:"+valid2Upper);
			 */
		}

		// storage for edge computation
		int leftEdge1 = Integer.MAX_VALUE;
		int leftEdge2 = Integer.MAX_VALUE;
		int rightEdge1 = Integer.MIN_VALUE;
		int rightEdge2 = Integer.MIN_VALUE;

		// count only the shifts in the correct place
		int validCount = 0;
		for (int iter = 0; iter < count; iter++)
		{
			int pos1 = pos1Index[iter];
			int pos2 = pos2Index[iter];

			// take only valid values
			if (Math.abs(posShift[iter] - medianShift) > absMaxShiftInOverlap)
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

		if (validCount <= 1)
			return new OverlapInfo(0.0, 0, 0, 0, 0, 0);

		// compute the score
		double score = (double) validCount / (double) (overlapSize);

		// get edge info uniformly minimum variance unbiased (UMVU) estimators
		// a = (n*a-b)/(n-1)
		// b = (n*b-a)/(n-1)
		int a1 = Math.max(0, (int) Math.round((validCount * leftEdge1 - rightEdge1) / (double) (validCount - 1)));
		int b1 = Math.max(0, (int) Math.round((validCount * leftEdge2 - rightEdge2) / (double) (validCount - 1)));
		int a2 = Math.min(this.seqLength,
				(int) Math.round((validCount * rightEdge1 - leftEdge1) / (double) (validCount - 1)));
		int b2 = Math.min(s.seqLength,
				(int) Math.round((validCount * rightEdge2 - leftEdge2) / (double) (validCount - 1)));

		// int ahang = a1-a2;
		// int bhang = (this.size()-b1>s.size()-b2) ? b1-this.size() : s.size()
		// - b2;

		// if (score>0.06)
		// {
		// int[] test = Arrays.copyOf(posShift, count);
		// int[] test2 = Arrays.copyOf(pos1Index, count);

		// System.err.println("Start = "+Math.max(0,
		// -medianShift)+", Overlap="+overlapSize+" Maxshift="+absMaxShiftInOverlap+": ["+Arrays.toString(test)+"; "+Arrays.toString(test2)+"];");
		// System.err.println("Overlap="+overlapSize+", Shift/overlap="+(double)(test[test.length-10]-test[10])/(double)overlapSize);
		// }

		// the hangs are adjusted by the rate of slide*distance traveled
		// relative to median, -medianShift-(a1-a2)
		// return new OverlapInfo(score, ahang, bhang);

		return new OverlapInfo(score * (double) REDUCTION, validCount, a1, a2, b1, b2);
	}

	public int size()
	{
		return this.orderedHashes.length;
	}
}
