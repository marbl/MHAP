package com.secret.fastalign.general;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;
import com.secret.fastalign.utils.Utils;

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

	public final static double SHIFT_CONSENSUS_PERCENTAGE = 0.50;
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
	
	public Pair<Double, Integer> getFullScore(OrderKmerHashes s, double maxShiftPercent)
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

			// perform merge operation to get the shift and the kmer count
			while (true)
			{
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
				int hash1 = allKmerHashes[ii1][i1][0];
				int pos1 = allKmerHashes[ii1][i1][1];

				int hash2 = sAllKmerHashes[ii2][i2][0];
				int pos2 = sAllKmerHashes[ii2][i2][1];

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
					}

					//adjust array size if needed
					if (posShift.length<=count)
						posShift = Arrays.copyOf(posShift, posShift.length*2);
					
					// compute the shift
					posShift[count] = currShift;					
					
					//if (repeat>0)
					//	System.err.println(""+pos1+" "+pos2);

					count++;
					i1++;
					i2++;
				}
			}

			// get the median
			if (count > 0)
			{
				medianShift = Utils.quickSelect(posShift, count / 2, count);
			}
			else
				medianShift = 0;

			/*
			int[] test = Arrays.copyOf(posShift, count);
			Arrays.sort(test);
			System.err.println(Arrays.toString(test));
			*/

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
		
		// count percent valid shift, there must be a consensus
		int validCount = 0;
		for (int iter = 0; iter < count; iter++)
		{
			if (Math.abs(posShift[iter] - medianShift) <= absMaxShiftInOverlap)
				validCount++;
		}
		double validShiftPercent = (double) validCount / (double) count;
		
		//System.err.println(overlapSize);

		double score = 0;
		if (overlapSize > 0 && validShiftPercent > SHIFT_CONSENSUS_PERCENTAGE)
			score = (double) count / (double) (overlapSize);

		return new Pair<Double, Integer>(score, medianShift);
	}

	public int size()
	{
		//make sure you get the last element even if its different size
		return (this.orderedHashes.length-1)*this.orderedHashes[0].length+this.orderedHashes[this.orderedHashes.length-1].length;
	}
}
