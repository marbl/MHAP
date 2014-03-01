package com.secret.fastalign.minhash;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;
import com.secret.fastalign.utils.Utils;

public final class SequenceMinHashes 
{
	private static final class SortableIntPair implements Comparable<SortableIntPair>
	{
		public final int x;
		public final int y;
				
		public SortableIntPair(int x, int y)
		{
			this.x = x;
			this.y = y;
		}


		@Override
		public int compareTo(SortableIntPair p)
		{
			return Integer.compare(this.x, p.x);
			
			//int result = Integer.compare(this.x, p.x);
			
			//if (result!=0)
			//	return result;
			
			//return Integer.compare(this.y, p.y);
		}		
	}
	
	private final int[][] completeHash;
	//private final int subKmerSize;
	//private final Sequence seq;
	
	private final SequenceId id;
	private final MinHash mainHashes;
	
	public final static double SHIFT_CONSENSUS_PERCENTAGE = 0.75;
	
	public static SequenceMinHashes fromByteStream(DataInputStream input, int offset) throws IOException
	{
		try
		{
			//dos.writeInt(this.id.getHeaderId());
			//dos.writeBoolean(this.id.isForward());
			SequenceId id = new SequenceId(input.readInt()+offset, input.readBoolean());
			
			//dos.write(this.mainHashes.getAsByteArray());
			MinHash mainHashes = MinHash.fromByteStream(input);
			
			if (mainHashes==null)
				throw new FastAlignRuntimeException("Unexpected hash read error.");
			
			//dos.writeInt(this.completeHash.length);
			int hashLength = input.readInt();			
			
			int[][] completeHash = new int[hashLength][]; 
			for (int iter=0; iter<hashLength; iter++)
			{
				//dos.writeInt(this.completeHash[iter][iter2]);
				completeHash[iter] = new int[2];
				completeHash[iter][0] = input.readInt();
				completeHash[iter][1] = input.readInt();					
			}
			
			return new SequenceMinHashes(id, mainHashes, completeHash);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numHashes, int subSequenceSize, int subKmerSize, boolean storeHashes, HashSet<Integer> filter)
	{
		this.id = seq.getId();
		this.mainHashes = new MinHash(seq, kmerSize, numHashes, subSequenceSize, filter);
			
		if (storeHashes)
		{
			this.completeHash = getFullHashes(seq, subKmerSize);
			//this.seq = null;
		}
		else
		{
			this.completeHash = getFullHashes(seq, subKmerSize);
			//this.seq = null;
		}
	}
	
	private SequenceMinHashes(SequenceId id, MinHash mainHashes, int[][] completeHash)
	{
		this.id = id;
		this.mainHashes = mainHashes;
		this.completeHash = completeHash;
	}
	
	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
    DataOutputStream dos = new DataOutputStream(bos);
    
    try
		{
      dos.writeInt(this.id.getHeaderId());
      dos.writeBoolean(this.id.isForward());
			dos.write(this.mainHashes.getAsByteArray());
			
			dos.writeInt(this.completeHash.length);
	    for (int iter=0; iter<this.completeHash.length; iter++)
	    	for (int iter2=0; iter2<2; iter2++)
	    		dos.writeInt(this.completeHash[iter][iter2]);
    
	    dos.flush();
	    return bos.toByteArray();
		}
    catch (IOException e)
    {
    	throw new FastAlignRuntimeException("Unexpected IO error.");
    }
	}
	
	public int[][] getFullHashes()
	{
		//if (this.completeHash!=null)
		return this.completeHash;
		
		//return getFullHashes(this.seq, this.subKmerSize);
	}
	
	private int[][] getFullHashes(Sequence seq, int subKmerSize)
	{
		//compute just direct hash of sequence
		int[][] hashes = Utils.computeKmerHashesInt(seq, subKmerSize, 1);	
		
		SortableIntPair[] completeHashAsPair = new SortableIntPair[hashes.length];	
		for (int iter=0; iter<hashes.length; iter++)
			completeHashAsPair[iter] = new SortableIntPair(hashes[iter][0],iter);
		
		//sort the results, sort is in place so no need to look at second
		Arrays.sort(completeHashAsPair);
		
		//store in array to reduce memory
		int[][] completeHash = new int[completeHashAsPair.length][];
		for (int iter=0; iter<completeHashAsPair.length; iter++)
		{
			completeHash[iter] = new int[2];
			completeHash[iter][0] = completeHashAsPair[iter].x;
			completeHash[iter][1] = completeHashAsPair[iter].y;
		}		
		
		return completeHash;
	}
	
	public Pair<Double, Integer> getFullScore(int[][] allKmerHashes, SequenceMinHashes s, int maxShift)
	{	
		if (allKmerHashes==null)
			throw new FastAlignRuntimeException("Hash input cannot be null.");
		
		//get the kmers of the second sequence
		int[][] sAllKmerHashes = s.getFullHashes();
		
		//init the ok regions
		int valid1Lower = 0;
		int valid1Upper = allKmerHashes.length;
		int valid2Lower = 0;
		int valid2Upper = sAllKmerHashes.length;
		int overlapSize = 0;
		int border = maxShift;
		
		int count = 0;
		int shift = 0;
		int[] posShift = new int[Math.min(allKmerHashes.length, sAllKmerHashes.length)];
		
		int numScoringRepeats = 2;
		if (maxShift<=0)
			numScoringRepeats = 1;
		
		//make it positive
		maxShift = Math.abs(maxShift);
		
		//refine multiple times to get better interval estimate
		for (int repeat=0; repeat<numScoringRepeats; repeat++)
		{
			count = 0;
			int iter1 = 0;
			int iter2 = 0;
			
			//perform merge operation to get the shift and the kmer count
			while (iter1<allKmerHashes.length && iter2<sAllKmerHashes.length)
			{
				int[] s1 = allKmerHashes[iter1];
				int[] s2 = sAllKmerHashes[iter2];
				
				if (s1[0] < s2[0] || s1[1]<valid1Lower || s1[1]>=valid1Upper)
					iter1++;
				else
				if (s2[0] < s1[0] || s2[1]<valid2Lower || s2[1]>=valid2Upper)
					iter2++;
				else
				{
					//compute the shift
					posShift[count] = s2[1]-s1[1];
					
					count++;
					iter1++;
					iter2++;
				}
			}
			
			//get the median
			if (count>0)
			{
				shift = Utils.quickSelect(posShift, count/2, count);
			}
			else
				shift = 0;
			
			//int[] test = Arrays.copyOf(posShift, count);	
			//Arrays.sort(test);
			//System.err.println(Arrays.toString(test));
			
			//get the updated borders
			valid1Lower = Math.max(0, -shift-border);
			valid1Upper = Math.min(allKmerHashes.length, sAllKmerHashes.length-shift+border);
			valid2Lower = Math.max(0, shift-border);
			valid2Upper = Math.min(sAllKmerHashes.length, allKmerHashes.length+shift+border);

			//get the actual overlap size
			int valid2LowerBorder = Math.max(0, shift);
			int valid2UpperBorder = Math.min(sAllKmerHashes.length, allKmerHashes.length+shift);
			overlapSize = valid2UpperBorder-valid2LowerBorder;
			
			//System.err.println(overlapSize);
			//System.err.println("Size1= "+allKmerHashes.length+" Lower:"+ valid1Lower+" Upper:"+valid1Upper+" Shift="+shift);
			//System.err.println("Size2= "+sAllKmerHashes.length+" Lower:"+ valid2Lower+" Upper:"+valid2Upper);			
		}
		
		//count percent valid shift, there must be a consensus 
		int validCount = 0;
		for (int iter=0; iter<count; iter++)
		{
			if (Math.abs(posShift[iter]-shift)<=maxShift)
				validCount++;
		}
		double validShiftPercent = (double)validCount/(double)count;
		
		double score = 0;
		if (overlapSize>0 && validShiftPercent>SHIFT_CONSENSUS_PERCENTAGE)
			score = (double)count/(double)(overlapSize);

		return new Pair<Double, Integer>(score, shift);
	}

	public Pair<Double, Integer> getFullScore(SequenceMinHashes s, int maxShift)
	{
		return getFullScore(getFullHashes(), s, maxShift);
	}
	
	public MinHash getMainHashes()
	{
		return this.mainHashes;
	}
	
	public SequenceId getSequenceId()
	{
		return this.id;
	}

	public int getSequenceLength()
	{
		return this.mainHashes.getSequenceLength();
	}
}
