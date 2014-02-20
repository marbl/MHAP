package com.secret.fastalign.minhash;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import com.secret.fastalign.general.AbstractReducedSequence;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Pair;
import com.secret.fastalign.utils.Utils;

public final class SequenceMinHashes extends AbstractReducedSequence<MinHash,SequenceMinHashes>
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
			int result = Integer.compare(this.x, p.x);
			
			if (result!=0)
				return result;
			
			return Integer.compare(this.y, p.y);
		}		
	}
	
	public final static int MAX_SHIFT_ALLOWED = 200;
	public final static double SHIFT_CONSENSUS_PERCENTAGE = 0.8;
	
	private final int completeHash[][];
	private final int subKmerSize;
	private final Sequence seq;
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numWords, int subSequenceSize, int subKmerSize, boolean storeHashes, HashSet<Integer> filter)
	{
		super(seq.getId(), new MinHash(seq, kmerSize, numWords, subSequenceSize, filter));
		this.subKmerSize = subKmerSize;
		
		if (storeHashes)
		{
			this.completeHash = getFullHashes(seq, subKmerSize);
			this.seq = null;
		}
		else
		{
			this.completeHash = null;
			this.seq = seq;
		}
	}

	private int[][] getFullHashes(Sequence seq, int subKmerSize)
	{
		//compute just direct hash of sequence
		int[][] hashes = Utils.computeKmerHashesInt(seq, subKmerSize, 1);	
		
		SortableIntPair[] completeHashAsPair = new SortableIntPair[hashes.length];	
		for (int iter=0; iter<hashes.length; iter++)
			completeHashAsPair[iter] = new SortableIntPair(hashes[iter][0],iter);
		
		//sort the results
		Arrays.sort(completeHashAsPair);
		
		//store in array to reduce memory
		int[][] completeHash = new int[completeHashAsPair.length][2];
		for (int iter=0; iter<completeHashAsPair.length; iter++)
		{
			completeHash[iter][0] = completeHashAsPair[iter].x;
			completeHash[iter][1] = completeHashAsPair[iter].y;
		}		
		
		return completeHash;
	}
	
	public int[][] getFullHashes()
	{
		if (this.completeHash!=null)
			return this.completeHash;
		
		return getFullHashes(this.seq, this.subKmerSize);
	}

	public Pair<Double, Integer> getFullScore(SequenceMinHashes s)
	{
		return getFullScore(getFullHashes(), s);
	}

	public Pair<Double, Integer> getFullScore(int[][] allKmerHashesw, SequenceMinHashes s)
	{		
		//init the ok regions
		int valid1Lower = 0;
		int valid1Upper = this.getSequenceLength();
		int valid2Lower = 0;
		int valid2Upper = s.getSequenceLength();
		int overlapSize = 0;
		int border = 200;
		
		int[][] allKmerHashes = getFullHashes();
		int[][] sAllKmerHashes = s.getFullHashes();
		
		int count = 0;
		int shift = 0;
		ArrayList<Integer> posShift = null;
		for (int repeat=0; repeat<1; repeat++)
		{
			posShift = new ArrayList<Integer>(allKmerHashes.length);
			
			count = 0;
			int iter1 = 0;
			int iter2 = 0;
			
			while (iter1<allKmerHashes.length && iter2<sAllKmerHashes.length)
			{
				int[] s1 = allKmerHashes[iter1];
				int[] s2 = sAllKmerHashes[iter2];
				
				if (s1[0] < s2[0] || s1[1]<valid1Lower || s1[1]>valid1Upper)
					iter1++;
				else
				if (s2[0] < s1[0] || s1[1]<valid2Lower || s2[1]>valid2Upper)
					iter2++;
				else
				{
					count++;
					
					posShift.add(s2[1]-s1[1]);
					
					iter1++;
					iter2++;
				}
			}
			
			//TODO replace with a selection algorithm maybe
			Collections.sort(posShift);
			
			//get the median
			if (!posShift.isEmpty())
				shift = posShift.get(posShift.size()/2);
			else
				shift = 0;
			
			//System.out.println(posShift);

			valid1Lower = Math.max(0, -shift-border);
			valid1Upper = Math.min(getSequenceLength(), s.getSequenceLength()-shift+border);
			valid2Lower = Math.max(0, shift-border);
			valid2Upper = Math.min(s.getSequenceLength(), getSequenceLength()+shift+border);

			int valid2LowerBorder = Math.max(0, shift);
			int valid2UpperBorder = Math.min(s.getSequenceLength(), getSequenceLength()+shift);
			overlapSize = valid2UpperBorder-valid2LowerBorder;
			
			//System.out.println("Size1= "+getSequenceLength()+" Lower:"+ valid1Lower+" Upper:"+valid1Upper+" Shift="+shift);
			//System.out.println("Size2= "+s.getSequenceLength()+" Lower:"+ valid2Lower+" Upper:"+valid2Upper);
		}
		
		//count percent valid shift, there must be a consensus 
		int percentValid = 0;
		for (int currShift : posShift)
		{
			if (Math.abs(currShift-shift)<=MAX_SHIFT_ALLOWED)
				percentValid++;
		}
		double validShiftPercent = (double)percentValid/(double)posShift.size();
		
		double score = 0;
		if (overlapSize>0 && validShiftPercent>SHIFT_CONSENSUS_PERCENTAGE)
			score = (double)count/(double)(overlapSize);

		return new Pair<Double, Integer>(score, shift);
	}
}
