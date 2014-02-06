package com.secret.fastalign.minhash;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import com.secret.fastalign.general.AbstractReducedSequence;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Pair;
import com.secret.fastalign.utils.SortablePair;
import com.secret.fastalign.utils.Utils;

public final class SequenceMinHashes extends AbstractReducedSequence<MinHash,SequenceMinHashes>
{
	private final class SortableIntPair extends SortablePair<Integer, Integer>
	{
		/**
		 * 
		 */
		private static final long serialVersionUID = -8427321872814544309L;

		public SortableIntPair(Integer x, Integer y)
		{
			super(x, y);
		}		
	}
	
	private final int completeHash[][];
	private final int subKmerSize;
	private final Sequence seq;
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numWords, int subKmerSize, boolean storeHashes)
	{
		super(seq.getId(), new MinHash(seq, kmerSize, numWords));
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
		int border = 0;
		
		int[][] allKmerHashes = getFullHashes();
		int[][] sAllKmerHashes = s.getFullHashes();
		
		int count = 0;
		int shift = 0;
		for (int repeat=0; repeat<1; repeat++)
		{
			ArrayList<Integer> posShift = new ArrayList<Integer>(allKmerHashes.length);
			
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
			overlapSize = valid2Upper-valid2Lower;
			
			//System.out.println("Size1= "+getSequenceLength()+" Lower:"+ valid1Lower+" Upper:"+valid1Upper+" Shift="+shift);
			//System.out.println("Size2= "+s.getSequenceLength()+" Lower:"+ valid2Lower+" Upper:"+valid2Upper);
		}
		
		double score = 0;
		if (overlapSize>0)
			score = (double)count/(double)(overlapSize);
		
		return new Pair<Double, Integer>(score, shift);
	}
}
