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
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numWords, int subKmerSize)
	{
		super(seq.getId(), new MinHash(seq, kmerSize, numWords));
		
		//compute just direct hash of sequence
		int[][] hashes = Utils.computeKmerHashesInt(seq, subKmerSize, 1);	
		
		SortableIntPair[] completeHashAsPair = new SortableIntPair[hashes.length];	
		for (int iter=0; iter<hashes.length; iter++)
			completeHashAsPair[iter] = new SortableIntPair((int)hashes[iter][0],iter);
		
		//sort the results
		Arrays.sort(completeHashAsPair);
		
		//store in array to reduce memory
		this.completeHash = new int[completeHashAsPair.length][2];
		for (int iter=0; iter<completeHashAsPair.length; iter++)
		{
			this.completeHash[iter][0] = completeHashAsPair[iter].x;
			this.completeHash[iter][1] = completeHashAsPair[iter].y;
		}		
	}

	public Pair<Double, Integer> getFullScore(SequenceMinHashes s)
	{		
		//init the ok regions
		int valid1Lower = 0;
		int valid1Upper = this.getSequenceLength();
		int valid2Lower = 0;
		int valid2Upper = s.getSequenceLength();
		int overlapSize = 0;
		int border = 0;
		
		int count = 0;
		int shift = 0;
		for (int repeat=0; repeat<1; repeat++)
		{
			ArrayList<Integer> posShift = new ArrayList<Integer>(this.completeHash.length);
			
			count = 0;
			int iter1 = 0;
			int iter2 = 0;
			
			while (iter1<this.completeHash.length && iter2<s.completeHash.length)
			{
				int[] s1 = this.completeHash[iter1];
				int[] s2 = s.completeHash[iter2];
				
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
