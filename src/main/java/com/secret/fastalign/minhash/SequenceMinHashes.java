package com.secret.fastalign.minhash;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import com.secret.fastalign.general.AbstractReducedSequence;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Pair;
import com.secret.fastalign.utils.SortablePair;
import com.secret.fastalign.utils.Utils;

public final class SequenceMinHashes extends AbstractReducedSequence<MinHash,SequenceMinHashes>
{
	private final ArrayList<SortablePair<Integer, Integer>> completeHash;
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numWords, int subKmerSize)
	{
		super(seq.getId(), new MinHash(seq, kmerSize, numWords));
		
		//compute just direct hash of sequence
		int[][] hashes = Utils.computeKmerHashesInt(seq, subKmerSize, 1);		
		this.completeHash = new ArrayList<SortablePair<Integer, Integer>>(hashes.length);
		
		for (int iter=0; iter<hashes.length; iter++)
			this.completeHash.add(new SortablePair<Integer, Integer>((int)hashes[iter][0],iter));
		
		//sort the results
		Collections.sort(this.completeHash);
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
			ArrayList<Integer> posShift = new ArrayList<Integer>(this.completeHash.size());
			
			count = 0;
			Iterator<SortablePair<Integer, Integer>> iter1 = this.completeHash.iterator();
			Iterator<SortablePair<Integer, Integer>> iter2 = s.completeHash.iterator();
			
			SortablePair<Integer, Integer> s1 = iter1.next();
			SortablePair<Integer, Integer> s2 = iter2.next();			
			while (iter1.hasNext() && iter2.hasNext())
			{
				if (s1.x < s2.x || s1.y<valid1Lower || s1.y>valid1Upper)
					s1 = iter1.next();
				else
				if (s2.x < s1.x || s2.y<valid2Lower || s2.y>valid2Upper)
					s2 = iter2.next();
				else
				{
					count++;
					
					posShift.add(s2.y-s1.y);
					
					s1 = iter1.next();
					s2 = iter2.next();
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
