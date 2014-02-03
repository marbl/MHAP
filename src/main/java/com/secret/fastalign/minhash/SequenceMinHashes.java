package com.secret.fastalign.minhash;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;

import com.secret.fastalign.general.AbstractOrderedSequenceHashes;
import com.secret.fastalign.general.AbstractSequenceHashes;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Pair;
import com.secret.fastalign.utils.SortablePair;

public final class SequenceMinHashes extends AbstractOrderedSequenceHashes<MinHash,SequenceMinHashes>
{
	private final ArrayList<SortablePair<Integer, Integer>> completeHash;
	
	public SequenceMinHashes(Sequence seq, int kmerSize, int numWords, int subKmerSize)
	{
		//super(new MinHash(seq, kmerSize, numWords), generateSubHashes(seq, subStringSize, subKmerSize, subWordSize));
		super(new MinHash(seq, kmerSize, numWords), null);
		
		//compute just direct hash of sequence
		long[][] hashes = AbstractSequenceHashes.computeKmerHashes(seq, subKmerSize, 2);		
		this.completeHash = new ArrayList<SortablePair<Integer, Integer>>(hashes.length);
		
		for (int iter=0; iter<hashes.length; iter++)
			this.completeHash.add(new SortablePair<Integer, Integer>((int)hashes[iter][0],iter));
		
		//sort the results
		Collections.sort(this.completeHash);
	}
	
	/*
	private static ArrayList<MinHash> generateSubHashes(Sequence seq, int subStringSize, int subKmerSize, int subWordSize)
	{
		//generate the array of simhashes
		ArrayList<MinHash> subHashes = new ArrayList<MinHash>(seq.length()-subKmerSize+1);
		for (int iter=0; iter<seq.length()-subKmerSize+1; iter+=subStringSize)
		{
			String subString = seq.getString().substring(iter, Math.min(iter+subStringSize, seq.length()));
			Sequence subSequence = new Sequence(subString, seq.getId());
			
			subHashes.add(new MinHash(subSequence, subKmerSize, subWordSize));
		}
		
		return subHashes;
	}
	*/

	/* (non-Javadoc)
	 * @see com.secret.fastalign.general.AbstractOrderedSequenceHashes#orderedScore(com.secret.fastalign.general.AbstractOrderedSequenceHashes)
	 */
	@Override
	public Pair<Double, Integer> orderedScore(SequenceMinHashes s)
	{		
		//init the ok regions
		int valid1Lower = 0;
		int valid1Upper = this.getSequenceLength();
		int valid2Lower = 0;
		int valid2Upper = s.getSequenceLength();
		
		int count = 0;
		int shift = 0;
		int border = 0;
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
			
			//System.out.println("Size1= "+getSequenceLength()+" Lower:"+ valid1Lower+" Upper:"+valid1Upper+" Shift="+shift);
			//System.out.println("Size2= "+s.getSequenceLength()+" Lower:"+ valid2Lower+" Upper:"+valid2Upper);
		}
		
		int overlapSize = valid2Upper-valid2Lower-border*2;
		//overlapSize = Math.min(getSequenceLength(), s.getSequenceLength());
		
		System.out.println(overlapSize);
		
		double score = 0;
		if (overlapSize>0)
			score = (double)count/(double)(overlapSize);
		
		return new Pair<Double, Integer>(score, shift);
	}
}
