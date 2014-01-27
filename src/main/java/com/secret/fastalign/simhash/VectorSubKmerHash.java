package com.secret.fastalign.simhash;

import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.data.SequenceId;

public final class VectorSubKmerHash implements VectorHash<VectorSubKmerHash>
{
	private static final int NUM_LETTERS = 4;
	
	private final int[] kmerVector;
	private final double kmerNorm;

	private final SequenceId sequenceId;
	

	public VectorSubKmerHash(Sequence seq, int kmerSize)
	{
		this.sequenceId = seq.getId();
		this.kmerVector = vectorCounts(seq.getString(), kmerSize);
		
		int sum = 0;
		for (int val : this.kmerVector)
			sum+=val*val;
		
		this.kmerNorm = Math.sqrt((double)sum);
	}
	
	private int[] vectorCounts(String seq, int kmerSize)
	{
		int[] counts = new int[NUM_LETTERS*kmerSize];
		
		char[] charSeq = seq.toCharArray();

		int size = seq.length()-kmerSize+1;
		for (int kmerIndex=0; kmerIndex<size;kmerIndex++)
		{
			int[] vals = kmerValue(charSeq, kmerIndex, kmerSize);

			for (int iter=0; iter<counts.length; iter++)
				counts[iter] += vals[iter];
		}
		
		return counts;
	}
	
	private static int[] kmerValue(char[] seq, int start, int kmerSize)
	{
		int subKmerSize = 6;
		
		int numWords = (int)Math.pow(NUM_LETTERS, subKmerSize);
		int[] vals = new int[numWords*kmerSize];
		
		for (int letterIndex=0; letterIndex<kmerSize-subKmerSize+1; letterIndex++)
		{
			int indexPosition = 0;
			for (int subIndex=0; subIndex<subKmerSize; subIndex++)
			{
				char c = seq[start+letterIndex+subIndex];

				//ACGT
				switch (c)
				{
					case 'A': indexPosition = (indexPosition<<2)+0; break;
					case 'C': indexPosition = (indexPosition<<2)+1; break;
					case 'G': indexPosition = (indexPosition<<2)+2; break;
					case 'T': indexPosition = (indexPosition<<2)+3; break;
				}	
			}
			
			//record the letter
			vals[letterIndex*numWords+indexPosition] = 1;
		}
		
		return vals;
	}

	@Override
	public double correlation(VectorSubKmerHash val)
	{
		int sum = 0;
		for (int iter=0; iter<this.kmerVector.length; iter++)
			sum += this.kmerVector[iter]*val.kmerVector[iter];
		
		return ((double)sum)/(this.kmerNorm*val.kmerNorm)*1000;
		//return (double)sum;
	}

	@Override
	public SequenceId getSequenceId()
	{
		return this.sequenceId;
	}

}
