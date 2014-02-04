package com.secret.fastalign.minhash;

import java.util.Arrays;

import com.secret.fastalign.general.AbstractSequenceHashes;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Utils;

public final class MinHash extends AbstractSequenceHashes<MinHash>
{
	public final int[] minHashes;
	public final int[] hashPositions;
	
	public MinHash(Sequence seq, int kmerSize, int numWords)
	{
		super(seq);
		
		long[][] hashes = Utils.computeKmerHashes(seq, kmerSize, numWords/2);
				
		int[] minHashes = new int[numWords];
		int[] pos = new int[numWords];
		
		Arrays.fill(minHashes, Integer.MAX_VALUE);
		for (int kmerIndex=0; kmerIndex<hashes.length; kmerIndex++)
			for (int hashIndex=0; hashIndex < numWords; hashIndex+=2)
			{
				//convert to int
				long val = hashes[kmerIndex][hashIndex/2];
				int i1 = (int)(val >> 32);
				int i2 = (int)val;
				
				if (i1<minHashes[hashIndex])
					minHashes[hashIndex] = i1;
				if (i2<minHashes[hashIndex+1])
					minHashes[hashIndex+1] = i2;
			}
		
		this.minHashes = minHashes;
		this.hashPositions = pos;
	}

	@Override
	public double jaccord(MinHash seqHashes)
	{
		int count = 0;
		for (int iter=0; iter<numHashes(); iter++)
			if (this.minHashes[iter]==seqHashes.minHashes[iter])
				count++;
		
		return (double)count/(double)numHashes();
	}

	/**
	 * @return the minHashes
	 */
	public final int[] getMinHashes()
	{
		return this.minHashes;
	}
	
	public final int numHashes()
	{
		return this.minHashes.length;
	}

	/**
	 * @return the hashPositions
	 */
	public final int[] getHashPositions()
	{
		return this.hashPositions;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "MinHash "+Arrays.toString(this.minHashes) + "";
	}
}
