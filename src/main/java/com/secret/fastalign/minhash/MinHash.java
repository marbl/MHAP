package com.secret.fastalign.minhash;

import java.util.Arrays;

import com.secret.fastalign.general.AbstractSequenceHashes;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Utils;

public final class MinHash extends AbstractSequenceHashes<MinHash>
{
	public final int[] minHashes;
	
	public MinHash(Sequence seq, int kmerSize, int numWords)
	{
		super(seq);
		
		//get the hashes
		this. minHashes = Utils.computeKmerMinHashes(seq, kmerSize, numWords);
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

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "MinHash "+Arrays.toString(this.minHashes) + "";
	}
}
