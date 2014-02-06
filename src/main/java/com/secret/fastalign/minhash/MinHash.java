package com.secret.fastalign.minhash;

import java.util.Arrays;

import com.secret.fastalign.general.AbstractSequenceHashes;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Utils;

public final class MinHash extends AbstractSequenceHashes<MinHash>
{
	private final int[][] minHashes;
	
	public MinHash(Sequence seq, int kmerSize, int numWords)
	{
		super(seq);
		
		int numberSubSeq = seq.length()/MinHashSearch.SUB_SEQUENCE_SIZE+1;
		if (seq.length()%MinHashSearch.SUB_SEQUENCE_SIZE<kmerSize)
			numberSubSeq--;
		
		this.minHashes = new int[numberSubSeq][];
		for (int iter=0; iter<numberSubSeq; iter++)
		{
			String subString = seq.getString().substring(iter*MinHashSearch.SUB_SEQUENCE_SIZE, Math.min(seq.length(), (iter+1)*MinHashSearch.SUB_SEQUENCE_SIZE));

			//get the hashes
			this.minHashes[iter] = Utils.computeKmerMinHashes(subString, kmerSize, numWords);
		}
		
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
	public final int[][] getSubSeqMinHashes()
	{
		return this.minHashes;
	}
	
	public final int numHashes()
	{
		return this.minHashes[0].length;
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
