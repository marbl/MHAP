package com.secret.fastalign.simhash;

import org.apache.lucene.util.OpenBitSet;

import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.data.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public final class SequenceSimHash implements Comparable<SequenceSimHash>
{
	private final OpenBitSet bits;
	private final SequenceId id;
	private final int length;
	
	private final static int NUM_LONG_BITS = 64;
	
	public SequenceSimHash(Sequence seq, int kmerSize, int numberWords)
	{
		this.bits = new OpenBitSet(NUM_LONG_BITS*numberWords);
		this.length = seq.length();
		this.id = seq.getId();
		
		int[] counts = new int[numberWords*NUM_LONG_BITS];
		
		//compute the hashes
		long[][] hashes = computeHash(seq, kmerSize);
		
		//perform count
		for (int kmerIndex=0; kmerIndex<hashes.length; kmerIndex++)
		{
			OpenBitSet val = new OpenBitSet(hashes[kmerIndex], numberWords);
			
			//for each hash go through all its bits and count
			for (int bitIndex=0; bitIndex<val.length(); bitIndex++)
			{
				if (val.fastGet(bitIndex))
					counts[bitIndex]++;
				else
					counts[bitIndex]--;
			}
		}
		
		for (int bitIndex=0; bitIndex<counts.length; bitIndex++)
		{
			if (counts[bitIndex]>0)
				this.bits.fastSet(bitIndex);
		}
	}
	
	public double adjScore(SequenceSimHash sh)
	{
		double score = score(sh);
		
		return score;
	}
	
	private long[][] computeHash(Sequence seq, int kmerSize)
	{
		int numberLongs = (int)this.bits.size()/NUM_LONG_BITS;
		
		String seqString = seq.getString();
		int size = seqString.length()-kmerSize;

		long[][] hashes = new long[size][numberLongs];
		

		HashFunction[] hf = new HashFunction[numberLongs];
		for (int hashIndex=0; hashIndex<numberLongs; hashIndex++)
		  hf[hashIndex] = Hashing.murmur3_128(hashIndex);

		//System.out.println(seqString.length());
		for (int kmerIndex=0; kmerIndex<size; kmerIndex++)
		{
			CharSequence subSeq = seqString.subSequence(kmerIndex, kmerIndex+kmerSize);
			for (int hashIndex=0; hashIndex<numberLongs; hashIndex++)
			{
				hashes[kmerIndex][hashIndex] = hf[hashIndex].newHasher().putUnencodedChars(subSeq).hash().asLong();
			}
		}
		
		return hashes;
	}
	
	public boolean getBit(int index)
	{
		return this.bits.fastGet(index);
	}
	
	public int getIntersectionCount(SequenceSimHash sh)
	{
		if (this.bits.size()!=sh.bits.size())
			throw new FastAlignRuntimeException("Size of bits in tables must match.");
		
		return (int)(this.bits.size()-OpenBitSet.xorCount(this.bits, sh.bits));
	}
	
	public SequenceId getSequenceId()
	{
		return this.id;
	}
	
	public double score(SequenceSimHash sh)
	{
		int count = getIntersectionCount(sh);
		
		return (double)count/(double)this.bits.size();
	}
	
	public int seqLength()
	{
		return this.length;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		StringBuilder s = new StringBuilder();
		for (int index=0; index<this.bits.length(); index++)
			if (this.bits.fastGet(index))
				s.append("1");
			else
				s.append("0");
		
		return s.toString();
	}

	@Override
	public int compareTo(SequenceSimHash sim)
	{
		for (int bitIndex=0; bitIndex<this.bits.length(); bitIndex++)
		{
			boolean v1 = this.bits.fastGet(bitIndex);
			boolean v2 = sim.bits.fastGet(bitIndex);	
			
			if (!v1 && v2)
				return -1;
			if (v1 && !v2)
				return 1;
		}
		
		return 0;
	}
}
