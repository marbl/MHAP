package com.secret.fastalign.simhash;

import org.apache.lucene.util.OpenBitSet;

import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.data.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public class AbstractSequenceBitHash implements VectorHash<AbstractSequenceBitHash>, Comparable<AbstractSequenceBitHash>
{
	protected OpenBitSet bits;
	protected final SequenceId id;
	protected final int length;

	public AbstractSequenceBitHash(Sequence seq, int kmerSize)
	{
		this.id = seq.getId();
		this.length = seq.length();
	}

	public double adjScore(SequenceSimHash sh)
	{
		double score = correlation(sh);
		
		return score;
	}

	@Override
	public int compareTo(AbstractSequenceBitHash sim)
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

	public boolean getBit(int index)
	{
		return this.bits.fastGet(index);
	}

	public int getIntersectionCount(AbstractSequenceBitHash sh)
	{
		if (this.bits.size()!=sh.bits.size())
			throw new FastAlignRuntimeException("Size of bits in tables must match.");
		
		return (int)(this.bits.size()-OpenBitSet.xorCount(this.bits, sh.bits));
	}

	@Override
	public SequenceId getSequenceId()
	{
		return this.id;
	}

	@Override
	public double correlation(AbstractSequenceBitHash sh)
	{
		int count = getIntersectionCount(sh);
		
		return (double)count/(double)this.bits.size();
	}

	public int seqLength()
	{
		return this.length;
	}

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
}