package com.secret.fastalign.simhash;

import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.data.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public class AbstractSequenceBitHash implements VectorHash<AbstractSequenceBitHash>, Comparable<AbstractSequenceBitHash>
{
	protected long bits[];
	protected final SequenceId id;
	protected final int length;

	public AbstractSequenceBitHash(Sequence seq, int kmerSize)
	{
		this.id = seq.getId();
		this.length = seq.length();
	}

	public final double adjScore(SequenceSimHash sh)
	{
		double score = correlation(sh);
		
		return score;
	}

	@Override
	public int compareTo(final AbstractSequenceBitHash sim)
	{
		for (int bitIndex=0; bitIndex<this.bits.length; bitIndex++)
		{
			if (this.bits[bitIndex] < this.bits[bitIndex])
				return -1;
			if (this.bits[bitIndex] > this.bits[bitIndex])
				return 1;
		}
		
		return 0;
	}

	public final int getIntersectionCount(final AbstractSequenceBitHash sh)
	{
		if (this.bits.length!=sh.bits.length)
			throw new FastAlignRuntimeException("Size of bits in tables must match.");
		
		int count = 0;
	  for (int longIndex=0; longIndex<this.bits.length; longIndex++)
	  {	      
	  	final long xor = this.bits[longIndex]^sh.bits[longIndex];
	  	
      count += Long.bitCount(xor);
	  }
		
		return this.bits.length*64-count;
	}

	@Override
	public final SequenceId getSequenceId()
	{
		return this.id;
	}

	@Override
	public final double correlation(final AbstractSequenceBitHash sh)
	{
		int count = getIntersectionCount(sh);
		
		return ((double)count/(double)(this.bits.length*64)-0.5)*2.0;
	}

	public final int seqLength()
	{
		return this.length;
	}

	@Override
	public String toString()
	{
		StringBuilder s = new StringBuilder();
	  for (int longIndex=0; longIndex<this.bits.length; longIndex++)
	  {
	  	long mask = 1L << 63;
	  	
      for (int bit=63; bit<=0; bit++)
      {
        if ((this.bits[longIndex]&mask)==0)
  				s.append("0");
        else
  				s.append("1");
        
        mask = mask >>> 1;
      }
	  }		
		return s.toString();
	}
}