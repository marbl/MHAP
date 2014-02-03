package com.secret.fastalign.simhash;

import com.secret.fastalign.general.AbstractSequenceHashes;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public abstract class AbstractSequenceBitHash<T extends AbstractSequenceBitHash<T>> extends AbstractSequenceHashes<T> implements VectorHash<T>, Comparable<T>
{
	protected long bits[];
	
	public AbstractSequenceBitHash(Sequence seq)
	{
		super(seq);
	}

	public final double adjScore(T sh)
	{
		double score = jaccord(sh);
		
		return score;
	}
	
	public long[] getBits()
	{
		return this.bits;
	}

	@Override
	public int compareTo(final T sim)
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

	public final int getIntersectionCount(final T sh)
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
	public final double jaccord(final T sh)
	{
		int count = getIntersectionCount(sh);
		
		return ((double)count/(double)(this.bits.length*64)-0.5)*2.0;
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