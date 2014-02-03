package com.secret.fastalign.simhash;

import com.secret.fastalign.general.AbstractSequenceHashes;
import com.secret.fastalign.general.Sequence;

public final class SimHash extends AbstractSequenceBitHash<SimHash>
{
	public SimHash(Sequence seq, int kmerSize, int numberWords)
	{
		super(seq);
		
		//compute the hashes
		long[][] hashes = AbstractSequenceHashes.computeKmerHashes(seq, kmerSize, numberWords);
		
		recordHashes(hashes, kmerSize, numberWords);
	}
	
	private final void recordHashes(final long[][] hashes, final int kmerSize, final int numWords)
	{
		final int[] counts = new int[numWords*64];

		//perform count for each kmer
		for (int kmerIndex=0; kmerIndex<hashes.length; kmerIndex++)
		{
			long[] kmerHashes = hashes[kmerIndex];
			
		  for (int longIndex=0; longIndex<numWords; longIndex++)
		  {	      
		  	final long val = kmerHashes[longIndex];
		  	final int offset = longIndex*64;

		  	long mask = 0b1;
		  	
	      for (int bit=0; bit<64; bit++)
	      {
	        /* if not different then increase count */
	        if ((val&mask)==0b0)
	          counts[offset+bit]--;
	        else
	        	counts[offset+bit]++;
	        	
	        mask = mask << 1;
	      }
		  }		  
		}
		
		this.bits = new long[numWords];
	  for (int longIndex=0; longIndex<numWords; longIndex++)
	  {	      
	  	final int offset = longIndex*64;
	  	long val = 0b0;
	  	long mask = 0b1;
	  	
      for (int bit=0; bit<64; bit++)
      {
      	if (counts[offset+bit]>0)
      		val = val | mask;
      	
      	//adjust the mask
        mask = mask << 1;
      }
      
      this.bits[longIndex] = val;
	  }	  
	}
}
