package com.secret.fastalign.simhash;

import java.nio.ByteBuffer;
import java.nio.LongBuffer;
import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.RabinKarpSeqHash;

public final class SequenceSimHash extends AbstractSequenceBitHash 
{
	public final static long[][] computeHash(final Sequence seq, final int kmerSize, final int numWords)
	{
		if (numWords%2!=0)
			throw new FastAlignRuntimeException("Number of words must be a multiple of 2.");

		String seqString = seq.getString();
		final int numberKmers = seqString.length()-kmerSize+1;
		
		final long[][] hashes = new long[numberKmers][numWords];

		//might want to change to city hash if it comes out
		HashFunction hf = Hashing.murmur3_128(0);
		
		RabinKarpSeqHash rabinHash = new RabinKarpSeqHash(kmerSize);

		final int[] rabinHashes = rabinHash.hashInt(seqString);
		
		for (int iter=0; iter<rabinHashes.length; iter++)
		{
			for (int word128=0; word128<numWords/2; word128++)
			{
				final Hasher hasher = hf.newHasher(0);
				final HashCode code = hasher.putInt(rabinHashes[iter]).putInt(word128).hash();

				//store the code
				LongBuffer bb = ByteBuffer.wrap(code.asBytes()).asLongBuffer();
				hashes[iter][word128*2+0] = bb.get(0);
				hashes[iter][word128*2+1] = bb.get(1);
			}
		}
		
		return hashes;
	}
	
	public SequenceSimHash(Sequence seq, int kmerSize, int numberWords)
	{
		super(seq, kmerSize);
		
		//compute the hashes
		long[][] hashes = computeHash(seq, kmerSize, numberWords);
		
		recordHashes(hashes, kmerSize, numberWords);
	}
	
	private final void recordHashes(final long[][] hashes, final int kmerSize, final int numWords)
	{
		final int[] counts = new int[numWords*64];

		//perform count for each kmer
		for (int kmerIndex=0; kmerIndex<hashes.length; kmerIndex++)
		{
		  for (int longIndex=0; longIndex<numWords; longIndex++)
		  {	      
		  	final long val = hashes[kmerIndex][longIndex];
		  	final int offset = longIndex*64;

		  	long mask = 0b1;
		  	
	      for (int bit=0; bit<64; bit++)
	      {
	        /* if not different then increase count */
	        if ((val&mask)==0b0)
	          counts[offset+bit]++;
	        else
	        	counts[offset+bit]--;
	        	
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
