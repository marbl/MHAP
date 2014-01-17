package com.secret.fastalign.simhash;

import org.apache.lucene.util.OpenBitSet;

import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;
import com.secret.fastalign.data.Sequence;

public final class SequenceSimHash extends AbstractSequenceBitHash 
{
	private final static int NUM_LONG_BITS = 64;
	
	public static long[][] computeHash(Sequence seq, int kmerSize, int numWords, int subKmerSize)
	{
		String seqString = seq.getString();
		int numberKmers = seqString.length()-kmerSize+1;

		int numerSubKmers = kmerSize-subKmerSize+1;
		long[][] hashes = new long[numberKmers][numWords*numerSubKmers];
		//long[][] hashes = new long[numberKmers][numWords];
		

		HashFunction[] hf = new HashFunction[numWords];
		for (int hashIndex=0; hashIndex<numWords; hashIndex++)
		  hf[hashIndex] = Hashing.murmur3_128(hashIndex);

		//System.out.println(seqString.length());
		for (int kmerIndex=0; kmerIndex<numberKmers; kmerIndex++)
		{
			String subSeq = seqString.substring(kmerIndex, kmerIndex+kmerSize);
			
			for (int hashIndex=0; hashIndex<numWords; hashIndex++)
			{
				for (int i=0; i<numerSubKmers; i++)
				{
					String is = subSeq.substring(i, i+subKmerSize);
					long val = hf[hashIndex].newHasher().putUnencodedChars(is).putInt(1).hash().asLong();
					
					hashes[kmerIndex][hashIndex*numerSubKmers+i] = val;
				}
			}
		}
		
		return hashes;
	}
	
	public SequenceSimHash(Sequence seq, int kmerSize, int numberWords)
	{
		super(seq, kmerSize);
		
		int subKmerSize = kmerSize;
		
		//compute the hashes
		long[][] hashes = computeHash(seq, kmerSize, numberWords, subKmerSize);
		
		recordHashes(hashes, kmerSize, numberWords, subKmerSize);
	}
	
	private void recordHashes(long[][] hashes, int kmerSize, int numWords, int subKmerSize)
	{
		this.bits = new OpenBitSet(NUM_LONG_BITS*numWords);
		int[] counts = new int[this.bits.length()];

		int numerSubKmers = kmerSize-subKmerSize+1;

		//perform count for each kmer
		for (int kmerIndex=0; kmerIndex<hashes.length; kmerIndex++)
		{
			OpenBitSet val = new OpenBitSet(hashes[kmerIndex], hashes[kmerIndex].length);
			
			//for each hash go through all its bits and count
			for (int wordsBitIndex=0; wordsBitIndex<counts.length; wordsBitIndex++)
			{
				for (int subIndex=0; subIndex<numerSubKmers; subIndex++)
				{
					if (val.fastGet(wordsBitIndex*numerSubKmers+subIndex))
						counts[wordsBitIndex]++;
					else
						counts[wordsBitIndex]--;
				}
			}
		}
		
		for (int bitIndex=0; bitIndex<counts.length; bitIndex++)
		{
			if (counts[bitIndex]>0)
				this.bits.fastSet(bitIndex);
		}		
	}
}
