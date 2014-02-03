package com.secret.fastalign.general;

import java.nio.ByteBuffer;
import java.nio.LongBuffer;

import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;
import com.secret.fastalign.utils.FastAlignRuntimeException;


public abstract class AbstractSequenceHashes<T extends AbstractSequenceHashes<T>>
{
	final private SequenceId id;
	final private int length;

	public AbstractSequenceHashes(Sequence seq)
	{
		this.id = seq.getId();
		this.length = seq.length();
	}
	
	public abstract double jaccord(T seqHashes);
	
	public SequenceId getSequenceId()
	{
		return this.id;
	}

	public int sequenceLength()
	{
		return this.length;
	}

	public final static long[][] computeKmerHashes(final Sequence seq, final int kmerSize, final int numWords)
	{
		if (numWords%2!=0)
			throw new FastAlignRuntimeException("Number of words must be a multiple of 2.");
	
		String seqString = seq.getString();
		final int numberKmers = seqString.length()-kmerSize+1;
		
		if (numberKmers<1)
			throw new FastAlignRuntimeException("Kmer size bigger than string length.");
	
		final long[][] hashes = new long[numberKmers][numWords];
	
		//might want to change to city hash if it comes out
		HashFunction hf = Hashing.murmur3_128(0);
		HashFunction hfInt = Hashing.murmur3_32(0);
		
		//RabinKarpSeqHash rabinHash = new RabinKarpSeqHash(kmerSize);
		//final int[] rabinHashes = rabinHash.hashInt(seqString);
		
		final int[] rabinHashes = new int[seq.numKmers(kmerSize)];
		for (int iter=0; iter<seq.numKmers(kmerSize); iter++)
		{
			String kmer = seq.getKmer(iter, kmerSize);
			
			rabinHashes[iter] = hfInt.newHasher(0).putUnencodedChars(kmer).hash().asInt();
		}
		
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

}