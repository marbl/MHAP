package com.secret.fastalign.utils;

public final class RollingSequenceHash
{
	private final int kmerSize;

	public RollingSequenceHash(int kmerSize)
	{
		this.kmerSize = kmerSize;
	}
	
	public final int[] hashInt(final String seq)
	{
		//get the string
		final char[] seqArray = seq.toCharArray();
		
		//RabinKarpHash ch = new RabinKarpHash(this.kmerSize);
		CyclicHash ch = new CyclicHash(this.kmerSize);
		
		//allocate hashes
		int[] hashes = new int[seqArray.length-this.kmerSize+1];
		
		int kmerIndex;
		for(kmerIndex = 0; kmerIndex<this.kmerSize-1; kmerIndex++) 
		{
			ch.eat(seqArray[kmerIndex]);
		}
		
		int count = 0;
		int hash = ch.eat(seqArray[kmerIndex++]);		
		hashes[count++] = hash;
		
		//start rolling
		for(; kmerIndex<seqArray.length; kmerIndex++) 
		{
			//System.out.println(k-this.kmerSize);
			hash = ch.update(seqArray[kmerIndex-this.kmerSize], seqArray[kmerIndex]);			
			hashes[count++] = hash;
		}
		
		return hashes;
	}
}
