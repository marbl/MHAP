package com.secret.fastalign.utils;

import java.util.Random;

public class RabinKarpSeqHash
{
	private final int kmerSize;

	public RabinKarpSeqHash(int kmerSize)
	{
		this.kmerSize = kmerSize;
	}
	
	public int[] hashInt(String seq)
	{
		return hashInt(seq, 0);
	}

	public int[] hashInt(String seq, int seed)
	{
		Random rn = new Random(seed);
		seed = rn.nextInt();
		
		//get the string
		char[] seqArray = seq.toCharArray();
		
		RabinKarpHash ch = new RabinKarpHash(this.kmerSize);
		
		//allocate hashes
		int[] hashes = new int[seqArray.length-this.kmerSize+1];
		
		int kmerIndex;
		for(kmerIndex = 0; kmerIndex<this.kmerSize-1; kmerIndex++) 
		{
			ch.eat(seqArray[kmerIndex]);
		}
		
		int hash = ch.eat(seqArray[kmerIndex++]);
		
		int count = 0;
		
		hashes[count++] = hash*22695477+seed;
		
		//start rolling
		for(; kmerIndex<seqArray.length; kmerIndex++) 
		{
			//System.out.println(k-this.kmerSize);
			hash = ch.update(seqArray[kmerIndex-this.kmerSize], seqArray[kmerIndex]);
			
			hashes[count++] = hash*22695477+seed;
		}
		
		return hashes;
	}
}
