package com.secret.fastalign.hash;

import com.secret.fastalign.utils.Utils;

public final class Sequence
{
	private final String sequence;
	private final boolean[] ori;
	private final String[] kmers;
	private final SequenceId id;
	//private final int kmerSize;
	
	public Sequence(String sequence, int kmerSize, SequenceId id)
	{
		this.sequence = sequence;
		this.id = id;
		
		int numKmers = this.sequence.length()-kmerSize+1;
		
		//this.kmerSize = kmerSize;
		this.kmers = new String[numKmers];
		this.ori = new boolean[numKmers];
		
    for (int i = 0; i <= this.sequence.length() - kmerSize; i++) {
      String fmer = this.sequence.substring(i, i+kmerSize);
      String rmer = Utils.rc(fmer);
      String currMer = null;
      boolean isFwd = true;
      if (fmer.compareTo(rmer) <= 0) {
         currMer = fmer;
         isFwd = true;
      } else {
         currMer = rmer;
         isFwd = false;
      }

      // store the canonical kmer and the orientation (i.e. whether cannonical was fwd or rev)
      this.kmers[i] = currMer;
      this.ori[i] = isFwd;
   }
	}
	
	public String getString()
	{
		return this.sequence;
	}
	
	public SequenceId getId()
	{
		return this.id;
	}
	
	public String getKmer(int index)
	{
		return this.kmers[index];
	}
	
	public int numKmers()
	{
		return this.kmers.length;
	}
	
	public boolean getOrientation(int index)
	{
		return this.ori[index];
	}
	
	public int length()
	{
		return this.sequence.length();
	}
}
