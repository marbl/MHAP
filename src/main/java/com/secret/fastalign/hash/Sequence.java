package com.secret.fastalign.hash;

public final class Sequence
{
	private final String sequence;
	private final SequenceId id;
	
	public Sequence(String sequence, int kmerSize, SequenceId id)
	{
		this.sequence = sequence;
		this.id = id;
	}
	
	public String getString()
	{
		return this.sequence;
	}
	
	public SequenceId getId()
	{
		return this.id;
	}
	
	public String getKmer(int index, int kmerSize)
	{
		return this.sequence.substring(index, index+kmerSize);
	}
	
	public int numKmers(int kmerSize)
	{
		return this.sequence.length() - kmerSize;
	}
	
	public int length()
	{
		return this.sequence.length();
	}
}
