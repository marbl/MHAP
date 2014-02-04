package com.secret.fastalign.general;

import com.secret.fastalign.utils.Utils;

public final class Sequence
{
	private final String sequence;
	private final SequenceId id;
	
	public Sequence(int[] sequence, SequenceId id)
	{
		this.id = id;
		
		StringBuilder s = new StringBuilder();
		for (int iter=0; iter<sequence.length; iter++)
		{
			switch(sequence[iter])
			{
				case 0 : s.append("U"); break;
				case 1 : s.append("C"); break;
				case 2 : s.append("G"); break;
				case 3 : s.append("T"); break;
				default : throw new RuntimeException("Uknown integer value.");
			}
		}
		
		this.sequence = s.toString();
	}
	
	public Sequence(String sequence, SequenceId id)
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
	
	public Sequence getReverseCompliment()
	{
		return new Sequence(Utils.rc(this.sequence), this.id.complimentId());
	}
	
	public String getKmer(int index, int kmerSize)
	{
		return this.sequence.substring(index, index+kmerSize);
	}
	
	public int numKmers(int kmerSize)
	{
		return this.sequence.length()-kmerSize+1;
	}
	
	public int length()
	{
		return this.sequence.length();
	}
}
