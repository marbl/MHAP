package com.secret.fastalign.general;

public abstract class AbstractSequenceHashes<T extends AbstractSequenceHashes<T>>
{
	final private int length;

	public AbstractSequenceHashes(Sequence seq)
	{
		this.length = seq.length();
	}
	
	public abstract double jaccord(T seqHashes);
	
	public int sequenceLength()
	{
		return this.length;
	}

}