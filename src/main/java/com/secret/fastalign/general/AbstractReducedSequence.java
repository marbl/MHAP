package com.secret.fastalign.general;

public class AbstractReducedSequence<H extends AbstractSequenceHashes<H>, T extends AbstractReducedSequence<H,T>>
{
	private final SequenceId id;
	private final H mainHashes;
	
	public AbstractReducedSequence(SequenceId id, H mainHashes)
	{
		this.mainHashes = mainHashes;
		this.id = id;
	}
	
	public SequenceId getSequenceId()
	{
		return this.id;
	}
	
	public int getSequenceLength()
	{
		return this.mainHashes.sequenceLength();
	}
	
	public H getMainHashes()
	{
		return this.mainHashes;
	}
	
	public double jaccord(T s)
	{
		return getMainHashes().jaccord(s.getMainHashes());
	}
}