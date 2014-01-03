package com.secret.fastalign.hash;

public final class SequenceId
{
	private final Long id;
	
	public SequenceId(long val)
	{
		this.id = val;
	}
	
	public long getLongId()
	{
		return this.id;
	}
}
