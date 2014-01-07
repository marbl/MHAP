package com.secret.fastalign.minhash;

import com.secret.fastalign.data.SequenceId;

public final class KmerInfo
{
	private final SequenceId id;
	private final int pos;
	
	public KmerInfo(SequenceId id, int pos)
	{
		this.id = id;
		this.pos = pos;
	}
	
	public SequenceId getId()
	{
		return this.id;		
	}
	
	public int getPosition()
	{
		return this.pos;		
	}

}
