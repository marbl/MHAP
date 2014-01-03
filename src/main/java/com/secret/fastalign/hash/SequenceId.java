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

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return this.id.toString();
	}
}
