package com.secret.fastalign.data;

import com.secret.fastalign.utils.HashCodeUtil;

public final class SequenceId
{
	private final int id;
	
	public SequenceId(int val)
	{
		this.id = val;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SequenceId other = (SequenceId) obj;
		
		return this.id == other.id;
	}

	public long getLongId()
	{
		return this.id;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode()
	{
		return HashCodeUtil.hash(HashCodeUtil.SEED, this.id);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return ""+this.id;
	}
}
