package com.secret.fastalign.data;

import com.secret.fastalign.utils.HashCodeUtil;

public final class SequenceId
{
	private final int id;
	
	private final boolean isFwd;
	
	public SequenceId(int val, boolean isFwd)
	{
		this.id = val;
		this.isFwd = isFwd;
	}
	
	public SequenceId complimentId()
	{
		return new SequenceId(this.id, !this.isFwd);
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
		
		return (this.id == other.id) && (this.isFwd == other.isFwd);
	}
	
	public boolean isForward()
	{
		return this.isFwd;
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
		int hash = HashCodeUtil.hash(HashCodeUtil.SEED, this.id);
		return HashCodeUtil.hash(hash, this.id);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return ""+this.id+(this.isFwd ? "(fwd)" : "(rev)");
	}
	
	public String toStringInt()
	{
		return ""+this.id+(this.isFwd ? " 1 " : " 0 ");
	}
}
