package com.secret.fastalign.general;

import java.io.Serializable;

public final class SequenceId implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2181572437818064822L;
	private final int id;
	private final boolean isFwd;
	
	public SequenceId(int id)
	{
		this(id, true);
	}
	
	public SequenceId(int id, boolean isFwd)
	{
		this.id = id;
		this.isFwd = isFwd;
	}
	
	public SequenceId createOffset(int offset)
	{
		return new SequenceId(this.id+offset, this.isFwd);
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
		
		return (this.id==other.id) && (this.isFwd == other.isFwd);
	}
	
	public boolean isForward()
	{
		return this.isFwd;
	}
	
	public int getHeaderId()
	{
		return this.id;
	}

	public String getHeader()
	{
		//String s = indicies.get(this.id);
		//if (s!=null)
		//	return s;
		
		return String.valueOf(this.id);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode()
	{
		return this.isFwd? this.id : -this.id;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return ""+getHeader()+(this.isFwd ? "(fwd)" : "(rev)");
	}
	
	public String toStringInt()
	{
		return ""+getHeader()+(this.isFwd ? " 1" : " 0");
	}
}
