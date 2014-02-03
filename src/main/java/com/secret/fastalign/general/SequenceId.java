package com.secret.fastalign.general;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.utils.HashCodeUtil;

public final class SequenceId
{
	private static ConcurrentHashMap<Long, String> indicies = new ConcurrentHashMap<Long, String>();
	private static AtomicLong globalCounter = new AtomicLong(-1);
	
	private final long id;
	private final boolean isFwd;
	private final int hash;
	
	public SequenceId(String id, boolean isFwd)
	{
		this.id = globalCounter.addAndGet(1);
		//indicies.put(this.id, id);
		
		this.isFwd = isFwd;
		this.hash = myHashCode();
	}
	
	private SequenceId(long id, boolean isFwd)
	{
		this.id = id;
		this.isFwd = isFwd;
		this.hash = myHashCode();
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
		
		if (other.hash!=this.hash)
			return false;
		
		return (this.id==other.id) && (this.isFwd == other.isFwd);
	}
	
	public boolean isForward()
	{
		return this.isFwd;
	}
	
	public long getHeaderId()
	{
		return this.id;
	}

	public String getHeader()
	{
		String s = indicies.get(this.id);
		if (s!=null)
			return s;
		
		return String.valueOf(this.id);
	}

	private int myHashCode()
	{
		int hash = HashCodeUtil.hash(HashCodeUtil.SEED, this.id);
		return HashCodeUtil.hash(hash, this.isFwd);		
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode()
	{
		return this.hash;
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
