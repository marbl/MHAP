package com.secret.fastalign.general;

public final class SubSequenceId
{
	private final SequenceId id;
	private final short sequenceIndex;

	public SubSequenceId(SequenceId fromId, short fromSubSeq)
	{
		this.id = fromId;
		this.sequenceIndex = fromSubSeq;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + ((this.id == null) ? 0 : this.id.hashCode());
		result = prime * result + this.sequenceIndex;
		return result;
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
		SubSequenceId other = (SubSequenceId) obj;

		if (this.sequenceIndex != other.sequenceIndex)
			return false;
		
		if (this.id == null)
		{
			if (other.id != null)
				return false;
		}
		else if (!this.id.equals(other.id))
			return false;
		return true;
	}

	public SequenceId getId()
	{
		return this.id;
	}

	public short getSequenceIndex()
	{
		return this.sequenceIndex;
	}
	
	
}