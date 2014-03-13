package com.secret.fastalign.general;


public final class MatchResult implements Comparable<MatchResult>
{
	private final SequenceId fromId;
	private final SequenceId toId;
	private final int a;
	private final int b;
	private final double score;
	
	public MatchResult(SequenceId fromId, SequenceId toId, double score, int a, int b)
	{
		this.fromId = fromId;
		this.toId = toId;
		
		this.a = a;
		this.b = b;
		
		if (score>1.0)
			this.score = 	1.0;
		else
			this.score = score;
	}

	/**
	 * @return the fromId
	 */
	public SequenceId getFromId()
	{
		return this.fromId;
	}

	/**
	 * @return the toId
	 */
	public SequenceId getToId()
	{
		return this.toId;
	}

	/**
	 * @return the score
	 */
	public double getScore()
	{
		return this.score;
	}

	@Override
	public int compareTo(MatchResult o)
	{
		return -Double.compare(this.score, o.score);
	}
	
	public int getAShift()
	{
		return this.a;
	}
	
	public int getBShift()
	{
		return this.b;
	}

	@Override
	public String toString()
	{
		return String.format("%s %s %s %d %d %.5f", getFromId().getHeaderId(), getToId().getHeaderId(),
				getFromId().isForward()&&getToId().isForward() ? 'N' : 'I', 
				getAShift(), getBShift(), (1.0-getScore())*100.0);
	}


}
