package com.secret.fastalign.hash;

public final class MatchResult implements Comparable<MatchResult>
{
	private final SequenceId fromId;
	private final SequenceId toId;
	private final double score;
	
	public MatchResult(SequenceId fromId, SequenceId toId, double score)
	{
		this.fromId = fromId;
		this.toId = toId;
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
		return Double.compare(this.score, o.score);
	}
}
