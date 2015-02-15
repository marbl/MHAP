package edu.umd.marbl.mhap.impl;

public final class OverlapInfo
{
	public final double score;
	public final double rawScore;
	public final int a1;
	public final int b1;
	public final int a2;
	public final int b2;
	
	public OverlapInfo(double score, double rawScore, int a1, int a2, int b1, int b2)
	{
		this.score = score;
		this.rawScore = rawScore;
		this.a1 = a1;
		this.a2 = a2;
		this.b1 = b1;
		this.b2 = b2;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "[score="+this.score+", a1="+this.a1+" a2="+this.a2+", b1="+this.b1+" b2="+this.b2+"]";
	}
	
	public String toBlasrString()
	{
		return ""+this.score+", "+this.a1+", "+this.a2+", "+this.b1+", "+this.b2;
	}

}
