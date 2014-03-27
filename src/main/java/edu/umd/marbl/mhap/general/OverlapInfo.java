package edu.umd.marbl.mhap.general;

public final class OverlapInfo
{
	public final double score;
	public final int a;
	public final int b;
	
	public OverlapInfo(double score, int a, int b)
	{
		this.score = score;
		this.a = a;
		this.b = b;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "[score="+this.score+", a="+this.a+", b="+this.b+"]";
	}
}
