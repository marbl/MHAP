package edu.umd.marbl.mhap.general;

public final class OverlapInfo
{
	public final double score;
	public final int kmerCount;
	public final int a1;
	public final int b1;
	public final int a2;
	public final int b2;
	
	public OverlapInfo(double score, int kmerCount, int a1, int a2, int b1, int b2)
	{
		this.score = score;
		this.kmerCount = kmerCount;
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
}
