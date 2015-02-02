package edu.umd.marbl.mhap.align;

public class AlignElementString implements AlignElement<AlignElementString>
{
	private final double EXACT_MATCH_SCORE = 1.0;
	
	private final double MISMATCH_SCORE = -1.0;
	private final char[] s;
	
	public AlignElementString(String s)
	{
		this.s = s.toCharArray();
	}
	
	@Override
	public int length()
	{
		return s.length;
	}

	@Override
	public double similarityScore(AlignElementString e, int i, int j)
	{
		if (this.s[i]==e.s[j])
			return EXACT_MATCH_SCORE;
		else
			return MISMATCH_SCORE;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return new String(s);
	}

	@Override
	public String toString(AlignElementString match, int i, int j)
	{
		return ""+s[i];
	}

	@Override
	public String toString(int i)
	{
		return ""+s[i];
	}
}
