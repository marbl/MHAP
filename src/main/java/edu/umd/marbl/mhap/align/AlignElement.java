package edu.umd.marbl.mhap.align;

public interface AlignElement<S extends AlignElement<S>>
{
	public int length();
	
	public double similarityScore(S e, int i, int j);
	
	@Override
	public String toString();
	public String toString(int i);
	public String toString(S match, int i, int j);
}
