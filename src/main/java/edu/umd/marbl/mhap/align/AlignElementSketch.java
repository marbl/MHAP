package edu.umd.marbl.mhap.align;

import edu.umd.marbl.mhap.sketch.Sketch;

public class AlignElementSketch<T extends Sketch<T>> implements AlignElement<AlignElementSketch<T>>
{
	private final T[] elements; 
	
	public AlignElementSketch(T[] sketchArray)
	{
		this.elements = sketchArray;
	}
	
	@Override
	public int length()
	{
		return this.elements.length;
	}

	@Override
	public double similarityScore(AlignElementSketch<T> e, int i, int j)
	{
		return this.elements[i].similarity(e.elements[j]);
	}

	@Override
	public String toString(AlignElementSketch<T> match, int i, int j)
	{
		return toString();
	}

	@Override
	public String toString(int i)
	{
		return this.elements[i].toString();
	}

}
