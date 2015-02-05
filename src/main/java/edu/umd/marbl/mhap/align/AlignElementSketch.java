package edu.umd.marbl.mhap.align;

import edu.umd.marbl.mhap.sketch.Sketch;

public class AlignElementSketch<T extends Sketch<T>> implements AlignElement<AlignElementSketch<T>>
{
	private final T[] elements;
	private final double simOffset;
	
	public AlignElementSketch(T[] sketchArray, double simOffset)
	{
		this.elements = sketchArray;
		this.simOffset = simOffset;
	}
	
	@Override
	public int length()
	{
		return this.elements.length;
	}

	@Override
	public double similarityScore(AlignElementSketch<T> e, int i, int j)
	{
		return this.elements[i].similarity(e.elements[j])+this.simOffset;
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

	@Override
	public double getSimOffset()
	{
		return this.simOffset;
	}

}
