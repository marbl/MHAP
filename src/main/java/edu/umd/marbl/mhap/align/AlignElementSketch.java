package edu.umd.marbl.mhap.align;

import edu.umd.marbl.mhap.impl.OverlapInfo;
import edu.umd.marbl.mhap.sketch.Sketch;

public final class AlignElementSketch<T extends Sketch<T>> implements AlignElement<AlignElementSketch<T>>
{
	private final T[] elements;
	private final int seqLength;
	private final int stepSize;
	
	public AlignElementSketch(T[] sketchArray, int stepSize, int seqLength)
	{
		this.elements = sketchArray;
		this.stepSize = stepSize;
		this.seqLength = seqLength;
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementSketch<T>> aligner, AlignElementSketch<T> b)
	{
		Alignment<AlignElementSketch<T>> aligment = localAlignOneSkip(aligner, b);
		
		int a1 = aligment.getA1();
		int a2 = aligment.getA2();
		int b1 = aligment.getB1();
		int b2 = aligment.getB2();
		
		a1 = Math.min(getSequenceLength()-1, aligment.getA1()*getStepSize());		
		if (a2>=length()-1)
			a2 = getSequenceLength()-1;
		else
			a2 = aligment.getA2()*getStepSize()+getStepSize();
			
		b1 = Math.min(b.getSequenceLength()-1, aligment.getB1()*b.getStepSize());		
		if (b2>=b.length()-1)
			b2 = b.getSequenceLength()-1;
		else
			b2 = aligment.getB2()*b.getStepSize()+b.getStepSize();
		
		return new OverlapInfo(aligment.getScore()/100000.0, aligment.getScore(), a1, a2, b1, b2);
	}

	public int getSequenceLength()
	{
		return this.seqLength;
	}

	public T getSketch(int index)
	{
		return this.elements[index];
	}
	
	public int getStepSize()
	{
		return this.stepSize;
	}
	
	@Override
	public int length()
	{
		return this.elements.length;
	}
	
	public Alignment<AlignElementSketch<T>> localAlignOneSkip(Aligner<AlignElementSketch<T>> aligner, AlignElementSketch<T> b)
	{
		return aligner.localAlignOneSkip(this, b);
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
