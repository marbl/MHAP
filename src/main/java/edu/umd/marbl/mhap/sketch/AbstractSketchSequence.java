package edu.umd.marbl.mhap.sketch;

import edu.umd.marbl.mhap.align.AlignElementSketch;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.align.Alignment;
import edu.umd.marbl.mhap.impl.OverlapInfo;


public abstract class AbstractSketchSequence<T extends AbstractSketchSequence<T,S>, S extends Sketch<S>>
{
	protected final S[] sequence;
	private final int stepSize;
	
	protected AbstractSketchSequence(S[] sequence, int stepSize)
	{
		this.sequence = sequence;
		this.stepSize = stepSize;
	}
	
	public Alignment<AlignElementSketch<S>> customAlignSmithWaterGotoh(Aligner<AlignElementSketch<S>> aligner, T b)
	{
		return aligner.customAlignSmithWaterGotoh(new AlignElementSketch<S>(this.sequence), new AlignElementSketch<S>(b.sequence));
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementSketch<S>> aligner, T b)
	{
		Alignment<AlignElementSketch<S>> aligment = customAlignSmithWaterGotoh(aligner, b);
		
		return new OverlapInfo(aligment.getOverlapScore(), -1, aligment.getA1()*stepSize, aligment.getA2()*stepSize, aligment.getB1()*stepSize, aligment.getB2()*stepSize);
	}
	
	public int length()
	{
		return this.sequence.length;
	}
}
