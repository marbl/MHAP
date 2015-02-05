package edu.umd.marbl.mhap.sketch;

import edu.umd.marbl.mhap.align.AlignElementSketch;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.align.Alignment;
import edu.umd.marbl.mhap.impl.OverlapInfo;


public abstract class AbstractSketchSequence<T extends AbstractSketchSequence<T,S>, S extends Sketch<S>>
{
	protected final S[] sequence;
	private final int stepSize;
	private final double simOffset;
	private final int minOverlap;
	
	protected AbstractSketchSequence(S[] sequence, double simOffset, int stepSize, int minOverlap)
	{
		this.sequence = sequence;
		this.stepSize = stepSize;
		this.simOffset = simOffset;
		this.minOverlap = minOverlap;
	}
	
	public Alignment<AlignElementSketch<S>> localAlignSmithWaterGotoh(Aligner<AlignElementSketch<S>> aligner, T b)
	{
		return aligner.localAlignSmithWaterGotoh(new AlignElementSketch<S>(this.sequence, this.simOffset), new AlignElementSketch<S>(b.sequence, this.simOffset));
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementSketch<S>> aligner, T b)
	{
		Alignment<AlignElementSketch<S>> aligment = localAlignSmithWaterGotoh(aligner, b);
		
		return new OverlapInfo(aligment.getOverlapScore(Math.max(1, this.minOverlap/stepSize)), -1, aligment.getA1()*stepSize, aligment.getA2()*stepSize, aligment.getB1()*stepSize, aligment.getB2()*stepSize);
	}
	
	public int length()
	{
		return this.sequence.length;
	}
}
