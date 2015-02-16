package edu.umd.marbl.mhap.sketch;

import edu.umd.marbl.mhap.align.AlignElementSketch;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.align.Alignment;
import edu.umd.marbl.mhap.impl.OverlapInfo;


public abstract class AbstractSequenceSubSketches<T extends AbstractSequenceSubSketches<T,S>, S extends Sketch<S>>
{
	protected final S[] sequence;
	private final int stepSize;
	private final int seqLength;
	
	protected AbstractSequenceSubSketches(S[] sequence, int stepSize, int seqLength)
	{
		this.sequence = sequence;
		this.stepSize = stepSize;
		this.seqLength = seqLength;
	}
	
	public Alignment<AlignElementSketch<S>> localAlignOneSkip(Aligner<AlignElementSketch<S>> aligner, T b)
	{
		return aligner.localAlignOneSkip(new AlignElementSketch<S>(this.sequence), new AlignElementSketch<S>(b.sequence));
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementSketch<S>> aligner, T b)
	{
		Alignment<AlignElementSketch<S>> aligment = localAlignOneSkip(aligner, b);
		
		int a1 = Math.min(getSequenceLength(), aligment.getA1()*this.stepSize);
		int a2 = Math.min(this.getSequenceLength(), aligment.getA2()*this.stepSize);
		int b1 = Math.min(b.getSequenceLength(), aligment.getB1()*this.stepSize);
		int b2 = Math.min(b.getSequenceLength(), aligment.getB2()*this.stepSize);
		
		return new OverlapInfo(aligment.getScore()/100000.0, aligment.getScore(), a1, a2, b1, b2);
	}
	
	public int getSequenceLength()
	{
		return this.seqLength;
	}
	
	public int getStepSize()
	{
		return this.stepSize;
	}
	
	public int length()
	{
		return this.sequence.length;
	}
}
