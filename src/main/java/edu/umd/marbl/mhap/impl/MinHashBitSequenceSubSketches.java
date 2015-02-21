package edu.umd.marbl.mhap.impl;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;

import edu.umd.marbl.mhap.align.AlignElementSketch;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.sketch.MinHashBitSketch;

public final class MinHashBitSequenceSubSketches
{
	private final AlignElementSketch<MinHashBitSketch> alignmentSketch;
	
	private final static MinHashBitSketch[] computeSequences(String seq, int kmerSize, int stepSize, int numWords)
	{
		int remainder = seq.length()%stepSize;
		
		//get number of sequence
		int numSequence = seq.length()/stepSize;
		
		if (remainder>=kmerSize)
			numSequence++;
				
		//make sketches out of them
		MinHashBitSketch[] sequence = new MinHashBitSketch[numSequence];
		int start = 0;
		for (int iter=0; iter<sequence.length; iter++)
		{
			int end = Math.min(seq.length(), start+stepSize);
			
			sequence[iter] = new MinHashBitSketch(seq.substring(start, end), kmerSize, numWords);
			start += stepSize;
		}

		return sequence;
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementSketch<MinHashBitSketch>> aligner, MinHashBitSequenceSubSketches b)
	{
		return this.alignmentSketch.getOverlapInfo(aligner, b.alignmentSketch);
	}
	
	public final static MinHashBitSequenceSubSketches fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			int numSketches = input.readInt();
			int numWordsPerSketch = input.readInt();
			int stepSize = input.readInt();		
			int seqLength = input.readInt();
						
			MinHashBitSketch[] sequence = new MinHashBitSketch[numSketches];
			
			for (int iter=0; iter<numSketches; iter++)
			{
				long[] bits = new long[numWordsPerSketch];
				for (int word=0; word<numWordsPerSketch; word++)
					bits[word] = input.readLong();
				
				sequence[iter] = new MinHashBitSketch(bits);
			}
			
			return new MinHashBitSequenceSubSketches(sequence, stepSize, seqLength);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	protected MinHashBitSequenceSubSketches(MinHashBitSketch[] sketches, int stepSize, int seqLength)
	{
		this.alignmentSketch = new AlignElementSketch<>(sketches, stepSize, seqLength);
	}
	
	public MinHashBitSequenceSubSketches(String seq, int kmerSize, int stepSize, int numWords)
	{
		this.alignmentSketch = new AlignElementSketch<>(computeSequences(seq, kmerSize, stepSize, numWords), stepSize, seq.length());
	}
	
	public byte[] getAsByteArray()
	{
		int numSketches = this.alignmentSketch.length();
		int numWordsPerSketch = this.alignmentSketch.getSketch(0).numberOfWords();
		
		ByteBuffer bb = ByteBuffer.allocate(8*numWordsPerSketch*numSketches+4*4);
		
		//store the size
		bb.putInt(numSketches);
		bb.putInt(numWordsPerSketch);
		bb.putInt(this.alignmentSketch.getStepSize());
		bb.putInt(this.alignmentSketch.getSequenceLength());
		
		//store the array
		for (int sketchIndex=0; sketchIndex<numSketches; sketchIndex++)
		{
			MinHashBitSketch sketch = this.alignmentSketch.getSketch(sketchIndex);
			for (int word=0; word<numWordsPerSketch; word++)
				bb.putLong(sketch.getWord(word)); 
		}
    
		return bb.array();
	}
}
