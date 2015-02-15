package edu.umd.marbl.mhap.sketch;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;

public final class MinHashBitSequenceSubSketches extends AbstractSequenceSubSketches<MinHashBitSequenceSubSketches, MinHashBitSketch>
{
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
	
	public final static MinHashBitSequenceSubSketches fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			int numSubSequences = input.readInt();
			int numWords = input.readInt();
			int stepSize = input.readInt();			
						
			MinHashBitSketch[] sequence = new MinHashBitSketch[numSubSequences];
			
			for (int iter=0; iter<numSubSequences; iter++)
			{
				long[] bits = new long[numWords];
				for (int word=0; word<numWords; word++)
					bits[word] = input.readLong();
				
				sequence[iter] = new MinHashBitSketch(bits);
			}
			
			return new MinHashBitSequenceSubSketches(sequence, stepSize);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	protected MinHashBitSequenceSubSketches(MinHashBitSketch[] sketches, int stepSize)
	{
		super(sketches, stepSize);
	}
	
	public MinHashBitSequenceSubSketches(String seq, int kmerSize, int stepSize, int numWords)
	{
		super(computeSequences(seq, kmerSize, stepSize, numWords), stepSize);
	}
	
	public byte[] getAsByteArray()
	{
		int numWords = this.sequence[0].numberOfWords();
		
		ByteBuffer bb = ByteBuffer.allocate(8*numWords*this.sequence.length+4+4+4);
		
		//store the size
		bb.putInt(this.sequence.length);
		bb.putInt(numWords);
		bb.putInt(getStepSize());
		
		//store the array
		for (int hash=0; hash<this.sequence.length; hash++)
		{
			MinHashBitSketch sketch = this.sequence[hash];
			for (int word=0; word<numWords; word++)
				bb.putLong(sketch.getWord(word)); 
		}
    
		return bb.array();
	}
}
