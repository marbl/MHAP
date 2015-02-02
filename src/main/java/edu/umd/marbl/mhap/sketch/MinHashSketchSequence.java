package edu.umd.marbl.mhap.sketch;

public final class MinHashSketchSequence extends AbstractSketchSequence<MinHashSketchSequence, MinHashBitSketch>
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
	
	public MinHashSketchSequence(String seq, int kmerSize, int stepSize, int numWords)
	{
		super(computeSequences(seq, kmerSize, stepSize, numWords), stepSize);
	}
}
