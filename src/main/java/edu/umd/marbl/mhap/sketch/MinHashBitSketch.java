package edu.umd.marbl.mhap.sketch;

public final class MinHashBitSketch extends AbstractBitSketch<MinHashBitSketch>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -44448450811302477L;

	private final static long[] getAsBits(int[] minHashes)
	{
		int numWords = minHashes.length/64;
		
		//now convert them to bits
		long[] bits = new long[numWords];
		
		//take only the last bit
		long mask = 0b1;
		
		int bitCount = 0;
		int wordCount = 0;
		for (int word = 0; word<numWords; word++)
		{
			long currWord = 0b0;
			
			for (int bit=0; bit<64; bit++)
			{
				currWord = (currWord << 1) | (minHashes[bitCount] & mask);				
						
				bitCount++;
			}
			
			bits[wordCount] = currWord;			
			wordCount++;
		}

		return bits;
	}
	
	public MinHashBitSketch(long[] bits)
	{
		super(bits);
	}
	
	public MinHashBitSketch(int[] minHashes)
	{
		super(getAsBits(minHashes));
	}
	
	public MinHashBitSketch(String seq, int nGramSize, int numWords)
	{
		super(getAsBits(new MinHashSketch(seq, nGramSize, numWords*64).getMinHashArray()));
	}
	
	public final double jaccard(final MinHashBitSketch sh)
	{
		int count = getIntersectionCount(sh);
		
		double sim = (double)count/(double) this.numberOfBits();
		double jaccard = (sim- 0.5) * 2.0;

		return Math.max(0.0, jaccard);
	}
}
