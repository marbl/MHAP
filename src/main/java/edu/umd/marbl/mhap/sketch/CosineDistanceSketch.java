package edu.umd.marbl.mhap.sketch;

import edu.umd.marbl.mhap.math.BasicMath;

public final class CosineDistanceSketch extends AbstractBitSketch<CosineDistanceSketch>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -6501138603779963996L;

	private static long[] getCuts(double[] vector, int numWords, int seed)
	{
		long[] bitVector = new long[numWords];
		
		for (int word=0; word<numWords; word++)
		{
			long currBitLong = 0b0;
			
			long mask = 0b1;
			for (int bit=0; bit<64; bit++)
			{
				double[] rVec = HashUtils.randomGuassianVector(vector.length, seed+(word+1)*bit);
				double proj = BasicMath.dotProduct(vector, rVec);
				
				if (proj>0.0)
					currBitLong = currBitLong | mask;
				
				mask = mask<<1;
			}
			
			bitVector[word] = currBitLong;
		}
		
		return bitVector;
	}
	
	public CosineDistanceSketch(double[] vector, int numWords, int seed)
	{
		super(getCuts(vector,numWords,seed));
	}

}
