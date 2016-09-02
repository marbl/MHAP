package edu.umd.marbl.mhap.sketch;

import it.unimi.dsi.fastutil.ints.IntArrays;

public class BottomSketch implements Sketch<BottomSketch>
{
	private final int[] hashPositions;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 9035607728472270206L;

	public BottomSketch(String str, int nGramSize, int k)
	{
		int[] hashes = HashUtils.computeSequenceHashes(str, nGramSize);
		
		k = Math.min(k, hashes.length);
		
		int[] perm = new int[hashes.length];
		for (int iter=0; iter<hashes.length; iter++)
			perm[iter] = iter;

		//sort the array
		IntArrays.radixSortIndirect(perm, hashes, true);
		
		hashPositions = new int[k];
		
		for (int iter=0; iter<k; iter++)
		{
			int index = perm[iter];
			hashPositions[iter] = hashes[index]; 
		}

	}
	
	public double jaccard(BottomSketch sh)
	{
		//make sure you look at same number
		int k = Math.min(this.hashPositions.length, sh.hashPositions.length);
		
		int i = 0;
		int j = 0;
		int intersectCount = 0;
		int unionCount = 0;
		while (unionCount<k)
		{
			if (this.hashPositions[i]<sh.hashPositions[j])
				i++;
			else
			if (this.hashPositions[i]>sh.hashPositions[j])
				j++;
			else
			{
				intersectCount++;
				i++;
				j++;
			}
			
			unionCount++;
		}
		
		return ((double)intersectCount)/(double)k;
	}

	@Override
	public double similarity(BottomSketch sh)
	{
		return jaccard(sh);
	}

}
