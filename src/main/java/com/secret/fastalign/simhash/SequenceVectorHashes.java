package com.secret.fastalign.simhash;

import java.util.ArrayList;
import java.util.Arrays;

import com.secret.fastalign.data.SequenceId;
import com.secret.fastalign.utils.Pair;

public final class SequenceVectorHashes<T extends VectorHash<T>>
{
	private final int length;
	private final T mainHash;
	private final ArrayList<? extends T> subHashes;
	
	protected SequenceVectorHashes(int length, T mainHash, ArrayList<? extends T> subHashes)
	{
		this.length = length;
		this.mainHash = mainHash;
		this.subHashes = subHashes;
	}
	
	public SequenceId getSequenceId()
	{
		return this.mainHash.getSequenceId();
	}
	
	public double jaccord(SequenceVectorHashes<T> s)
	{
		return this.mainHash.jaccord(s.mainHash);
	}
	
	public Pair<Double,Integer> orderedScore(SequenceVectorHashes<T> s)
	{
		double[] scoreMatrix = new double[this.subHashes.size()*s.subHashes.size()];
		Arrays.fill(scoreMatrix, Double.NaN);
		
		double bestAlignmentScore = 0.0;
		int bestShift = 0;
		int overlapSize = 0;
		
		for (int shift=-this.subHashes.size(); shift<this.subHashes.size(); shift++)
		{
			double score = 0.0;
			
			int minIndex = Math.max(0, shift);
			int maxIndex = Math.min(s.subHashes.size(), this.subHashes.size()+shift);
			
			for (int toIndex = minIndex; toIndex<maxIndex; toIndex++)
			{
				int fromIndex = toIndex-shift;
				T hash = this.subHashes.get(fromIndex);
				
				int matrixIndex = fromIndex*s.subHashes.size()+toIndex-1;
				if (toIndex-1>=0)
				{
					double val = scoreMatrix[matrixIndex];
					if (Double.isNaN(val))
					{
						val = hash.jaccord(s.subHashes.get(toIndex-1));
						scoreMatrix[matrixIndex] = val;
					}
					else
						val = scoreMatrix[matrixIndex];
					
					score += val;
				}

				matrixIndex++;
				if (toIndex>=0)
				{
					double val = scoreMatrix[matrixIndex];
					if (Double.isNaN(val))
					{
						val = hash.jaccord(s.subHashes.get(toIndex));
						scoreMatrix[matrixIndex] = val;
					}
					else
						val = scoreMatrix[matrixIndex];
					
					score += val;
				}
	
				matrixIndex++;
				if (toIndex+1<s.subHashes.size())
				{
					double val = scoreMatrix[matrixIndex];
					if (Double.isNaN(val))
					{
						val = hash.jaccord(s.subHashes.get(toIndex+1));
						scoreMatrix[matrixIndex] = val;
					}
					else
						val = scoreMatrix[matrixIndex];
					
					score += val;
				}
			}
			
			if (score>bestAlignmentScore)
			{
				bestAlignmentScore = score;
				bestShift = shift;
				overlapSize = maxIndex-minIndex;
			}
		}
		
		//compute the average jaccord score
		double relativeScore = bestAlignmentScore/(double)overlapSize;
		
		return new Pair<Double,Integer>(relativeScore, bestShift);
	}
	
	public int sequenceLength()
	{
		return this.length;
	}
}
