package com.secret.fastalign.general;

import java.util.ArrayList;
import java.util.Arrays;

import com.secret.fastalign.utils.Pair;

public class AbstractOrderedSequenceHashes<H extends AbstractSequenceHashes<H>, T extends AbstractOrderedSequenceHashes<H,T>>
{
	private final H mainHashes;
	private final ArrayList<H> orderedHashes;
	
	public AbstractOrderedSequenceHashes(H mainHashes, ArrayList<H> orderedHashes)
	{
		this.mainHashes = mainHashes;
		this.orderedHashes = orderedHashes;
	}
	
	public SequenceId getSequenceId()
	{
		return this.mainHashes.getSequenceId();
	}
	
	public int getSequenceLength()
	{
		return this.mainHashes.sequenceLength();
	}
	
	public H getMainHashes()
	{
		return this.mainHashes;
	}
	
	public ArrayList<H> getOrderedHashes()
	{
		return this.orderedHashes;
	}


	public double jaccord(T s)
	{
		return getMainHashes().jaccord(s.getMainHashes());
	}

	public Pair<Double,Integer> orderedScore(T s)
	{
		ArrayList<H> myOrderedHashes = getOrderedHashes();
		ArrayList<H> sOrderedHashes = s.getOrderedHashes();
		
		double[] scoreMatrix = new double[myOrderedHashes.size()*sOrderedHashes.size()];
		Arrays.fill(scoreMatrix, Double.NaN);
		
		double bestAlignmentScore = Double.NEGATIVE_INFINITY;
		int bestShift = 0;
		//int bestOverlapSize = 0;
		
		for (int shift=-myOrderedHashes.size()+1; shift<sOrderedHashes.size(); shift++)
		{
			double score = 0.0;
			
			int minIndex = Math.max(0, shift);
			int maxIndex = Math.min(sOrderedHashes.size(), myOrderedHashes.size()+shift);
			int overlapSize = maxIndex-minIndex;
			
			for (int toIndex = minIndex; toIndex<maxIndex; toIndex++)
			{
				double currIndexScore = 0.0;
				
				int fromIndex = toIndex-shift;
				H hash = myOrderedHashes.get(fromIndex);
				
				int matrixIndex = fromIndex*sOrderedHashes.size()+toIndex-1;
				if (toIndex-1>=0)
				{
					double val = scoreMatrix[matrixIndex];
					if (Double.isNaN(val))
					{
						val = hash.jaccord(sOrderedHashes.get(toIndex-1));
						scoreMatrix[matrixIndex] = val;
					}
					else
						val = scoreMatrix[matrixIndex];
					
					currIndexScore += val;
				}
	
				matrixIndex++;
				if (toIndex>=0)
				{
					double val = scoreMatrix[matrixIndex];
					if (Double.isNaN(val))
					{
						val = hash.jaccord(sOrderedHashes.get(toIndex));
						scoreMatrix[matrixIndex] = val;
					}
					else
						val = scoreMatrix[matrixIndex];
					
					currIndexScore += val;
				}
				
				//store the lowest score
				//score = Math.min(score, currIndexScore);
				score += currIndexScore;
			}
			
			score = score/(double)overlapSize;
			
			if (score>bestAlignmentScore)
			{
				bestAlignmentScore = score;
				bestShift = shift;
				//bestOverlapSize = overlapSize;
			}
		}
		
		//compute the average jaccord score
		double relativeScore = bestAlignmentScore;
				
		return new Pair<Double,Integer>(relativeScore, bestShift);
	}
}