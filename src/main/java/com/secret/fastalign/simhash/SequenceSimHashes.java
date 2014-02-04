package com.secret.fastalign.simhash;

import java.util.ArrayList;
import java.util.Arrays;

import com.secret.fastalign.general.AbstractReducedSequence;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.Pair;

public final class SequenceSimHashes extends AbstractReducedSequence<SimHash,SequenceSimHashes>
{
	private final ArrayList<SimHash> orderedHashes;
	
	public SequenceSimHashes(Sequence seq, int kmerSize, int numWords, int subStringSize, int subKmerSize, int subWordSize)
	{
		super(seq.getId(), new SimHash(seq, kmerSize, numWords));
		this.orderedHashes = generateSubHashes(seq, subStringSize, subKmerSize, subWordSize);
	}
	
	private static ArrayList<SimHash> generateSubHashes(Sequence seq, int subStringSize, int subKmerSize, int subWordSize)
	{
		//generate the array of simhashes
		ArrayList<SimHash> subHashes = new ArrayList<SimHash>(seq.length()-subKmerSize+1);
		for (int iter=0; iter<seq.length()-subKmerSize+1; iter+=subStringSize)
		{
			String subString = seq.getString().substring(iter, Math.min(iter+subStringSize, seq.length()));
			Sequence subSequence = new Sequence(subString, seq.getId());
			
			subHashes.add(new SimHash(subSequence, subKmerSize, subWordSize));
		}
		
		return subHashes;
	}
	
	private ArrayList<SimHash> getOrderedHashes()
	{
		return this.orderedHashes;
	}
	
	public Pair<Double,Integer> orderedScore(SequenceSimHashes s)
	{
		ArrayList<SimHash> myOrderedHashes = getOrderedHashes();
		ArrayList<SimHash> sOrderedHashes = s.getOrderedHashes();
		
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
				SimHash hash = myOrderedHashes.get(fromIndex);
				
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
