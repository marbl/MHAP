package edu.umd.marbl.mhap.align;

import java.util.ArrayList;
import java.util.Collections;

import edu.umd.marbl.mhap.align.Alignment.Operation;

public final class Aligner<S extends AlignElement<S>>
{
	private final float gapOpen;
	private final float gapExtend;
	private final boolean storePath;
	private final float scoreOffset;
	
	public Aligner(boolean storePath, double gapOpen, double gapExtend, double scoreOffset)
	{
		this.gapOpen = (float)gapOpen;
		this.gapExtend = (float)gapExtend;
		this.storePath = storePath;
		this.scoreOffset = (float)scoreOffset;
	}
	
	/*
	public Alignment<S> localAlignSmithWater(S a, S b)
	{
		if (a.length()==0 && b.length()==0)
			return null;
		else
		if (a.length()==0 || b.length()==0)
			return null;
		
		float[][] scores = new float[a.length()+1][b.length()+1];
		
		for (int i=1; i<=a.length(); i++)
			for (int j=1; j<=b.length(); j++)
			{
				float hNext = scores[i-1][j-1]+Math.min(0.0f, (float)a.similarityScore(b, i-1, j-1));
				
				float hDeletion = scores[i-1][j]+this.gapOpen;				
				float hInsertion = scores[i][j-1]+this.gapOpen;
				
				//adjustments for end
				//if (i==a.length())
				//	hInsertion = hInsertion-this.gapOpen;
				//if (j==b.length())
				//	hDeletion = hDeletion-this.gapOpen;
				
				float value = Math.max(Math.max(Math.max(0.0f, hNext), hDeletion), hInsertion);
				
				scores[i][j] = value;
			}
		
		double bestValue = scores[a.length()-1][b.length()-1];
		double score = bestValue/(double)Math.max(a.length(), b.length());
		
		//if (a.length()<500)
		//	System.err.println(edu.umd.marbl.mhap.utils.Utils.toString(scores));
		
		if (storePath)
		{
			//figure out the path
			ArrayList<Alignment.Operation> backOperations = new ArrayList<>(a.length()+b.length());
		
		
			int i = a.length();
			int j = b.length();
			while (i>0 || j>0)
			{
				if (i==0)
				{
					backOperations.add(Operation.DELETE);
					j--;					
				}
				else
				if (j==0)
				{
					backOperations.add(Operation.INSERT);
					i--;					
				}
				else
				if (scores[i-1][j-1]>=scores[i-1][j] && scores[i-1][j-1]>=scores[i][j-1])
				{
					backOperations.add(Operation.MATCH);
					i--;
					j--;
				}
				else
				if (scores[i-1][j]>=scores[i-1][j-1])
				{
					backOperations.add(Operation.INSERT);
					i--;
				}
				else
				{
					backOperations.add(Operation.DELETE);
					j--;
				}
			}
			
			return new Alignment<S>(a, b, score, this.gapOpen, backOperations);
		}
		
		return new Alignment<S>(a, b, score, this.gapOpen, null);		
	}
	*/
	
	public Alignment<S> localAlignSmithWaterGotoh(S a, S b)
	{		
		float[][] D = new float[a.length()+1][b.length()+1];
		float[][] P = new float[a.length()+1][b.length()+1];
		float[][] Q = new float[a.length()+1][b.length()+1];
		
		for (int i=1; i<=a.length(); i++)
		{
			D[i][0] = 0.0f;
			P[i][0] = Float.NEGATIVE_INFINITY;
			Q[i][0] = Float.NEGATIVE_INFINITY;
		}
		for (int j=1; j<=b.length(); j++)
		{
			D[0][j] = 0.0f;
			P[0][j] = Float.NEGATIVE_INFINITY;
			Q[0][j] = Float.NEGATIVE_INFINITY;
		}
		
		float maxValue = 0.0f;
		int maxI = 0;
		int maxJ = 0;
		for (int i=1; i<=a.length(); i++) {
			for (int j=1; j<=b.length(); j++)
			{			
				P[i][j] = Math.max(D[i-1][j]+this.gapOpen, P[i-1][j]+this.gapExtend);
				Q[i][j] = Math.max(D[i][j-1]+this.gapOpen, Q[i][j-1]+this.gapExtend);
								
				float score = D[i-1][j-1]+(float)a.similarityScore(b, i-1, j-1)+this.scoreOffset;
				
				//compute the actual score
				D[i][j] = Math.max(score, Math.max(P[i][j], Q[i][j]));
				
				if (D[i][j] > maxValue) {
					maxValue = D[i][j];
					maxI = i;
					maxJ = j;
				}
			}
		}
		
		float score = maxValue;
				
		int a1 = 0;
		int a2 = maxI;
		int b1 = 0;
		int b2 = maxJ;

		if (storePath)
		{
			//figure out the path
			ArrayList<Alignment.Operation> backOperations = new ArrayList<>(a.length()+b.length());
			int i = a.length();
			while (i > maxI) {
				backOperations.add(Operation.DELETE);
				i--;
			}
			
			i = maxI;
			int j = maxJ;
			while (i>0 && j>0)
			{
				if ((P[i][j]>=Q[i][j] && P[i][j]==D[i][j]) || j==0)
				{
					backOperations.add(Operation.DELETE);
					i--;
				}
				else
				if (Q[i][j]==D[i][j] || i==0)
				{
					backOperations.add(Operation.INSERT);
					j--;
				}
				else
				{
					backOperations.add(Operation.MATCH);
					i--;
					j--;
				}
			}
			a1 = i+1;
			b1 = j+1;
			while (i > 0) {
				backOperations.add(Operation.DELETE);
				i--;
			}
			
			//reverse the direction
			Collections.reverse(backOperations);
		
			return new Alignment<S>(a, b, a1, a2, b1, b2, score, this.gapOpen, backOperations);
		}
		
		return new Alignment<S>(a, b, a1, a2, b1, b2, score, this.gapOpen, null);
	}
	
	public Alignment<S> localAlignOneSkip(S a, S b)
	{		
		float[][] D = new float[a.length()+1][b.length()+1];
		float[][] P = new float[a.length()+1][b.length()+1];
		float[][] S = new float[a.length()+1][b.length()+1];
		
		float maxValue = 0.0f;
		int maxI = 0;
		int maxJ = 0;
		for (int i=1; i<=a.length(); i++) {
			for (int j=1; j<=b.length(); j++)
			{	
				float sim = (float)a.similarityScore(b, i-1, j-1)+this.scoreOffset;
				
				P[i][j] = Math.max(D[i-1][j]+this.gapOpen, D[i][j-1]+this.gapOpen);
				D[i][j] = D[i-1][j-1]+sim;
				
				S[i][j] = Math.max(P[i][j], D[i][j]);
				if (i==a.length())
					S[i][j] = Math.max(S[i][j], P[i][j-1]+this.gapOpen);
				if (j==b.length())
					S[i][j] = Math.max(S[i][j], P[i-1][j]+this.gapOpen);
				
				
				if (S[i][j] > maxValue && (i==a.length() || j==b.length())) 
				{
					maxValue = S[i][j];
					maxI = i;
					maxJ = j;
				}
			}
		}
		
		float score = maxValue;
				
		int a1 = 0;
		int a2 = maxI;
		int b1 = 0;
		int b2 = maxJ;

		if (this.storePath)
		{
			//figure out the path
			ArrayList<Alignment.Operation> backOperations = new ArrayList<>(a.length()+b.length());

			int i = maxI;
			int j = maxJ;
			while (i>0 && j>0)
			{
				if (S[i][j]==D[i-1][j]+this.gapOpen)
				{
					backOperations.add(Operation.DELETE);
					i--;
				}
				else
				if (S[i][j]==D[i][j-1]+this.gapOpen)
				{
					backOperations.add(Operation.INSERT);
					j--;
				}
				else
				{
					backOperations.add(Operation.MATCH);
					i--;
					j--;
				}
			}
			
			a1 = i;
			b1 = j;
			while (i > 0) 
			{
				backOperations.add(Operation.DELETE);
				i--;
			}
			while (j > 0) 
			{
				backOperations.add(Operation.INSERT);
				j--;
			}
			
			//reverse the direction
			Collections.reverse(backOperations);
		
			return new Alignment<S>(a, b, a1, a2, b1, b2, score, this.gapOpen, backOperations);
		}
		
		return new Alignment<S>(a, b, a1, a2, b1, b2, score, this.gapOpen, null);
	}
}