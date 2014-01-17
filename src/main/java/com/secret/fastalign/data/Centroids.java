package com.secret.fastalign.data;

import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import com.secret.fastalign.utils.Pair;

public final class Centroids
{
  public static int hammingDistance(String s1, String s2) 
  {

    // check preconditions
    if (s1 == null || s2 == null || s1.length() != s2.length()) 
    {
       throw new IllegalArgumentException();
    }

    int distance = 0;
    for (int iter = 0; iter < s1.length(); iter++) 
    {
      if (s1.charAt(iter) != s2.charAt(iter)) 
      {
        distance++;
      }
    }
    
    return distance;

  }

	private static String generateRandomString(Random generator, int size)
	{
		char[] s = new char[size];
		
		//ACGT
		
		//generate the random sequence
		for (int iter=0; iter<s.length; iter++)
		{
			switch (generator.nextInt(4))
			{
				case 0: s[iter] = 'A'; break;
				case 1: s[iter] = 'C'; break;
				case 2: s[iter] = 'G'; break;
				case 3: s[iter] = 'T'; break;
			}
		}
		
		return new String(s);		
	}
	
	private final HashMap<Pair<Integer,String>,String> lookup;
	private final ArrayList<ArrayList<String>> centroids;
	private final Matrix matrix;
		
	public Centroids(int kmerSize, int num) throws MatrixLoaderException
	{
		Random generator = new Random(1);
		
		this.matrix = MatrixLoader.load("/Users/kberlin/Dropbox/Projects/fast-align/src/test/resources/com/secret/fastalign/matrix/score_matrix.txt");

		//generate centroids
		this.centroids = new ArrayList<ArrayList<String>>();
		int n = (int)(Math.pow(4, kmerSize)/Math.pow(4, kmerSize/8));
		
		for (int index=0; index<num; index++)
		{
			ArrayList<String> currCentroids = new ArrayList<String>();
			this.centroids.add(currCentroids);
			for (int iter=0; iter<n; iter++)
				currCentroids.add(generateRandomString(generator, kmerSize));
		}
		
		this.lookup = new HashMap<Pair<Integer,String>,String>();
	}
	
	public String findCentroid(int index, String s)
	{
		Pair<Integer,String> pair = new Pair<Integer,String>(index, s);
		String centroid = this.lookup.get(pair);
		
		if (centroid==null)
		{
			double bestScore = Double.NEGATIVE_INFINITY;
			
			//go through all the centroids and align
			//compute the actual match
			for (String currCenter : this.centroids.get(index))
			{
				//double score = jaligner.SmithWatermanGotoh.align(new jaligner.Sequence(s), new jaligner.Sequence(currCenter), this.matrix, 5, 5).getScore();
				double score = -hammingDistance(s, currCenter);
	
				if (score>bestScore)
				{
					centroid = currCenter;
					bestScore = score;
				}
			}
			
			//hash the result
			this.lookup.put(pair, centroid);
		}
		
		return centroid;
	}
}
