/* 
 * MHAP package
 * 
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2014 by Konstantin Berlin and Sergey Koren
 * University Of Maryland
 * 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */
package edu.umd.marbl.mhap.main;
import jaligner.Alignment;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.LogManager;

import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.MatchResult;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.sketch.MinHashSearch;
import edu.umd.marbl.mhap.sketch.SequenceSketchStreamer;

public class AlignmentHashRun 
{	
	private static final int DEFAULT_NUM_WORDS = 1024;

	private static final int DEFAULT_KMER_SIZE = 16;
	
	protected static final int DEFAULT_NUM_MIN_MATCHES = 1;

	protected static final int DEFAULT_SUB_SEQUENCE_SIZE = 5000;

	private static final double DEFAULT_THRESHOLD = 0.00;

	public static void main(String[] args) throws Exception {
		String inFile = null;
		int kmerSize = DEFAULT_KMER_SIZE;
		double threshold = DEFAULT_THRESHOLD;
		int numWords = DEFAULT_NUM_WORDS; 
		int numThreads = Runtime.getRuntime().availableProcessors()*2;

		for (int i = 0; i < args.length; i++) {
			if (args[i].trim().equalsIgnoreCase("-k")) {
				kmerSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("-s")) {
				inFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("--num-hashes")) {
				numWords = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--threshold")) {
				threshold = Double.parseDouble(args[++i]);
			}
		}

		System.err.println("Running with input fasta: " + inFile);
		System.err.println("kmer size:\t" + kmerSize);
		System.err.println("threshold:\t" + threshold);
		System.err.println("num hashes\t" + numWords);

		LogManager.getLogManager().reset();
				
		// read and index the kmers
		long startTime = System.nanoTime();

		FastaData data = new FastaData(inFile, 0);
		data.enqueueFullFile();
			
		//SimHashSearch hashSearch = new SimHashSearch(kmerSize, numWords);
		SequenceSketchStreamer seqStreamer = new SequenceSketchStreamer(inFile, kmerSize, numWords, DEFAULT_SUB_SEQUENCE_SIZE, 12, null, null, true, 0);

		MinHashSearch hashSearch = new MinHashSearch(seqStreamer, numWords, DEFAULT_NUM_MIN_MATCHES, numThreads, true, 0, 800, DEFAULT_THRESHOLD);
		System.err.println("Processed "+seqStreamer.getNumberProcessed()+" sequences.");
		System.err.println("Time (s) to hash: " + (System.nanoTime() - startTime)*1.0e-9);
		
		System.out.println("Read in and processed "+seqStreamer.getNumberProcessed()+" sequences.");


		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		startTime = System.nanoTime();
		
		ArrayList<MatchResult> results = hashSearch.findMatches();
		
		System.err.println("Time (s) to score: " + (System.nanoTime() - startTime)*1.0e-9);
		
		//sort to get the best scores on top
		ArrayList<MatchResult> mixedResults = new ArrayList<MatchResult>();
		
		Collections.sort(results);		
		mixedResults.addAll(results.subList(0, Math.min(results.size(), 100)));
		//Collections.shuffle(results);
		mixedResults.addAll(results.subList(Math.max(0,results.size()-50), results.size()));
		
		//Collections.shuffle(mixedResults);

		System.out.println("Found "+results.size()+" matches:");
		
		Matrix matrix = MatrixLoader.load("/Users/kberlin/Dropbox/Projects/fast-align/src/test/resources/com/secret/fastalign/matrix/score_matrix.txt");
		
		//output result
		int count = 0;
		double mean = 0;
		for (MatchResult match : mixedResults)
		{
			//this already computes the reverse compliment
			Sequence s1 = data.getSequence(match.getFromId());
			Sequence s2 = data.getSequence(match.getToId());
			
			//compute the actual match
			double score = computeAlignment(s1, s2, matrix);
			
			//System.out.format("Sequence match (%s - %s) with identity score %f (SW=%f).\n", match.getFromId(), match.getToId(), match.getScore(), score);
			System.out.format("%s %f %f\n", match, match.getScore(), score);
			
			mean += match.getScore();
			
			count++;
			
			if (count>200)
				break;
		}
		
		mean = mean/count;
		
		System.out.println("Mean: "+mean);		
	}
	
	public static double computeAlignment(Sequence s1, Sequence s2, Matrix matrix)
	{		
		//compute the actual match
		Alignment alignment = jaligner.SmithWatermanGotoh.align(new jaligner.Sequence(s1.getString()), new jaligner.Sequence(s2.getString()), matrix, 5, 3);
		
		//double score = alignment.getScore();
		double score = getScoreWithNoTerminalGaps(alignment);
		
		//score = score/(double)Math.min(s1.length(), s2.length());
		
		return score;
	}
	
	public static float getScoreWithNoTerminalGaps(Alignment alignment)
	{
		char[] sequence1 = alignment.getSequence1();
		char[] sequence2 = alignment.getSequence2();
		char GAP = '-';
		//float extend = alignment.getExtend();
		//float open = alignment.getOpen();
		float extend = 1.0f;
		float open = 1.0f;
		//Matrix matrix = alignment.getMatrix();

		// The calculated score
		float calcScore = 0;

		// In the previous step there was a gap in the first sequence
		boolean previous1wasGap = false;

		// In the previous step there was a gap in the second sequence
		boolean previous2wasGap = false;

		int start = 0;
		int end = sequence1.length - 1;

		if (sequence1[start] == GAP)
		{
			while (sequence1[start] == GAP)
			{
				start++;
			}
		}
		else if (sequence2[start] == GAP)
		{
			while (sequence2[start] == GAP)
			{
				start++;
			}
		}

		if (sequence1[end] == GAP)
		{
			while (sequence1[end] == GAP)
			{
				end--;
			}
		}
		else if (sequence2[end] == GAP)
		{
			while (sequence2[end] == GAP)
			{
				end--;
			}
		}

		char c1, c2; // the next character
		for (int i = start; i <= end; i++)
		{
			c1 = sequence1[i];
			c2 = sequence2[i];
			// the next character in the first sequence is a gap
			if (c1 == GAP)
			{
				if (previous1wasGap)
				{
					calcScore -= extend;
					//calcScore -= 0;
				}
				else
				{
					calcScore -= open;
					//calcScore -= 0;
				}
				previous1wasGap = true;
				previous2wasGap = false;
			}
			// the next character in the second sequence is a gap
			else if (c2 == GAP)
			{
				if (previous2wasGap)
				{
					calcScore -= extend;
					//calcScore -= 0;
				}
				else
				{
					calcScore -= open;
					//calcScore -= 0;
				}
				previous1wasGap = false;
				previous2wasGap = true;
			}
			// the next characters in boths sequences are not gaps
			else
			{
				//calcScore += matrix.getScore(c1, c2);
				calcScore += 1;
				previous1wasGap = false;
				previous2wasGap = false;
			}
		}
		
		calcScore = calcScore/(float)(end-start+1);
		
		return calcScore;
	}
}
