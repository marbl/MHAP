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
 * Copyright (c) 2015 by Konstantin Berlin and Sergey Koren
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
package edu.umd.marbl.mhap.impl;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import edu.umd.marbl.mhap.align.AlignElementDoubleSketch;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.sketch.HashUtils;
import edu.umd.marbl.mhap.sketch.MinHashBitSketch;
import edu.umd.marbl.mhap.sketch.SketchRuntimeException;
import edu.umd.marbl.mhap.utils.HitCounter;

public final class MinHashBitSequenceSubSketches
{
	private final AlignElementDoubleSketch<MinHashBitSketch> alignmentSketch;
	
	private final static int[] computeNgramMinHashesWeighted(String seq, final int nGramSize, final int numHashes)
	{
		final int numberNGrams = seq.length() - nGramSize + 1;
	
		if (numberNGrams < 1)
			throw new SketchRuntimeException("N-gram size bigger than string length.");
	
		// get the kmer hashes
		final long[] kmerHashes = HashUtils.computeSequenceHashesLong(seq, nGramSize, 0);
		
		//now compute the counts of occurance
		HashMap<Long, HitCounter> hitMap = new LinkedHashMap<>(kmerHashes.length);
		int maxCount = 0;
		for (long kmer : kmerHashes)
		{
			HitCounter counter = hitMap.get(kmer);
			if (counter==null)
			{
				counter = new HitCounter(1);
				hitMap.put(kmer, counter);
			}
			else
				counter.addHit();

			if (maxCount<counter.count)
				maxCount = counter.count;
		}
	
		int[] best = new int[numHashes];
		Arrays.fill(best, Integer.MAX_VALUE);

		for (Entry<Long, HitCounter> kmer : hitMap.entrySet())
		{
			long key = kmer.getKey();
			int weight = kmer.getValue().count;
					
			//set the initial shift value
			int x = (int)key;
			for (int word = 0; word < numHashes; word++)
			{
				for (int count = 0; count<weight; count++)
				{				
					// XORShift Random Number Generators
					x ^= (x << 21);
					x ^= (x >>> 35);
					x ^= (x << 4);
					
					int intX = (int)x; 
		
					if (intX < best[word])
						best[word] = intX;
				}
			}
		}
		
		return best;
	}
	
	public final static MinHashBitSketch[] computeSequences(String seq, int nGramSize, int stepSize, int numWords)
	{
		int remainder = seq.length()%stepSize;
		
		//get number of sequence
		int numSequence = (seq.length()-remainder)/stepSize;
		
		if (remainder>0)
			numSequence++;
				
		//make sketches out of them
		int start = 0;		
		MinHashBitSketch[] sequence = new MinHashBitSketch[numSequence];
		for (int iter=0; iter<numSequence; iter++)
		{
			int end = Math.min(seq.length(), start+stepSize);
			int currStart = Math.max(0, end-stepSize);			

			//compute minhashes
			int[] sketch = computeNgramMinHashesWeighted(seq.substring(currStart, end), nGramSize, numWords*64);
			
			sequence[iter] = new MinHashBitSketch(sketch);
			
			start += stepSize;
		}

		return sequence;
	}
	
	public final static MinHashBitSketch[] computeSequencesDouble(String seq, int nGramSize, int stepSize, int numWords)
	{
		int remainder = seq.length()%stepSize;
		
		//get number of sequence
		int numSequence = (seq.length()-remainder)/stepSize;
		
		if (remainder>0)
			numSequence++;
				
		//make sketches out of them
		int start = 0;		
		int[][] sketches = new int[numSequence][numWords*64];
		for (int iter=0; iter<numSequence; iter++)
		{
			int end = Math.min(seq.length(), start+stepSize);
			int currStart = Math.max(0, end-stepSize);			

			//compute minhashes
			sketches[iter] = computeNgramMinHashesWeighted(seq.substring(currStart, end), nGramSize, numWords*64);
			
			start += stepSize;
		}
		
		MinHashBitSketch[] sequence = new MinHashBitSketch[numSequence];
		for (int iter=0; iter<sketches.length; iter++)
		{
			//now convert in sequence double the length
			if ((iter+1)<sketches.length)
			{
				sequence[iter] = new MinHashBitSketch(union(sketches[iter], sketches[iter+1]));
				if ((iter+2)<sketches.length)
					sequence[iter+1] = new MinHashBitSketch(union(sketches[iter+1], sketches[iter+2]));
				else
					sequence[iter+1] = new MinHashBitSketch(sketches[iter+1]);
			}
			else
				sequence[iter] = new MinHashBitSketch(sketches[iter]);			

		}

		return sequence;
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementDoubleSketch<MinHashBitSketch>> aligner, MinHashBitSequenceSubSketches b)
	{
		return this.alignmentSketch.getOverlapInfo(aligner, b.alignmentSketch);
	}
	
	public final static MinHashBitSequenceSubSketches fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			int numSketches = input.readInt();
			int numWordsPerSketch = input.readInt();
			int stepSize = input.readInt();		
			int seqLength = input.readInt();
						
			MinHashBitSketch[] sequence = new MinHashBitSketch[numSketches];
			
			for (int iter=0; iter<numSketches; iter++)
			{
				long[] bits = new long[numWordsPerSketch];
				for (int word=0; word<numWordsPerSketch; word++)
					bits[word] = input.readLong();
				
				sequence[iter] = new MinHashBitSketch(bits);
			}
			
			return new MinHashBitSequenceSubSketches(sequence, stepSize, seqLength);
			
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	protected MinHashBitSequenceSubSketches(MinHashBitSketch[] sketches, int stepSize, int seqLength)
	{
		this.alignmentSketch = new AlignElementDoubleSketch<>(sketches, stepSize, seqLength);
	}
	
	public MinHashBitSequenceSubSketches(String seq, int kmerSize, int stepSize, int numWords)
	{
		this.alignmentSketch = new AlignElementDoubleSketch<>(computeSequencesDouble(seq, kmerSize, stepSize, numWords), stepSize, seq.length());
	}
	
	private static int[] union(int[] minHashes1, int[] minHashes2)
	{
		int[] newHashes = new int[minHashes1.length]; 
		
		for (int iter=0; iter<newHashes.length; iter++)
			newHashes[iter] = Math.min(minHashes1[iter], minHashes2[iter]);
		
		return newHashes;
	}
	
	public byte[] getAsByteArray()
	{
		int numSketches = this.alignmentSketch.length();
		int numWordsPerSketch = this.alignmentSketch.getSketch(0).numberOfWords();
		
		ByteBuffer bb = ByteBuffer.allocate(8*numWordsPerSketch*numSketches+4*4);
		
		//store the size
		bb.putInt(numSketches);
		bb.putInt(numWordsPerSketch);
		bb.putInt(this.alignmentSketch.getStepSize());
		bb.putInt(this.alignmentSketch.getSequenceLength());
		
		//store the array
		for (int sketchIndex=0; sketchIndex<numSketches; sketchIndex++)
		{
			MinHashBitSketch sketch = this.alignmentSketch.getSketch(sketchIndex);
			for (int word=0; word<numWordsPerSketch; word++)
				bb.putLong(sketch.getWord(word)); 
		}
    
		return bb.array();
	}
}
