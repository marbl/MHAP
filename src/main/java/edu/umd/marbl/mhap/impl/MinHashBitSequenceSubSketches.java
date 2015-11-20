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
import edu.umd.marbl.mhap.align.AlignElementDoubleSketch;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.sketch.MinHashBitSketch;
import edu.umd.marbl.mhap.sketch.MinHashSketch;

public final class MinHashBitSequenceSubSketches
{
	private final AlignElementDoubleSketch<MinHashBitSketch> alignmentSketch;
	
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
			int[] sketch = new MinHashSketch(seq.substring(currStart, end), nGramSize, numWords*64).getMinHashArray();
			
			sequence[iter] = new MinHashBitSketch(sketch);
			
			start += stepSize;
		}

		return sequence;
	}
	
	public final static MinHashBitSketch[] computeSequencesDouble(String seq, int nGramSize, int stepSize, int numWords)
	{
		int remainder = seq.length()%stepSize;
		
		//get number of sequence
		int numSequence = (seq.length()-remainder)/stepSize-1;
		
		//make sure big engough 
		if (remainder>=stepSize/2 && remainder>=nGramSize)
			numSequence++;
				
		//make sketches out of them
		int start = 0;		
		MinHashBitSketch[] sketches = new MinHashBitSketch[numSequence];
		for (int iter=0; iter<numSequence; iter++)
		{
			int end = Math.min(seq.length(), start+stepSize*2);
			int currStart = Math.max(0, end-stepSize*2);			

			//compute minhashes
			sketches[iter] = new MinHashBitSketch(new MinHashSketch(seq.substring(currStart, end), nGramSize, numWords*64).getMinHashArray());
			
			start += stepSize;
		}
		
		return sketches;
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
	
	/*
	private static int[] union(int[] minHashes1, int[] minHashes2)
	{
		int[] newHashes = new int[minHashes1.length]; 
		
		for (int iter=0; iter<newHashes.length; iter++)
			newHashes[iter] = Math.min(minHashes1[iter], minHashes2[iter]);
		
		return newHashes;
	}
	*/
	
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
