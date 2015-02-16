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
package edu.umd.marbl.mhap.sketch;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import edu.umd.marbl.mhap.utils.HitCounter;

public final class MinHashSketch implements Sketch<MinHashSketch>
{
	private final int[] minHashes;
	/**
	 * 
	 */
	private static final long serialVersionUID = 8846482698636860862L;
	
	public final static int[] computeNgramMinHashesWeighted(String seq, final int nGramSize, final int numHashes,
			HashSet<Long> filter, NGramCounts kmerCount, boolean weighted)
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
	
		int[] hashes = new int[Math.max(1,numHashes)];		
	
		long[] best = new long[numHashes];
		Arrays.fill(best, Long.MAX_VALUE);

		for (Entry<Long, HitCounter> kmer : hitMap.entrySet())
		{
			long key = kmer.getKey();
			int weight = kmer.getValue().count;
			
			if (kmerCount!=null || filter != null)
			{				
				if ((kmerCount!=null && kmerCount.documentFrequencyRatio(key)>kmerCount.getFilterCutoff()) || (filter != null && filter.contains(key)))
				{
					//System.err.println("Bad = "+kmerCount.inverseDocumentFrequency(key)+", "+kmerCount.weight(key, weight, maxCount));
					if (!weighted)
						weight = 0;
				}
				else
				{
					if (weighted)
						weight = weight*3;
					else
						weight = 1;
				}
			}
			//System.err.println("Good = "+kmerCount.inverseDocumentFrequency(key)+", "+kmerCount.weight(key, weight, maxCount));
			//int weight = Math.min(1, (int)Math.round(kmerCount.weight(key, kmer.getValue().count, maxCount)));
			
			if (weight<=0)
				continue;
		
			//set the initial shift value
			long x = key;
			for (int word = 0; word < numHashes; word++)
			{
				for (int count = 0; count<weight; count++)
				{				
					// XORShift Random Number Generators
					x ^= (x << 21);
					x ^= (x >>> 35);
					x ^= (x << 4);
		
					if (x < best[word])
					{
						best[word] = x;
						if (word%2==0)
							hashes[word] = (int)key;
						else
							hashes[word] = (int)key>>32;
					}
				}
			}
		}
		
		//now combine into super shingles
		/*
		HashFunction hf = Hashing.murmur3_32(0);
		
		int[] superShingles = new int[numHashes];
		for (int iter=0; iter<hashes.length; iter++)
		{
			int i1 = iter;
			int i2 = (iter+1)%numHashes;
			
			HashCode hc = hf.newHasher().
					putInt(hashes[i1]).
					putInt(hashes[i2]).
					hash();
			superShingles[iter] = hc.asInt();
		}
		hashes = superShingles;
		*/ 
		
		return hashes;
	}

	public static MinHashSketch fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			//store the size
			int hashNum = input.readInt();
			
			//store the array
			int[] minHashes = new int[hashNum];
			for (int hash=0; hash<hashNum; hash++)
			{
				minHashes[hash] = input.readInt();
			}
			
			return new MinHashSketch(minHashes);
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	private MinHashSketch(int[] minHashes)
	{
		this.minHashes = minHashes;
	}
	
	public MinHashSketch(String seq, int nGramSize, int numHashes, HashSet<Long> filter, NGramCounts kmerCount, boolean weighted)
	{
		this.minHashes = MinHashSketch.computeNgramMinHashesWeighted(seq, nGramSize, numHashes, filter, kmerCount, weighted);
	}
	
	public MinHashSketch(String str, int nGramSize, int numHashes)
	{
		this.minHashes = MinHashSketch.computeNgramMinHashesWeighted(str, nGramSize, numHashes, null, null, true);
	}

	public byte[] getAsByteArray()
	{
		ByteBuffer bb = ByteBuffer.allocate(4*(1+this.minHashes.length));
		
		//store the size
		bb.putInt(this.minHashes.length);
		
		//store the array
		for (int hash=0; hash<this.minHashes.length; hash++)
			bb.putInt(this.minHashes[hash]); 
    
		return bb.array();
	}
	
	/**
	 * @return the minHashes
	 */
	public final int[] getMinHashArray()
	{
		return this.minHashes;
	}

	public final double jaccard(MinHashSketch h)
	{
		int count = 0;
		int size = this.minHashes.length;
		
		if (h.minHashes.length!=size)
			throw new SketchRuntimeException("MinHashes must be of same length in order to be compared.");
		
		for (int iter=0; iter<size; iter++)
		{
			if (this.minHashes[iter]==h.minHashes[iter])
				count++;
		}
		
		return (double)count/(double)size;
	}
	
	public final int numHashes()
	{
		return this.minHashes.length;
	}
	
	@Override
	public double similarity(MinHashSketch sh)
	{
		return jaccard(sh);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "MinHash "+Arrays.toString(this.minHashes) + "";
	}
}
