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
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.HitCounter;
import edu.umd.marbl.mhap.utils.MersenneTwisterFast;
import edu.umd.marbl.mhap.utils.Utils;

public final class MinHash implements Serializable
{
	private final int[] minHashes;
	private final int seqLength;
	/**
	 * 
	 */
	private static final long serialVersionUID = 8846482698636860862L;
	

	
	public final static int[] computeKmerMinHashesWeightedIntSuper(String seq, final int kmerSize, final int numHashes,
			HashSet<Long> filter, KmerCounts kmerCount, boolean weighted)
	{
		final int numberKmers = seq.length() - kmerSize + 1;
	
		if (numberKmers < 1)
			throw new MhapRuntimeException("Kmer size bigger than string length.");
	
		// get the kmer hashes
		final long[] kmerHashes = Utils.computeSequenceHashesLong(seq, kmerSize, 0);
		
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
						hashes[word] = (int)key;
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
	
	public final static int[] computeKmerMinHashes(String seq, final int kmerSize, final int numHashes,
			HashSet<Integer> filter)
	{
		if (numHashes % 2 != 0)
			throw new MhapRuntimeException("Number of words must be multiple of 2.");
	
		final int numberKmers = seq.length() - kmerSize + 1;
	
		if (numberKmers < 1)
			throw new MhapRuntimeException("Kmer size bigger than string length.");
	
		// get the rabin hashes
		final int[] kmerHashes = Utils.computeSequenceHashes(seq, kmerSize);
	
		int[] hashes = new int[Math.max(1,numHashes)];
		
		Arrays.fill(hashes, Integer.MAX_VALUE);
	
		int numWordsBy2 = numHashes / 2;
	
		// Random rand = new Random(0);
		for (int iter = 0; iter < kmerHashes.length; iter++)
		{
			// do not compute minhash for filtered data, keep Integer.MAX_VALUE
			if (filter != null && filter.contains(kmerHashes[iter]))
				continue;
	
			// set it in case requesting 0
			if (numHashes==0)
			{
				hashes[0] = kmerHashes[iter];
				continue;
			}
	
			long x = kmerHashes[iter];
			for (int word = 0; word < numWordsBy2; word++)
			{
				// hashes[iter][word] = rand.nextLong();
	
				// XORShift Random Number Generators
				x ^= (x << 21);
				x ^= (x >>> 35);
				x ^= (x << 4);
	
				int val1 = (int) x;
				int val2 = (int) (x >> 32);
	
				if (val1 < hashes[2 * word])
					hashes[2 * word] = val1;
	
				if (val2 < hashes[2 * word + 1])
					hashes[2 * word + 1] = val2;
			}
		}
	
		return hashes;
	}
	
	public static MinHash fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			//store the size
			//bb.putInt(this.seqLength);
			int seqLength = input.readInt();
			
			//bb.putInt(this.minHashes.length);
			int hashNum = input.readInt();
			
			//store the array
			int[] minHashes = new int[hashNum];
			for (int hash=0; hash<hashNum; hash++)
			{
				//bb.putInt(this.minHashes[seq][hash]);
				minHashes[hash] = input.readInt();
			}
			
			return new MinHash(seqLength, minHashes);
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	private MinHash(int seqLength, int[] minHashes)
	{
		this.seqLength = seqLength;
		this.minHashes = minHashes;
	}

	public MinHash(Sequence seq, int kmerSize, int numHashes, HashSet<Long> filter, KmerCounts kmerCount, boolean weighted)
	{
		this.seqLength = seq.length();

		//this.minHashes = MinHash.computeKmerMinHashes(seq.getString(), kmerSize, numHashes, filter);
		//this.minHashes = MinHash.computeKmerMinHashesWeighted(seq.getString(), kmerSize, numHashes, filter);
		//this.minHashes = MinHash.computeKmerMinHashesWeightedInt(seq.getString(), kmerSize, numHashes, filter, kmerCount);
		this.minHashes = MinHash.computeKmerMinHashesWeightedIntSuper(seq.getString(), kmerSize, numHashes, filter, kmerCount, weighted);
	}
	
	public MinHash(String str, int kmerSize, int numHashes)
	{
		this.seqLength = str.length();
		this.minHashes = MinHash.computeKmerMinHashesWeightedIntSuper(str, kmerSize, numHashes, null, null, true);
	}
	
	public byte[] getAsByteArray()
	{
		ByteBuffer bb = ByteBuffer.allocate(4*(2+this.minHashes.length));
		
		//store the size
		bb.putInt(this.seqLength);
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

	public final int getSequenceLength()
	{
		return this.seqLength;
	}
	
	public final double jaccard(MinHash h)
	{
		int count = 0;
		int size = this.minHashes.length;
		
		if (h.minHashes.length!=size)
			throw new MhapRuntimeException("MinHashes must be of same length in order to be comapred.");
		
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

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "MinHash "+Arrays.toString(this.minHashes) + "";
	}
	
	public final static int[] computeKmerMinHashesWeighted(String seq, final int kmerSize, final int numHashes,
			HashSet<Integer> filter, KmerCounts kmerCount)
	{
		final int numberKmers = seq.length() - kmerSize + 1;
	
		if (numberKmers < 1)
			throw new MhapRuntimeException("Kmer size bigger than string length.");
	
		// get the rabin hashes
		final int[] kmerHashes = Utils.computeSequenceHashes(seq, kmerSize);
		
		//now compute the counts of occurance
		HashMap<Integer, HitCounter> hitMap = new LinkedHashMap<>(kmerHashes.length);
		int maxCount = 0;
		for (int kmer : kmerHashes)
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
	
		MersenneTwisterFast rand = new MersenneTwisterFast();

		//init the data
		int[] hashes = new int[Math.max(1,numHashes)];		
		Arrays.fill(hashes, Integer.MAX_VALUE);
		
		double[] best = new double[numHashes];
		Arrays.fill(best, Double.MAX_VALUE);

		for (Entry<Integer, HitCounter> kmer : hitMap.entrySet())
		{
			int key = kmer.getKey();
			int weight = kmer.getValue().count;
			
			//set the seed
			rand.setSeed(key);
			
			for (int word = 0; word < numHashes; word++)
			{
				double val = rand.nextDouble();
				
				//transform value based on 2008 Near Duplicate Image Detection: min-Hash and tf-idf Weighting paper, eq 8
				val = -Math.log(val)/weight;
				
				if (val < best[word])
				{
					hashes[word] = key;
					best[word] = val;
				}
			}
		}
	
		return hashes;
	}
	
	public final static int[] computeKmerMinHashesWeightedInt(String seq, final int kmerSize, final int numHashes,
			HashSet<Integer> filter, KmerCounts kmerCount)
	{
		if (numHashes % 2 != 0)
			throw new MhapRuntimeException("Number of words must be multiple of 2.");
	
		final int numberKmers = seq.length() - kmerSize + 1;
	
		if (numberKmers < 1)
			throw new MhapRuntimeException("Kmer size bigger than string length.");
	
		// get the rabin hashes
		final int[] kmerHashes = Utils.computeSequenceHashes(seq, kmerSize);
		
		//now compute the counts of occurance
		HashMap<Integer, HitCounter> hitMap = new LinkedHashMap<>(kmerHashes.length);
		int maxCount = 0;
		for (int kmer : kmerHashes)
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
		int numWordsBy2 = numHashes / 2;
	
		int[] best1 = new int[numWordsBy2];
		int[] best2 = new int[numWordsBy2];		
		Arrays.fill(best1, Integer.MAX_VALUE);
		Arrays.fill(best2, Integer.MAX_VALUE);

		for (Entry<Integer, HitCounter> kmer : hitMap.entrySet())
		{
			int key = kmer.getKey();
			int weight = kmer.getValue().count;
			
			if (kmerCount.documentFrequencyRatio(key)>1.0e-5)
			{
				continue;
			}

			/*
			if (kmerCount.documentFrequencyRatio(key)>1.0e-5)
			{
				System.err.println("Bad = "+kmerCount.inverseDocumentFrequency(key)+", "+kmerCount.weight(key, weight, maxCount));								
				continue;
			}
			System.err.println("Good = "+kmerCount.inverseDocumentFrequency(key)+", "+kmerCount.weight(key, weight, maxCount));
			*/
			//int weight = Math.min(1, (int)Math.round(kmerCount.weight(key, kmer.getValue().count, maxCount)));

			// do not compute minhash for filtered data, keep Integer.MAX_VALUE
			if (filter != null && filter.contains(key))
				continue;
		
			long x = key;
			
			for (int word = 0; word < numWordsBy2; word++)
			{
				for (int count = 0; count<weight; count++)
				{				
					// XORShift Random Number Generators
					x ^= (x << 21);
					x ^= (x >>> 35);
					x ^= (x << 4);
		
					int val1 = (int) x;
					int val2 = (int) (x >> 32);
		
					if (val1 < best1[word])
					{
						best1[word] = val1;
						hashes[2 * word] = key;
					}
		
					if (val2 < best2[word])
					{
						best2[word] = val2;
						hashes[2 * word + 1] = key;
					}
				}
			}
		}
		
		return hashes;
	}
}
