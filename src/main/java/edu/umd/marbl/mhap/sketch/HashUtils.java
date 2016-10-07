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
package edu.umd.marbl.mhap.sketch;

import java.nio.ByteBuffer;

import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;

import edu.umd.marbl.mhap.math.BasicMath;
import edu.umd.marbl.mhap.utils.MersenneTwisterFast;
import edu.umd.marbl.mhap.utils.Utils;

public class HashUtils
{
	public static long[] computeHashes(String item, int numWords, int seed)
	{
		long[] hashes = new long[numWords];
	
		for (int word = 0; word < numWords; word += 2)
		{
			HashFunction hashFunc = Hashing.murmur3_128(seed + word);
			Hasher hasher = hashFunc.newHasher();
			hasher.putUnencodedChars(item);
	
			// get the two longs out
			HashCode hc = hasher.hash();
			ByteBuffer bb = ByteBuffer.wrap(hc.asBytes());
			hashes[word] = bb.getLong(0);
			if (word + 1 < numWords)
				hashes[word + 1] = bb.getLong(8);
		}
	
		return hashes;
	}

	public final static int[] computeHashesInt(Object obj, int numWords, int seed)
	{
		if (obj instanceof Integer)
			return computeHashesIntInt((Integer) obj, numWords, seed);
		if (obj instanceof Long)
			return computeHashesIntLong((Long) obj, numWords, seed);
		if (obj instanceof Double)
			return computeHashesIntDouble((Double) obj, numWords, seed);
		if (obj instanceof Float)
			return computeHashesIntFloat((Float) obj, numWords, seed);
		if (obj instanceof String)
			return computeHashesIntString((String) obj, numWords, seed);
	
		throw new SketchRuntimeException("Cannot hash class type " + obj.getClass().getCanonicalName());
	}
	
	public final static int[] computeHashesIntDouble(double obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];
	
		HashFunction hf = Hashing.murmur3_32(seed);
	
		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putDouble(obj).putInt(iter).hash();
	
			hashes[iter] = hc.asInt();
		}
	
		return hashes;
	}

	public final static int[] computeHashesIntFloat(float obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];
	
		HashFunction hf = Hashing.murmur3_32(seed);
	
		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putFloat(obj).putInt(iter).hash();
	
			hashes[iter] = hc.asInt();
		}
	
		return hashes;
	}

	public final static int[] computeHashesIntInt(int obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];
	
		HashFunction hf = Hashing.murmur3_32(seed);
	
		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putInt(obj).putInt(iter).hash();
	
			hashes[iter] = hc.asInt();
		}
	
		return hashes;
	}

	public final static int[] computeHashesIntLong(long obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];
	
		HashFunction hf = Hashing.murmur3_32(seed);
	
		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putLong(obj).putInt(iter).hash();
	
			hashes[iter] = hc.asInt();
		}
	
		return hashes;
	}

	public final static int[] computeHashesIntString(String obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];
	
		HashFunction hf = Hashing.murmur3_32(seed);
	
		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putUnencodedChars(obj).putInt(iter).hash();
	
			hashes[iter] = hc.asInt();
		}
	
		return hashes;
	}

	public final static long[][] computeNGramHashes(final String seq, final int nGramSize, final int numWords, final int seed, boolean doReverseCompliment)
	{
		final int numberNGrams = seq.length()-nGramSize+1;
	
		if (numberNGrams < 1)
			throw new SketchRuntimeException("N-gram size bigger than string length.");
	
		// get the rabin hashes
		final long[] rabinHashes = computeSequenceHashesLong(seq, nGramSize, seed, doReverseCompliment);
	
		final long[][] hashes = new long[rabinHashes.length][numWords];
	
		// Random rand = new Random(0);
		for (int iter = 0; iter < rabinHashes.length; iter++)
		{
			// rand.setSeed(rabinHashes[iter]);
			long x = rabinHashes[iter];
	
			for (int word = 0; word < numWords; word++)
			{
				// hashes[iter][word] = rand.nextLong();
	
				// XORShift Random Number Generators
				x ^= (x << 21);
				x ^= (x >>> 35);
				x ^= (x << 4);
				hashes[iter][word] = x;
			}
		}
		
		return hashes;
	}

	public final static long[][] computeNGramHashesExact(final String seq, final int nGramSize, final int numWords, final int seed)
	{
		HashFunction hf = Hashing.murmur3_128(seed);
	
		long[][] hashes = new long[seq.length() - nGramSize + 1][numWords];
		for (int iter = 0; iter < hashes.length; iter++)
		{
			String subStr = seq.substring(iter, iter + nGramSize);
			
			for (int word=0; word<numWords; word++)
			{
				HashCode hc = hf.newHasher().putUnencodedChars(subStr).putInt(word).hash();
				hashes[iter][word] = hc.asLong();
			}
		}
		
		return hashes;
	}

	public final static int[] computeSequenceHashes(final String seq, final int nGramSize, boolean doReverseCompliment)
	{
		HashFunction hf = Hashing.murmur3_32(0);
	
		int[] hashes = new int[seq.length() - nGramSize + 1];
		for (int iter = 0; iter < hashes.length; iter++)
		{
			String str = seq.substring(iter, iter + nGramSize);
			
			String strReverse = null;
			if (doReverseCompliment)
			{
				strReverse  = Utils.rc(str);
				if (strReverse.compareTo(str)<0)
					str = strReverse;
			}

			HashCode hc = hf.newHasher().putUnencodedChars(str).hash();
			hashes[iter] = hc.asInt();
		}
	
		return hashes;
	}

	public final static long[] computeSequenceHashesLong(final String seq, final int nGramSize, final int seed, final boolean doReverseCompliment)
	{
		HashFunction hf = Hashing.murmur3_128(seed);
	
		long[] hashes = new long[seq.length() - nGramSize + 1];
		for (int iter = 0; iter < hashes.length; iter++)
		{
			String str = seq.substring(iter, iter + nGramSize);
			String strReverse = null;
			if (doReverseCompliment)
			{
				strReverse  = Utils.rc(str);
				if (strReverse.compareTo(str)<0)
					str = strReverse;
			}
			
			HashCode hc = hf.newHasher().putUnencodedChars(str).hash();
			hashes[iter] = hc.asLong();
		}
	
		return hashes;
	}
	
	public static double[] randomGuassianVector(int n, int seed)
	{
		//now generate the guassian
		MersenneTwisterFast rand = new MersenneTwisterFast(seed);
		
		double[] vec = new double[n];
		for (int iter=0; iter<n; iter++)
		{
			vec[iter] = rand.nextGaussian();
		}
		
		//normalize
		double norm = BasicMath.norm(vec);		
		if (norm<1.0e-10)
			return vec;
		
		return BasicMath.mult(vec, 1.0/norm);		
	}

	public static double[] randomStringGuassianVector(String str, int n, int seed)
	{
		int[] seeds = new int[4];
		for (int iter=0; iter<4; iter++)
		{
			HashFunction hf = Hashing.murmur3_32(seed*4+iter);
			HashCode hc = hf.newHasher().putUnencodedChars(str).hash();
			
			seeds[iter] = hc.asInt();
		}
		
		//now generate the guassian
		MersenneTwisterFast rand = new MersenneTwisterFast(seeds);
		
		double[] vec = new double[n];
		for (int iter=0; iter<n; iter++)
		{
			vec[iter] = rand.nextGaussian();
		}
		
		//normalize
		double norm = BasicMath.norm(vec);		
		if (norm<1.0e-10)
			return vec;
		
		return BasicMath.mult(vec, 1.0/norm);
	}

	private HashUtils()
	{
	}

}
