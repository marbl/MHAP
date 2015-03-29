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

public final class MinHashBitSketch extends AbstractBitSketch<MinHashBitSketch>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -44448450811302477L;

	private final static long[] getAsBits(int[] minHashes)
	{
		int numWords = minHashes.length/64;
		
		//now convert them to bits
		long[] bits = new long[numWords];
		
		//take only the last bit
		long mask = 0b1;
		
		int bitCount = 0;
		int wordCount = 0;
		for (int word = 0; word<numWords; word++)
		{
			long currWord = 0b0;
			
			for (int bit=0; bit<64; bit++)
			{
				currWord = (currWord << 1) | (minHashes[bitCount] & mask);				
						
				bitCount++;
			}
			
			bits[wordCount] = currWord;			
			wordCount++;
		}

		return bits;
	}
	
	public MinHashBitSketch(long[] bits)
	{
		super(bits);
	}
	
	public MinHashBitSketch(int[] minHashes)
	{
		super(getAsBits(minHashes));
	}
	
	public MinHashBitSketch(String seq, int nGramSize, int numWords)
	{
		super(getAsBits(new MinHashSketch(seq, nGramSize, numWords*64).getMinHashArray()));
	}
	
	public final double jaccard(final MinHashBitSketch sh)
	{
		int count = getIntersectionCount(sh);
		
		double sim = (double)count/(double) this.numberOfBits();
		double jaccard = (sim- 0.5) * 2.0;

		return Math.max(0.0, jaccard);
	}
}
