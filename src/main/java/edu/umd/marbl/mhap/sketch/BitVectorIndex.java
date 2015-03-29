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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import edu.umd.marbl.mhap.utils.MersenneTwisterFast;
import edu.umd.marbl.mhap.utils.Pair;
import edu.umd.marbl.mhap.utils.SortablePair;

public final class BitVectorIndex<T,B extends AbstractBitSketch<B>>
{
	private final long bitsUsed[][];
	private final ArrayList<HashMap<Integer,ArrayList<Pair<T,B>>>> hashList;
	private final HashMap<T,B> indexedWords;
	private final double minSimilarity;
	
	public BitVectorIndex(List<Pair<T,B>> valuePairs, double minSimilarity, double confidence)
	{
		this.minSimilarity = minSimilarity;
		
		//should go off the valuePairs list
		int b = 10;
		
		//probability of a hit in numIndexes when using b: confidence = 1-(1-minSimilarity^b)^(numIndexes)
		//solve for b, Step 1: root_numIndexes (1-confidence) = (1-minSimilarity^b)
		//             Step 2: b = log(1-root_numIndexes (1-confidence))/log(minSimilarity)
		
		//figure out b
		int numIndexes = (int)Math.ceil(Math.log(1.0-confidence)/Math.log(1.0-Math.pow(this.minSimilarity, (double)b)));
		
		//allocate the memory
		this.bitsUsed = new long[numIndexes][b];
		
		//now generate random permuations
		MersenneTwisterFast rand = new MersenneTwisterFast();

		//get number of bits
		long numBits = 1;
		if (!valuePairs.isEmpty())
			numBits = valuePairs.get(0).y.numberOfBits();

		//generate the bits
		for (int index=0; index<numIndexes; index++)
			for (int bit=0; bit<b; bit++)
				this.bitsUsed[index][bit] = rand.nextLong(numBits);

		//allocate the memory
		this.hashList = new ArrayList<>(numIndexes);
		for (int iter=0; iter<numIndexes; iter++)
			this.hashList.add(new HashMap<>(valuePairs.size()));
		
		this.indexedWords = new HashMap<>(valuePairs.size());
		
		//encode all data in parallel
		valuePairs.parallelStream().forEach(pair-> {
			
			//get the lookup positions
			int[] lookupPositions = lookupPositions(pair.y);
			
			int count = 0;
			for(HashMap<Integer,ArrayList<Pair<T,B>>> map : this.hashList)
			{
				//get the array list
				ArrayList<Pair<T,B>> list;
				synchronized (map)
				{
					list = map.computeIfAbsent(lookupPositions[count], key-> new ArrayList<>(1));
				}
				
				//add the pair to the index
				synchronized(list)
				{
					list.add(pair);
				}
				
				count++;
			}
			
			//add the word
			synchronized (this.indexedWords)
			{
				this.indexedWords.put(pair.x, pair.y);
			}
		});
	}
	
	public int getBitsPerHash()
	{
		return bitsUsed[0].length;
	}
	
	public Map<T,B> getIndexedItems()
	{
		return Collections.unmodifiableMap(this.indexedWords);
	}
	
	public List<SortablePair<Double,T>> getNeighbors(B sketch, double minSimilarity)
	{		
		if (minSimilarity<this.minSimilarity)
			throw new SketchRuntimeException("Similarity request threshold below the ability of the indexer to compute.");
		
		int[] lookupPositions = lookupPositions(sketch);
		
		//now get a large hashset of items
		HashSet<Pair<T,B>> set = new HashSet<>();
		
		int count = 0;
		for(HashMap<Integer,ArrayList<Pair<T,B>>> map : this.hashList)
		{
			ArrayList<Pair<T,B>> list = map.get(lookupPositions[count]);
			
			if (list==null)
				continue;
					
			//add all the elements
			set.addAll(list);
			
			count++;
		}
		
		ArrayList<SortablePair<Double,T>> returnList = new ArrayList<SortablePair<Double,T>>();
		
		//now do direct compare
		for (Pair<T,B> pair : set)
		{
			double score = pair.y.similarity(sketch);
			
			if (score>=minSimilarity)
				returnList.add(new SortablePair<>(score, pair.x));
		}
		
		return returnList;
	}
	
	public int getNumberOfIndexes()
	{
		return this.hashList.size();
	}
	
	public B getSketch(T word)
	{
		return this.indexedWords.get(word);
	}
	
	public boolean isEmpty()
	{
		return indexedWords.isEmpty();
	}

	private int[] lookupPositions(B bits)
	{
		int numIndexes = hashList.size();
		
		int[] returnValues = new int[numIndexes];
		for (int index=0; index<numIndexes; index++)
		{
			long[] usedBits = bitsUsed[index];
			
			int val = 0b0;
			int mask = 0b1;
			for (int bit=0; bit<usedBits.length; bit++)
			{
				if (bits.getBit(usedBits[bit]))
					val = val | mask;
				
				mask = mask<<1;
			}
			
			returnValues[index] = val;
		}
		
		return returnValues;
	}
}
