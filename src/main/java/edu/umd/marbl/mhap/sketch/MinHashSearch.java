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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicLong;

import edu.umd.marbl.mhap.general.AbstractMatchSearch;
import edu.umd.marbl.mhap.general.MatchResult;
import edu.umd.marbl.mhap.general.OverlapInfo;
import edu.umd.marbl.mhap.general.SequenceId;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.HitCounter;

public final class MinHashSearch extends AbstractMatchSearch
{
	private final double acceptScore;

	private final ArrayList<Map<Integer, ArrayList<SequenceId>>> hashes;
	private final double maxShift;
	private final AtomicLong minhashSearchTime;
	private final AtomicLong sortMergeSearchTime;
	private final int minStoreLength;
	private final AtomicLong numberElementsProcessed;

	private final AtomicLong numberSequencesFullyCompared;
	private final AtomicLong numberSequencesHit;
	private final AtomicLong numberSequencesMinHashed;
	private final AtomicLong numberSubSequences;
	private final AtomicLong numberSubSequencesHit;

	private final int numMinMatches;
	private final HashMap<SequenceId, SequenceSketch> sequenceVectorsHash;

	
	public MinHashSearch(SequenceSketchStreamer data, int numHashes, int numMinMatches, int numThreads, 
			boolean storeResults, int minStoreLength, double maxShift, double acceptScore) throws IOException
	{
		super(numThreads, storeResults);

		this.minStoreLength = minStoreLength;
		this.numMinMatches = numMinMatches;
		this.maxShift = maxShift;
		this.acceptScore = acceptScore;
		this.numberSubSequencesHit = new AtomicLong();
		this.numberSequencesHit = new AtomicLong();
		this.numberSequencesFullyCompared = new AtomicLong();
		this.numberSubSequences = new AtomicLong();
		this.numberSequencesMinHashed = new AtomicLong();
		this.numberElementsProcessed = new AtomicLong();
		this.minhashSearchTime = new AtomicLong();
		this.sortMergeSearchTime = new AtomicLong();
		
		// enqueue full file, since have to know full size
		data.enqueueFullFile(false, this.numThreads);

		this.sequenceVectorsHash = new HashMap<SequenceId, SequenceSketch>(data.getNumberProcessed() + 100, (float) 0.75);

		this.hashes = new ArrayList<Map<Integer, ArrayList<SequenceId>>>(numHashes);
		for (int iter = 0; iter < numHashes; iter++)
		{
			Map<Integer,ArrayList<SequenceId>> map = new HashMap<Integer, ArrayList<SequenceId>>(data.getNumberSubSequencesProcessed()+100);
			
			this.hashes.add(map);
		}
		
		addData(data);
		
		System.err.println("Stored "+this.sequenceVectorsHash.size()+" sequences in the index.");
	}

	@Override
	public boolean addSequence(SequenceSketch currHash)
	{
		int[] currMinHashes = currHash.getMinHashes().getMinHashArray();

		if (currMinHashes.length != this.hashes.size())
			throw new MhapRuntimeException("Number of MinHashes of the sequence does not match current settings.");

		// put the result into the hashmap
		synchronized (this.sequenceVectorsHash)
		{
			SequenceSketch minHash = this.sequenceVectorsHash.put(currHash.getSequenceId(), currHash);
			if (minHash != null)
			{
				this.sequenceVectorsHash.put(currHash.getSequenceId(), minHash);

				throw new MhapRuntimeException("Sequence ID already exists in the hash table.");
			}			
		}
		
		// add the hashes
		int count = 0;
		SequenceId id = currHash.getSequenceId();
		for (Map<Integer, ArrayList<SequenceId>> hash : this.hashes)
		{
			ArrayList<SequenceId> currList;
			final int hashVal = currMinHashes[count];

			// get the list
			synchronized (hash)
			{
				currList = hash.get(hashVal);

				if (currList == null)
				{
					currList = new ArrayList<SequenceId>(2);
					hash.put(hashVal, currList);
				}
			}

			// add the element
			synchronized (currList)
			{
				currList.add(id);
			}
			
			count++;
		}

		//increment the subsequence counter 
		this.numberSubSequences.getAndIncrement();
		
		//increment the counter
		this.numberSequencesMinHashed.getAndIncrement();
		
		return true;
	}

	@Override
	public List<MatchResult> findMatches(SequenceSketch seqHashes, boolean toSelf)
	{
		//for performance reasons might need to change
		long startTime = System.nanoTime();

		MinHash minHash = seqHashes.getMinHashes();

		if (this.hashes.size() != minHash.numHashes())
			throw new MhapRuntimeException("Number of hashes does not match. Stored size " + this.hashes.size()
					+ ", input size " + minHash.numHashes() + ".");

		HashMap<SequenceId, HitCounter> bestSequenceHit = new HashMap<SequenceId, HitCounter>(this.numberSequencesMinHashed.intValue()/5+1);
		int[] minHashes = minHash.getMinHashArray();
		
		int hashIndex = 0;
		for (Map<Integer,ArrayList<SequenceId>> currHash : this.hashes)
		{
			ArrayList<SequenceId> currentHashMatchList = currHash.get(minHashes[hashIndex]);

			// if some matches exist add them
			if (currentHashMatchList != null)
			{
				this.numberElementsProcessed.getAndAdd(currentHashMatchList.size());

				for (SequenceId sequenceId : currentHashMatchList)
				{
					// get current count in the list
					HitCounter currentHitInfo = bestSequenceHit.get(sequenceId);

					// increment the count
					if (currentHitInfo == null)
					{
						currentHitInfo = new HitCounter(1);
						bestSequenceHit.put(sequenceId, currentHitInfo);
					}
					else
						currentHitInfo.addHit();
				}
			}
			
			hashIndex++;
		}
		
		//record the search time
		long minHashEndTime = System.nanoTime();
		this.minhashSearchTime.getAndAdd(minHashEndTime - startTime);
		
		//record number of hash matches processed
		this.numberSequencesHit.getAndAdd(bestSequenceHit.size());
		this.numberSubSequencesHit.getAndAdd(bestSequenceHit.size());
		
		// compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>(32);
		
		for (Entry<SequenceId, HitCounter> match : bestSequenceHit.entrySet())
		{
			//get the match id
			SequenceId matchId = match.getKey();
			
			// do not store matches with smaller ids, unless its coming from a short read
			if (toSelf && matchId.getHeaderId() == seqHashes.getSequenceId().getHeaderId())
				continue;

			//see if the hit number is high enough			
			if (match.getValue().count >= this.numMinMatches)
			{
				SequenceSketch matchedHashes = this.sequenceVectorsHash.get(match.getKey());
				if (matchedHashes==null)
					throw new MhapRuntimeException("Hashes not found for given id.");
				
				//never process short to short
				if (matchedHashes.getSequenceLength()<this.minStoreLength && seqHashes.getSequenceLength()<this.minStoreLength)
					continue;
				
				//never process long to long in self, with greater id
				if (toSelf 
						&& matchId.getHeaderId() > seqHashes.getSequenceId().getHeaderId()
						&& matchedHashes.getSequenceLength()>=this.minStoreLength
						&& seqHashes.getSequenceLength()>=this.minStoreLength)
					continue;
				
				//never do short to long
				if (toSelf 
						&& matchedHashes.getSequenceLength()<this.minStoreLength
						&& seqHashes.getSequenceLength()>=this.minStoreLength)
					continue;
				
				//compute the direct hash score
				OverlapInfo result = seqHashes.getOrderedHashes().getFullScore(matchedHashes.getOrderedHashes(), this.maxShift);
				
				//increment the counter
				this.numberSequencesFullyCompared.getAndIncrement();

				//if score is good add
				if (result.score >= this.acceptScore)
				{
					//OverlapInfo result2 = seqHashes.getOrderedHashes().getFullScoreExperimental(matchedHashes.getOrderedHashes(), this.maxShift);				
					//System.err.println(result.score+"  "+result2.score);

					MatchResult currResult = new MatchResult(seqHashes.getSequenceId(), matchId, result, seqHashes.getSequenceLength(), matchedHashes.getSequenceLength());

					// add to list
					matches.add(currResult);
				}
			}
		}
		
		//record the search time
		//TODO not clear why not working. Perhaps everything is too fast?
		long endTime = System.nanoTime();
		this.sortMergeSearchTime.getAndAdd(endTime-minHashEndTime);

		return matches;
	}

	public double getMinHashSearchTime()
	{
		return this.minhashSearchTime.longValue() * 1.0e-9;
	}
	
	public double getSortMergeTime()
	{
		return this.sortMergeSearchTime.longValue() * 1.0e-9;
	}


	public long getNumberElementsProcessed()
	{
		return this.numberElementsProcessed.get();
	}

	public long getNumberSequenceHashed()
	{
		return this.numberSequencesMinHashed.get();
	}

	public long getNumberSequencesFullyCompared()
	{
		return this.numberSequencesFullyCompared.get();
	}
	
	public long getNumberSequencesHit()
	{
		return this.numberSequencesHit.get();
	}
	
	public long getNumberSubSequencesHit()
	{
		return this.numberSubSequencesHit.get();
	}
		
	@Override
	public List<SequenceId> getStoredForwardSequenceIds()
	{
		ArrayList<SequenceId> seqIds = new ArrayList<SequenceId>(this.sequenceVectorsHash.size());
		for (SequenceSketch hashes : this.sequenceVectorsHash.values())
			if (hashes.getSequenceId().isForward())
				seqIds.add(hashes.getSequenceId());
		
		return seqIds;
	}

	@Override
	public SequenceSketch getStoredSequenceHash(SequenceId id)
	{
		return this.sequenceVectorsHash.get(id);
	}

	@Override
	public int size()
	{
		return this.sequenceVectorsHash.size();
	}
}
