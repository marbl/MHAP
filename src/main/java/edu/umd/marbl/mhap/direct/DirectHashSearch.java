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
package edu.umd.marbl.mhap.direct;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import edu.umd.marbl.mhap.general.AbstractMatchSearch;
import edu.umd.marbl.mhap.general.AbstractSequenceHashStreamer;
import edu.umd.marbl.mhap.general.MatchResult;
import edu.umd.marbl.mhap.general.OrderKmerHashes;
import edu.umd.marbl.mhap.general.OverlapInfo;
import edu.umd.marbl.mhap.general.SequenceId;
import edu.umd.marbl.mhap.utils.FastAlignRuntimeException;
import edu.umd.marbl.mhap.utils.Utils;

public final class DirectHashSearch extends AbstractMatchSearch<SequenceDirectHashes>
{
	private final double acceptScore;

	private final HashMap<Integer, ArrayList<SequenceId>> hashes;
	private final double maxShift;
	private final int minStoreLength;
	private final AtomicLong numberElementsProcessed;
	private final AtomicLong numberSequencesFullyCompared;
	
	private final AtomicLong numberSequencesHashed;
	private final AtomicLong numberSequencesHit;

	private final int numHashes;
	private final int numMinMatches;

	private final ConcurrentHashMap<SequenceId, SequenceDirectHashes> sequenceVectorsHash;

	public DirectHashSearch(AbstractSequenceHashStreamer<SequenceDirectHashes> data, int numHashes, int numMinMatches, int numThreads, 
			boolean storeResults, int minStoreLength, double maxShift, double acceptScore) throws IOException
	{
		super(numThreads, storeResults);
		
		this.minStoreLength = minStoreLength;
		this.numMinMatches = numMinMatches;
		this.numHashes = numHashes;
		this.maxShift = maxShift;
		this.acceptScore = acceptScore;
		this.numberSequencesHit = new AtomicLong();
		this.numberSequencesFullyCompared = new AtomicLong();
		this.numberSequencesHashed = new AtomicLong();
		this.numberElementsProcessed = new AtomicLong();
		
		// enqueue full file, since have to know full size
		data.enqueueFullFile(false, this.numThreads);

		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, SequenceDirectHashes>(data.getNumberProcessed() + 100, (float) 0.75, this.numThreads);

		this.hashes = new HashMap<Integer, ArrayList<SequenceId>>(data.getNumberProcessed()*5000+100);
		
		addData(data);
	}

	@Override
	public boolean addSequence(SequenceDirectHashes currHash)
	{
		OrderKmerHashes currOrderedHashes = currHash.getMainHashes();
		SequenceId seqId = currHash.getSequenceId();

		// put the result into the hashmap
		SequenceDirectHashes minHash = this.sequenceVectorsHash.put(seqId, currHash);
		if (minHash != null)
		{
			// put back the original value, WARNING: not thread safe
			this.sequenceVectorsHash.put(seqId, minHash);

			throw new FastAlignRuntimeException("Sequence id already exists in the hashtable.");
		}
		
		// add the hashes
		ArrayList<SequenceId> currList;
		int prevValue = Integer.MAX_VALUE;
		for (int kmer = 0; kmer < currOrderedHashes.size(); kmer++)
		{
			final int hashVal = currOrderedHashes.getHash(kmer);
			
			//don't hash repeats
			if (hashVal==prevValue)
				continue;

			// get the list
			synchronized (this.hashes)
			{
				currList = this.hashes.get(hashVal);

				if (currList == null)
				{
					currList = new ArrayList<SequenceId>(2);
					this.hashes.put(hashVal, currList);
				}
			}

			// add the element
			synchronized (currList)
			{
				currList.add(seqId);
			}

			prevValue = hashVal;
		}
		
		//increment the counter
		this.numberSequencesHashed.getAndIncrement();
		
		return true;
	}

	@Override
	public List<MatchResult> findMatches(SequenceDirectHashes seqHashes, boolean toSelf)
	{
		//create a random generator
		Random rand = new Random(0);
		
		OrderKmerHashes mainHashes = seqHashes.getMainHashes();
		int numStoredHashes = mainHashes.size();
		
		HashMap<SequenceId, Integer> bestSequenceHit = new HashMap<SequenceId, Integer>(1024);
		int lookupSize = this.numHashes>0 ? this.numHashes : numStoredHashes;
		for (int kmer = 0; kmer < lookupSize; kmer++)
		{
			int randIndex = kmer;
			if (this.numHashes>0)
				randIndex = rand.nextInt(numStoredHashes);
			
			//get the random hash
			int hashValue = mainHashes.getHash(randIndex);
			ArrayList<SequenceId> currentHashMatchList = this.hashes.get(hashValue);

			// if some matches exist add them
			if (currentHashMatchList != null)
			{
				this.numberElementsProcessed.getAndAdd(currentHashMatchList.size());
				for (SequenceId id : currentHashMatchList)
				{
					Integer count = bestSequenceHit.get(id);
					if (count==null)
						count = new Integer(0);
					
					bestSequenceHit.put(id, count.intValue()+1);
				}
			}
		}
		
		this.numberSequencesHit.getAndAdd(bestSequenceHit.size());
		
		// compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>(32);

		for (Entry<SequenceId, Integer> match : bestSequenceHit.entrySet())
		{
			//get the match id
			SequenceId matchId = match.getKey();
			
			// do not store matches with smaller ids, unless its coming from a short read
			if (toSelf && matchId.getHeaderId() == seqHashes.getSequenceId().getHeaderId())
				continue;

			//see if the hit number is high enough			
			if (match.getValue() >= this.numMinMatches)
			{
				SequenceDirectHashes matchedHashes = this.sequenceVectorsHash.get(match.getKey());
				if (matchedHashes==null)
					throw new FastAlignRuntimeException("Hashes not found for given id.");
				
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
				
				OverlapInfo result = seqHashes.getOrderedHashes().getFullScore(matchedHashes.getOrderedHashes(), this.maxShift);
				
				//increment the counter
				this.numberSequencesFullyCompared.getAndIncrement();
				
				//if score is good add
				if (result.score >= this.acceptScore)
				{
					MatchResult currResult = new MatchResult(seqHashes.getSequenceId(), matchId, result.score, result.a, result.b);

					// add to list
					matches.add(currResult);
				}
			}
		}

		return matches;
	}

	public long getNumberElementsProcessed()
	{
		return this.numberElementsProcessed.get();
	}

	public long getNumberSequenceHashed()
	{
		return this.numberSequencesHashed.get();
	}

	public long getNumberSequencesFullyCompared()
	{
		return this.numberSequencesFullyCompared.get();
	}
	
	public long getNumberSequencesHit()
	{
		return this.numberSequencesHit.get();
	}
	
	@Override
	public List<SequenceId> getStoredForwardSequenceIds()
	{
		ArrayList<SequenceId> seqIds = new ArrayList<SequenceId>(this.sequenceVectorsHash.size());
		for (SequenceDirectHashes hashes : this.sequenceVectorsHash.values())
			if (hashes.getSequenceId().isForward())
				seqIds.add(hashes.getSequenceId());
		
		return seqIds;
	}
	
	@Override
	public SequenceDirectHashes getStoredSequenceHash(SequenceId id)
	{
		return this.sequenceVectorsHash.get(id);
	}

	public double hashTableNormalizedEnthropy()
	{
		return Utils.hashEfficiency(this.hashes);
	}

	@Override
	public int size()
	{
		return this.sequenceVectorsHash.size();
	}
}
