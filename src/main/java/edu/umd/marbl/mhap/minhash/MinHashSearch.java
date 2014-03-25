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
package edu.umd.marbl.mhap.minhash;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicLong;

import edu.umd.marbl.mhap.general.AbstractMatchSearch;
import edu.umd.marbl.mhap.general.AbstractSequenceHashStreamer;
import edu.umd.marbl.mhap.general.MatchResult;
import edu.umd.marbl.mhap.general.SequenceId;
import edu.umd.marbl.mhap.general.SubSequenceId;
import edu.umd.marbl.mhap.utils.FastAlignRuntimeException;
import edu.umd.marbl.mhap.utils.Pair;
import edu.umd.marbl.mhap.utils.Utils;

public final class MinHashSearch extends AbstractMatchSearch<SequenceMinHashes>
{
	public static final class HitInfo
	{
		public int count;

		public HitInfo()
		{
			this.count = 0;
		}

		public void addHit()
		{
			this.count++;
		}
	}
	
	public static final class MatchId
	{
		public final SubSequenceId id;
		public final int toSubSequenceNum;
		
		public MatchId(SubSequenceId id, int toSubSequenceNum)
		{
			this.id = id;
			this.toSubSequenceNum = toSubSequenceNum;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj)
		{
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			MatchId other = (MatchId) obj;
			if (this.toSubSequenceNum != other.toSubSequenceNum)
				return false;
			if (this.id == null)
			{
				if (other.id != null)
					return false;
			}
			else if (!this.id.equals(other.id))
				return false;
			return true;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode()
		{
			final int prime = 31;
			int result = 1;
			result = prime * result + ((this.id == null) ? 0 : this.id.hashCode());
			result = prime * result + this.toSubSequenceNum;
			return result;
		}
		
		
	}

	private final double acceptScore;

	private final ArrayList<HashMap<Integer, ArrayList<SubSequenceId>>> hashes;
	private final double maxShift;
	private final int minStoreLength;
	private final AtomicLong numberElementsProcessed;
	private final AtomicLong numberSequencesFullyCompared;

	private final AtomicLong numberSequencesHit;
	private final AtomicLong numberSequencesMinHashed;
	private final AtomicLong numberSubSequences;
	private final AtomicLong numberSubSequencesHit;

	private final int numMinMatches;
	private final HashMap<SequenceId, SequenceMinHashes> sequenceVectorsHash;

	
	public MinHashSearch(AbstractSequenceHashStreamer<SequenceMinHashes> data, int numHashes, int numMinMatches, int numThreads, 
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
		
		// enqueue full file, since have to know full size
		data.enqueueFullFile(false, this.numThreads);

		this.sequenceVectorsHash = new HashMap<SequenceId, SequenceMinHashes>(data.getNumberProcessed() + 100, (float) 0.75);

		this.hashes = new ArrayList<HashMap<Integer, ArrayList<SubSequenceId>>>(numHashes);
		for (int iter = 0; iter < numHashes; iter++)
			this.hashes.add(new HashMap<Integer, ArrayList<SubSequenceId>>(data.getNumberSubSequencesProcessed()+100));
		
		addData(data);
	}

	@Override
	public boolean addSequence(SequenceMinHashes currHash)
	{
		int[][] currMinHashes = currHash.getMinHashes().getSubSeqMinHashArray();

		if (currMinHashes[0].length != this.hashes.size())
			throw new FastAlignRuntimeException("Number of hashes does not match.");

		// put the result into the hashmap
		synchronized (this.sequenceVectorsHash)
		{
			SequenceMinHashes minHash = this.sequenceVectorsHash.put(currHash.getSequenceId(), currHash);
			if (minHash != null)
			{
				this.sequenceVectorsHash.put(currHash.getSequenceId(), minHash);

				throw new FastAlignRuntimeException("Sequence id already exists in the hashtable.");
			}			
		}
		
		// add the hashes
		for (int subSequences = 0; subSequences < currMinHashes.length; subSequences++)
		{
			SubSequenceId subId = new SubSequenceId(currHash.getSequenceId(), (short)subSequences);
			
			int count = 0;
			for (HashMap<Integer, ArrayList<SubSequenceId>> hash : this.hashes)
			{
				ArrayList<SubSequenceId> currList;
				final int hashVal = currMinHashes[subSequences][count];

				// get the list
				synchronized (hash)
				{
					currList = hash.get(hashVal);

					if (currList == null)
					{
						currList = new ArrayList<SubSequenceId>(2);
						hash.put(hashVal, currList);
					}
				}

				// add the element
				synchronized (currList)
				{
					currList.add(subId);
				}
				
				count++;
			}

			//increment the subsequence counter 
			this.numberSubSequences.getAndIncrement();
		}
		
		//increment the counter
		this.numberSequencesMinHashed.getAndIncrement();
		
		return true;
	}

	@Override
	public List<MatchResult> findMatches(SequenceMinHashes seqMinHashes, boolean toSelf)
	{
		MinHash minHash = seqMinHashes.getMinHashes();

		if (this.hashes.size() != minHash.numHashes())
			throw new FastAlignRuntimeException("Number of hashes does not match. Stored size " + this.hashes.size()
					+ ", input size " + minHash.numHashes() + ".");

		HashMap<MatchId, HitInfo> matchHitMap = new HashMap<MatchId, HitInfo>(this.numberSubSequences.intValue()/5+1);

		int[][] subSeqMinHashes = minHash.getSubSeqMinHashArray();
		for (int subSequence = 0; subSequence < subSeqMinHashes.length; subSequence++)
		{
			int hashIndex = 0;
			for (HashMap<Integer,ArrayList<SubSequenceId>> currHash : this.hashes)
			{
				ArrayList<SubSequenceId> currentHashMatchList = currHash.get(subSeqMinHashes[subSequence][hashIndex]);

				// if some matches exist add them
				if (currentHashMatchList != null)
				{
					this.numberElementsProcessed.getAndAdd(currentHashMatchList.size());

					for (SubSequenceId subSequenceId : currentHashMatchList)
					{
						MatchId matchedId = new MatchId(subSequenceId, subSequence);
						
						// get current count in the list
						HitInfo currentHitInfo = matchHitMap.get(matchedId);

						// increment the count
						if (currentHitInfo == null)
						{
							currentHitInfo = new HitInfo();
							matchHitMap.put(matchedId, currentHitInfo);
						}

						// record the match of the kmer hash
						currentHitInfo.addHit();
					}
				}
				
				hashIndex++;
			}
		}
		
		//record number of hash matches processed
		this.numberSubSequencesHit.getAndAdd(matchHitMap.size());
		
		//shrink the list to only the best hits per sequence
		HashMap<SequenceId, Integer> bestSequenceHit = new HashMap<SequenceId, Integer>(matchHitMap.size());
		for (Entry<MatchId, HitInfo> match : matchHitMap.entrySet())
		{
			//get the current number
			SequenceId currId = match.getKey().id.getId();
			int currValue = match.getValue().count;

			//get the current best count
			Integer prevBestCount = bestSequenceHit.get(currId);
			
			if (prevBestCount==null || prevBestCount<currValue)
				bestSequenceHit.put(currId, currValue);
		}
		
		this.numberSequencesHit.getAndAdd(bestSequenceHit.size());
		
		// compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>(32);
		
		for (Entry<SequenceId, Integer> match : bestSequenceHit.entrySet())
		{
			//get the match id
			SequenceId matchId = match.getKey();
			
			// do not store matches with smaller ids, unless its coming from a short read
			if (toSelf && matchId.getHeaderId() == seqMinHashes.getSequenceId().getHeaderId())
				continue;

			//see if the hit number is high enough			
			if (match.getValue() >= this.numMinMatches)
			{
				SequenceMinHashes matchedHashes = this.sequenceVectorsHash.get(match.getKey());
				if (matchedHashes==null)
					throw new FastAlignRuntimeException("Hashes not found for given id.");
				
				//never process short to short
				if (matchedHashes.getSequenceLength()<this.minStoreLength && seqMinHashes.getSequenceLength()<this.minStoreLength)
					continue;
				
				//never process long to long in self, with greater id
				if (toSelf 
						&& matchId.getHeaderId() > seqMinHashes.getSequenceId().getHeaderId()
						&& matchedHashes.getSequenceLength()>=this.minStoreLength
						&& seqMinHashes.getSequenceLength()>=this.minStoreLength)
					continue;
				
				//never do short to long
				if (toSelf 
						&& matchedHashes.getSequenceLength()<this.minStoreLength
						&& seqMinHashes.getSequenceLength()>=this.minStoreLength)
					continue;
				
				Pair<Double, Integer> result = seqMinHashes.getOrderedHashes().getFullScore(matchedHashes.getOrderedHashes(), this.maxShift);
				double matchScore = result.x;
				int shift = result.y;
				int shiftb = -shift - seqMinHashes.getSequenceLength() + matchedHashes.getSequenceLength();
				
				this.numberSequencesFullyCompared.getAndIncrement();
				
				if (matchScore >= this.acceptScore)
				{
					MatchResult currResult = new MatchResult(seqMinHashes.getSequenceId(), matchId, matchScore, -shift, shiftb);

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
	
	public double hashTableNormalizedEnthropy()
	{
		double sum = 0.0;
		for (HashMap<Integer, ArrayList<SubSequenceId>> map : this.hashes)
		{
			sum+=Utils.hashEfficiency(map);
		}
		
		return sum/(double)this.hashes.size();
	}
	
	@Override
	public List<SequenceId> getStoredForwardSequenceIds()
	{
		ArrayList<SequenceId> seqIds = new ArrayList<SequenceId>(this.sequenceVectorsHash.size());
		for (SequenceMinHashes hashes : this.sequenceVectorsHash.values())
			if (hashes.getSequenceId().isForward())
				seqIds.add(hashes.getSequenceId());
		
		return seqIds;
	}

	@Override
	public SequenceMinHashes getStoredSequenceHash(SequenceId id)
	{
		return this.sequenceVectorsHash.get(id);
	}

	@Override
	public int size()
	{
		return this.sequenceVectorsHash.size();
	}
}