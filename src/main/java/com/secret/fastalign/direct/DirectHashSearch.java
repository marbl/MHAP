package com.secret.fastalign.direct;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.general.AbstractMatchSearch;
import com.secret.fastalign.general.AbstractSequenceHashStreamer;
import com.secret.fastalign.general.MatchResult;
import com.secret.fastalign.general.OrderKmerHashes;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;

public final class DirectHashSearch extends AbstractMatchSearch<SequenceDirectHashes>
{
	private final HashMap<Integer, ArrayList<SequenceId>> hashes;

	private final int minStoreLength;
	private final AtomicLong numberSequencesFullyCompared;
	private final AtomicLong numberSequencesHit;
	private final AtomicLong numberSequencesHashed;
	private final AtomicLong numberElementsProcessed;
	
	private final int numMinMatches;
	private final ConcurrentHashMap<SequenceId, SequenceDirectHashes> sequenceVectorsHash;

	private final int maxShift;
	private final double acceptScore;

	private final int numHashes;

	public DirectHashSearch(AbstractSequenceHashStreamer<SequenceDirectHashes> data, int numHashes, int numMinMatches, int numThreads, 
			boolean storeResults, int minStoreLength, int maxShift, double acceptScore) throws IOException
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
				
				Pair<Double, Integer> result = seqHashes.getOrderedHashes().getFullScore(matchedHashes.getOrderedHashes(), this.maxShift);
				double matchScore = result.x;
				int shift = result.y;
				int shiftb = -shift - seqHashes.getSequenceLength() + matchedHashes.getSequenceLength();
				
				this.numberSequencesFullyCompared.getAndIncrement();
				
				if (matchScore >= this.acceptScore)
				{
					MatchResult currResult = new MatchResult(seqHashes.getSequenceId(), matchId, matchScore, -shift, shiftb);

					// add to list
					matches.add(currResult);
				}
			}
		}

		return matches;
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
	
	public long getNumberElementsProcessed()
	{
		return this.numberElementsProcessed.get();
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

	@Override
	public int size()
	{
		return this.sequenceVectorsHash.size();
	}
}
