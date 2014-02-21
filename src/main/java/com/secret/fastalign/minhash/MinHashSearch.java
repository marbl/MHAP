package com.secret.fastalign.minhash;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.general.AbstractHashSearch;
import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.MatchResult;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.general.SubSequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;

public final class MinHashSearch extends AbstractHashSearch<MinHash, SequenceMinHashes>
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
		
		
	}

	public static final int DEFAULT_SUB_KMER_SIZE = 12;
	public final static int INITIAL_HASH_SIZE = 10000;

	private final ArrayList<HashMap<Integer, ArrayList<SubSequenceId>>> hashes;
	private final ConcurrentHashMap<SequenceId, SequenceMinHashes> sequenceVectorsHash;
	private final boolean storeKmerInMemory;
	private final HashSet<Integer> filter;
	
	private final AtomicLong numberSubSequencesHit;
	private final AtomicLong numberSequencesHit;
	private final AtomicLong numberSequencesFullyCompared;

	private final int numMinMatches;
	private final int subSequenceSize;

	private static int[] errorString(int[] s, double readError, Random generator)
	{
		int[] snew = s.clone();

		for (int iter = 0; iter < s.length; iter++)
		{
			if (generator.nextDouble() < readError)
				while (snew[iter] == s[iter])
					snew[iter] = generator.nextInt(4);
		}

		return snew;
	}

	public static double probabilityKmerMatches(double readError, int kmerSize)
	{
		// allocate a random sequence of integers
		int[] s = new int[100000 + kmerSize];

		// allocate a random generator
		Random generator = new Random(1);

		// generate the random sequence
		for (int iter = 0; iter < s.length; iter++)
			s[iter] = generator.nextInt(4);

		// get the permuted matches
		int[] s1 = errorString(s, readError, generator);
		int[] s2 = errorString(s, readError, generator);

		// compute the expected number of matches
		int kmerCount = 0;
		int numKmers = s.length - kmerSize;
		for (int kmerIndex = 0; kmerIndex < numKmers; kmerIndex++)
		{
			boolean matches = true;
			for (int index = 0; index < kmerSize; index++)
				if (s1[kmerIndex + index] != s2[kmerIndex + index])
				{
					matches = false;
					break;
				}

			if (matches)
				kmerCount++;

		}

		double ratio = (double) kmerCount / (double) numKmers;

		return ratio;
	}

	public MinHashSearch(FastaData data, int kmerSize, int numHashes, int numMinMatches, int subSequenceSize, int numThreads,
			boolean storeKmerInMemory, boolean storeResults, HashSet<Integer> filter) throws IOException
	{
		super(kmerSize, numHashes, numThreads, storeResults);

		this.numMinMatches = numMinMatches;
		this.subSequenceSize = subSequenceSize;
		this.storeKmerInMemory = storeKmerInMemory;
		this.numberSubSequencesHit = new AtomicLong();
		this.numberSequencesHit = new AtomicLong();
		this.numberSequencesFullyCompared = new AtomicLong();
		
		//add the sequence filter
		this.filter = filter;

		// enqueue full file
		data.enqueueFullFile();

		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, SequenceMinHashes>(
				(int) data.getNumberProcessed() * 2 + 100, (float) 0.75, this.numThreads);

		this.hashes = new ArrayList<HashMap<Integer, ArrayList<SubSequenceId>>>(numHashes);
		for (int iter = 0; iter < numHashes; iter++)
			this.hashes.add(new HashMap<Integer, ArrayList<SubSequenceId>>((int) data.getNumberProcessed() * 2 * 10));

		addData(data);
	}

	@Override
	public boolean addDirectionalSequence(Sequence seq)
	{
		SequenceMinHashes currHash = getSequenceHash(seq);
		int[][] currMinHashes = currHash.getMainHashes().getSubSeqMinHashes();

		if (currMinHashes[0].length != this.hashes.size())
			throw new FastAlignRuntimeException("Number of hashes does not match.");

		// put the result into the hashmap
		SequenceMinHashes minHash = this.sequenceVectorsHash.put(seq.getId(), currHash);
		if (minHash != null)
		{
			// put back the original value, WARNING: not thread safe
			this.sequenceVectorsHash.put(seq.getId(), minHash);

			throw new FastAlignRuntimeException("Number of hashes does not match.");
		}

		// add the hashes
		int count = 0;
		for (HashMap<Integer, ArrayList<SubSequenceId>> hash : this.hashes)
		{
			for (int subSequences = 0; subSequences < currMinHashes.length; subSequences++)
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
					currList.add(new SubSequenceId(seq.getId(), (short)subSequences));
				}
			}

			count++;
		}

		return true;
	}

	@Override
	public List<MatchResult> findMatches(SequenceMinHashes seqMinHashes, double minScore, boolean allToAll)
	{
		MinHash minHash = seqMinHashes.getMainHashes();

		if (this.hashes.size() != minHash.numHashes())
			throw new FastAlignRuntimeException("Number of hashes does not match. Stored size " + this.hashes.size()
					+ ", input size " + minHash.numHashes() + ".");

		HashMap<MatchId, HitInfo> matchHitMap = new HashMap<MatchId, HitInfo>(size()/10+1);

		int[][] subSeqMinHashes = minHash.getSubSeqMinHashes();
		for (int subSequence = 0; subSequence < subSeqMinHashes.length; subSequence++)
		{
			int hashIndex = 0;
			for (HashMap<Integer,ArrayList<SubSequenceId>> currHash : this.hashes)
			{
				ArrayList<SubSequenceId> currentHashMatchList = currHash.get(subSeqMinHashes[subSequence][hashIndex]);

				// if some matches exist add them
				if (currentHashMatchList != null)
				{
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

		int[][] fullKmerMatch = null;
		for (Entry<SequenceId, Integer> match : bestSequenceHit.entrySet())
		{
			// do not store matches to yourself
			if (match.getKey().getHeaderId() == seqMinHashes.getSequenceId().getHeaderId())
				continue;
			// do not store matches smaller ids
			if (allToAll && match.getKey().getHeaderId() > seqMinHashes.getSequenceId().getHeaderId())
				continue;

			//see if the hit number is high enough
			if (match.getValue() >= this.numMinMatches)
			{
				// get the info for the id
				SequenceMinHashes matchedHash = this.sequenceVectorsHash.get(match.getKey());

				// if are not holding to the full kmer info load it
				if (fullKmerMatch == null)
					fullKmerMatch = seqMinHashes.getFullHashes();

				Pair<Double, Integer> result = seqMinHashes.getFullScore(fullKmerMatch, matchedHash);
				double matchScore = result.x;
				int shift = result.y;
				int shiftb = -shift - seqMinHashes.getSequenceLength() + matchedHash.getSequenceLength();
				
			  this.numberSequencesFullyCompared.getAndIncrement();

				if (matchScore >= minScore)
				{
					MatchResult currResult = new MatchResult(seqMinHashes.getSequenceId(), match.getKey(), matchScore, -shift, shiftb);

					// add to list
					matches.add(currResult);
				}
			}
		}

		return matches;
	}

	@Override
	public SequenceMinHashes getSequenceHash(Sequence seq)
	{
		return new SequenceMinHashes(seq, this.kmerSize, this.numWords, this.subSequenceSize, DEFAULT_SUB_KMER_SIZE,
				this.storeKmerInMemory, this.filter);
	}

	@Override
	public Collection<SequenceId> getStoredForwardSequenceIds()
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

	public AtomicLong getNumberSubSequencesHit()
	{
		return this.numberSubSequencesHit;
	}

	public long getNumberSequencesHit()
	{
		return this.numberSequencesHit.get();
	}

	public long getNumberSequencesFullyCompared()
	{
		return this.numberSequencesFullyCompared.get();
	}
}
