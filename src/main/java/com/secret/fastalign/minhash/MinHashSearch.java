package com.secret.fastalign.minhash;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import com.secret.fastalign.general.AbstractHashSearch;
import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.MatchResult;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;

public final class MinHashSearch extends AbstractHashSearch<MinHash, SequenceMinHashes>
{
	public static final class HitInfo
	{
		private int count;

		public HitInfo()
		{
			this.count = 0;
		}

		public void addHit()
		{
			this.count++;
		}
	}

	protected static final int DEFAULT_SUB_KMER_SIZE = 12;
	public final static int INITIAL_HASH_SIZE = 10000;

	private final ArrayList<HashMap<Integer, ArrayList<SequenceId>>> hashes;
	private final ConcurrentHashMap<SequenceId, SequenceMinHashes> sequenceVectorsHash;
	private final boolean storeKmerInMemory;

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
			boolean storeKmerInMemory, boolean storeResults) throws IOException
	{
		super(kmerSize, numHashes, numThreads, storeResults);

		this.numMinMatches = numMinMatches;
		this.subSequenceSize = subSequenceSize;
		this.storeKmerInMemory = storeKmerInMemory;

		// enqueue full file
		data.enqueueFullFile();

		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, SequenceMinHashes>(
				(int) data.getNumberProcessed() * 2 + 100, (float) 0.75, this.numThreads);

		this.hashes = new ArrayList<HashMap<Integer, ArrayList<SequenceId>>>(numHashes);
		for (int iter = 0; iter < numHashes; iter++)
			this.hashes.add(new HashMap<Integer, ArrayList<SequenceId>>((int) data.getNumberProcessed() * 2 * 20));

		addData(data);
	}

	@Override
	public boolean addDirectionalSequence(Sequence seq)
	{
		SequenceMinHashes currHash = getSequenceHash(seq);

		if (currHash.getMainHashes().numHashes() != this.hashes.size())
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
		for (HashMap<Integer, ArrayList<SequenceId>> hash : this.hashes)
		{
			int[][] minHashes = currHash.getMainHashes().getSubSeqMinHashes();

			for (int subSequences = 0; subSequences < minHashes.length; subSequences++)
			{
				ArrayList<SequenceId> currList;

				// get the list
				synchronized (hash)
				{
					final int hashVal = minHashes[subSequences][count];
					currList = hash.get(hashVal);

					if (currList == null)
					{
						currList = new ArrayList<SequenceId>(16);
						hash.put(hashVal, currList);
					}
				}

				// add the element
				synchronized (currList)
				{
					currList.add(seq.getId());
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

		HashMap<SequenceId, HitInfo> matchHitMap = new HashMap<SequenceId, HitInfo>(size()/10+1);

		int[][] subSeqMinHashes = minHash.getSubSeqMinHashes();
		for (int subSequences = 0; subSequences < subSeqMinHashes.length; subSequences++)
		{
			int hashIndex = 0;
			for (HashMap<Integer,ArrayList<SequenceId>> currHash : this.hashes)
			{
				ArrayList<SequenceId> currentHashMatchList = currHash.get(subSeqMinHashes[subSequences][hashIndex]);

				// if some matches exist add them
				if (currentHashMatchList != null)
				{
					for (SequenceId matchedId : currentHashMatchList)
					{
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

		// compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>(32);

		int[][] fullKmerMatch = null;
		for (SequenceId id : matchHitMap.keySet())
		{
			// do not store matches to yourself
			if (id.getHeaderId() == seqMinHashes.getSequenceId().getHeaderId())
				continue;
			// do not store matches smaller ids
			if (allToAll && id.getHeaderId() > seqMinHashes.getSequenceId().getHeaderId())
				continue;

			// get the hit info
			HitInfo hit = matchHitMap.get(id);

			if (hit.count >= this.numMinMatches)
			{
				// get the info for the id
				SequenceMinHashes matchedHash = this.sequenceVectorsHash.get(id);

				// if are not holding to the full kmer info load it
				if (fullKmerMatch != null)
					fullKmerMatch = seqMinHashes.getFullHashes();

				Pair<Double, Integer> result = seqMinHashes.getFullScore(fullKmerMatch, matchedHash);
				double matchScore = result.x;
				int shift = result.y;
				int shiftb = -shift - seqMinHashes.getSequenceLength() + matchedHash.getSequenceLength();

				if (matchScore >= minScore)
				{
					MatchResult currResult = new MatchResult(seqMinHashes.getSequenceId(), id, matchScore, -shift, shiftb);

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
				this.storeKmerInMemory);
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
}
