package com.secret.fastalign.minhash;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import com.secret.fastalign.general.AbstractHashSearch;
import com.secret.fastalign.general.MatchResult;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;

public final class MinHashSearch extends AbstractHashSearch<MinHash,SequenceMinHashes>
{
	public final class HitInfo
	{
		public int count;
		public int posSumFrom;
		public int posSumTo;
		public int uniqueHits;
		
		public HitInfo()
		{
			this.count = 1;
			this.uniqueHits = 0;
			this.posSumFrom = 0;
			this.posSumTo = 0;
		}
		
		public HitInfo(int posFrom, int posTo)
		{
			this.count = 1;
			this.uniqueHits = 1;
			this.posSumFrom = posFrom;
			this.posSumTo = posTo;
		}

		public void addNonUniqueHit()
		{
			this.count++;
		}
		
		public void addUniqueHit(int posFrom, int posTo)
		{
			this.count++;
			this.posSumFrom+=posFrom;
			this.posSumTo+=posTo;
			this.uniqueHits++;
		}
		
		public int fromShift()
		{
			double averageFrom = (double)this.posSumFrom/(double)this.uniqueHits;
			double averageTo = (double)this.posSumTo/(double)this.uniqueHits;
			
			return (int)Math.round(averageTo-averageFrom);
		}

	}
	
	protected static final int SUB_KMER_SIZE = 16;
	//protected static final int SUB_STRING_SIZE = 25;	
	//protected static final int SUB_WORD_SIZE = 16;
	
	private final ArrayList<HashMap<Integer, ArrayList<KmerInfo>>> hashes;
	//private final HashMap<SequenceId, Integer> seqLengths;
  private final ConcurrentHashMap<SequenceId, SequenceMinHashes> sequenceVectorsHash;
	
	private static int[] errorString(int[] s, double readError, Random generator)
	{
		int[] snew = s.clone();
		
		for (int iter=0; iter<s.length; iter++)
		{
			if (generator.nextDouble()<readError)
				while(snew[iter]==s[iter])
					snew[iter] = generator.nextInt(4);
		}
		
		return snew;
	}
	
	public static double probabilityKmerMatches(double readError, int kmerSize)
	{
		//allocate a random sequence of integers
		int[] s = new int[100000+kmerSize];

		//allocate a random generator
		Random generator = new Random(1);
		
		//generate the random sequence
		for (int iter=0; iter<s.length; iter++)
			s[iter] = generator.nextInt(4);
		
		//get the permuted matches
		int[] s1 = errorString(s, readError, generator);
		int[] s2 = errorString(s, readError, generator);
		
		//compute the expected number of matches
		int kmerCount = 0;
		int numKmers = s.length-kmerSize;
		for (int kmerIndex=0; kmerIndex<numKmers; kmerIndex++)
		{
			boolean matches = true;
			for (int index=0; index<kmerSize; index++)
				if (s1[kmerIndex+index]!=s2[kmerIndex+index])
				{
					matches = false;
					break;
				}
			
			if (matches)
				kmerCount++;
			
		}
		
		double ratio = (double)kmerCount/(double)numKmers;
						
		return ratio;
	}
	
	public MinHashSearch(int kmerSize, int numHashes)
	{
		super(kmerSize, numHashes);
		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, SequenceMinHashes>();
		
		this.hashes = new ArrayList<HashMap<Integer, ArrayList<KmerInfo>>>(numHashes);
		for (int iter=0; iter<numHashes; iter++)
			this.hashes.add(new HashMap<Integer, ArrayList<KmerInfo>>());
	}
	
	@Override
	public boolean addDirectionalSequence(Sequence seq)
	{
		SequenceMinHashes currHash = getSequenceHash(seq);
		
		//put the result into the hashmap
		SequenceMinHashes minHash = this.sequenceVectorsHash.put(seq.getId(), currHash);		
		if (minHash!=null)
		{
			//put back the original value, WARNING: not thread safe
			this.sequenceVectorsHash.put(seq.getId(), minHash);
			
			return false;
		}
		
		if (currHash.getMainHashes().numHashes()!=this.hashes.size())
			throw new FastAlignRuntimeException("Number of hashes does not match.");
		
		//add the hashes
		int count = 0;
		for (HashMap<Integer, ArrayList<KmerInfo>> hash : this.hashes)
		{
			synchronized (hash)
			{
				ArrayList<KmerInfo> currList = hash.get(currHash.getMainHashes().minHashes[count]);
				
				if (currList==null)
				{
					currList = new ArrayList<KmerInfo>();
					hash.put(currHash.getMainHashes().minHashes[count], currList);
				}
				
				currList.add(new KmerInfo(seq.getId(), currHash.getMainHashes().minHashes[count]));			
			}
			
			count++;
		}
		
		return true;
	}
	
	@Override
	public List<MatchResult> findMatches(SequenceMinHashes seqMinHashes, double score, boolean allToAll)
	{
		MinHash minHash = seqMinHashes.getMainHashes();
		
		if (this.hashes.size()!=minHash.numHashes())
			throw new FastAlignRuntimeException("Number of hashes does not match. Stored size "+this.hashes.size()+", input size "+minHash.numHashes()+".");
		
		HashMap<SequenceId, HitInfo> matchHitMap = new HashMap<SequenceId, HitInfo>();
			
		for (int hashIndex=0; hashIndex<minHash.numHashes(); hashIndex++)
		{
			ArrayList<KmerInfo> currentHashMatchList = this.hashes.get(hashIndex).get(minHash.minHashes[hashIndex]);
			
			//if some matches exist add them
			if (currentHashMatchList!=null)
			{
				for (KmerInfo matchedId : currentHashMatchList)
				{
					//get current count in the list
					HitInfo currentHitInfo = matchHitMap.get(matchedId.getId());
					
					//increment the count
					if (currentHitInfo==null)
					{
						//if this kmer is unique then use it for positioning
						if (minHash.getHashPositions()[hashIndex]>=0 && matchedId.getPosition()>=0)
							matchHitMap.put(matchedId.getId(), new HitInfo(minHash.getHashPositions()[hashIndex], matchedId.getPosition()));
						else
							matchHitMap.put(matchedId.getId(), new HitInfo());
					}
					else
					{
						//if this kmer is unique then use it for positioning
						if (minHash.getHashPositions()[hashIndex]>=0 && matchedId.getPosition()>=0)
							currentHitInfo.addUniqueHit(minHash.getHashPositions()[hashIndex], matchedId.getPosition());
						else
							currentHitInfo.addNonUniqueHit();
					}
				}
			}
		}
				
		//compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>();
		
		for (SequenceId id : matchHitMap.keySet())
		{
			//do not store matches to yourself
			if (id.getHeaderId()==seqMinHashes.getSequenceId().getHeaderId())
				continue;
			//do not store matches smaller ids
			if (allToAll && id.getHeaderId()<seqMinHashes.getSequenceId().getHeaderId())
				continue;
			
			//get the hit info
			HitInfo hit = matchHitMap.get(id);
			
			double matchScore = (double)hit.count/(double)this.numWords;
						
			if (hit.count>=2)
			{
				//get the info for the id
				SequenceMinHashes matchedHash = this.sequenceVectorsHash.get(id);
				
				Pair<Double,Integer> result = seqMinHashes.orderedScore(matchedHash);
				matchScore = result.x;
				int shift = result.y;
				
				//if (matchScore>0.25)
					matches.add(new MatchResult(minHash.getSequenceId(), id, matchScore, shift));
			}
		}
		
		return matches;
	}

	@Override
	public SequenceMinHashes getSequenceHash(Sequence seq)
	{
		return new SequenceMinHashes(seq, this.kmerSize, this.numWords, SUB_KMER_SIZE);
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
