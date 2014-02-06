package com.secret.fastalign.simhash;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

import com.secret.fastalign.general.AbstractHashSearch;
import com.secret.fastalign.general.MatchResult;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.Pair;

public final class SimHashSearch extends AbstractHashSearch<SimHash,SequenceSimHashes>
{
	protected final ConcurrentHashMap<SequenceId, SequenceSimHashes> sequenceVectorsHash;
	
	protected static final int SUB_KMER_SIZE = 6;
	protected static final int SUB_STRING_SIZE = 200;	
	protected static final int SUB_WORD_SIZE = 2;
	
	public SimHashSearch(int kmerSize, int numWords) 
	{
		super(kmerSize, numWords);

		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, SequenceSimHashes>();
	}
	
	@Override
	public boolean addDirectionalSequence(Sequence seq)
	{
		//put the result into the hashmap
		SequenceSimHashes simHash = this.sequenceVectorsHash.put(seq.getId(), getSequenceHash(seq));
		
		if (simHash!=null)
		{
			//put back the original value, WARNING: not thread safe
			this.sequenceVectorsHash.put(seq.getId(), simHash);
		}
		else
			return simHash==null;
		
		//now add the subsequences
		//int len = seq.getString().length();
		//ArrayList<T> subHashes = new ArrayList<T>();
		//for (int iter=0; iter<len; iter++)
			//do nothing
		
		return true;
	}
	
	@Override
	public List<MatchResult> findMatches(SequenceSimHashes seqHash, double acceptScore, boolean allToAll)
	{		
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (SequenceSimHashes hash : this.sequenceVectorsHash.values())
		{
			if (seqHash.getSequenceId().getHeaderId()==hash.getSequenceId().getHeaderId())
				continue;
			if (allToAll && seqHash.getSequenceId().getHeaderId()>hash.getSequenceId().getHeaderId())
				continue;
			
			if (hash.getSequenceLength()<300 || seqHash.getSequenceLength()<300)
				continue;

			
			//compute the initial score
			double score = seqHash.jaccord(hash);
			//score = score*(double)(hash.sequenceLength()+seqHash.sequenceLength())/(double)Math.min(hash.sequenceLength(),seqHash.sequenceLength()*2.0);
			
			//acceptScore = 0.14;
			
			if (score>=acceptScore)
			{
				Pair<Double,Integer> result = seqHash.orderedScore(hash);
				//score = result.x;
				int shift = result.y*SUB_STRING_SIZE;
				int shiftb = 0;
				
				results.add(new MatchResult(seqHash.getSequenceId(), hash.getSequenceId(), score, -shift, shiftb));
			}
		}
		
		//not decided if should sort
		Collections.sort(results);
		
		return results;
	}


	@Override
	public SequenceSimHashes getSequenceHash(Sequence seq)
	{
		return new SequenceSimHashes(seq, this.kmerSize, this.numWords, SUB_STRING_SIZE, SUB_KMER_SIZE, SUB_WORD_SIZE);
	}


	@Override
	public Collection<SequenceId> getStoredForwardSequenceIds()
	{
		ArrayList<SequenceId> seqIds = new ArrayList<SequenceId>(this.sequenceVectorsHash.size());
		for (SequenceSimHashes hashes : this.sequenceVectorsHash.values())
			if (hashes.getSequenceId().isForward())
				seqIds.add(hashes.getSequenceId());
		
		return seqIds;
	}


	@Override
	public SequenceSimHashes getStoredSequenceHash(SequenceId id)
	{
		return this.sequenceVectorsHash.get(id);
	}


	@Override
	public int size()
	{
		return this.sequenceVectorsHash.size();
	}
}
