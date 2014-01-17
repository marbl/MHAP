package com.secret.fastalign.simhash;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import com.secret.fastalign.data.FastaData;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.main.MatchResult;

public abstract class AbstractVectorStore<T extends VectorHash<T>>
{
	protected final int kmerSize;
	protected final ArrayList<T> sequenceVectors;
	protected final int numWords;
	
	public AbstractVectorStore(int kmerSize, int numWords) 
	{
		this.kmerSize = kmerSize;
		this.sequenceVectors = new ArrayList<T>();
		this.numWords = numWords;
	}
	

	public void addData(FastaData data)
	{
		for (Sequence seq : data.getSequences())
			addSequence(seq);
	}
	
	public void addSequence(Sequence seq)
	{
		T simHash = getVectorHash(seq);
		
		this.sequenceVectors.add(simHash);
	}
	
	public abstract T getVectorHash(Sequence seq);
	
	public List<MatchResult> findMatches(Sequence seq, double acceptScore)
	{
		T simHash = getVectorHash(seq);
		
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (T hash : this.sequenceVectors)
		{
			if (seq.getId().equals(hash.getSequenceId()))
				continue;
			
			double score = simHash.correlation(hash);
			
			if (score>acceptScore)
				results.add(new MatchResult(seq.getId(), hash.getSequenceId(), score, 0));
		}
		
		Collections.sort(results);
		
		return results;
	}
}
