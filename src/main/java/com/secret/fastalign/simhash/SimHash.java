package com.secret.fastalign.simhash;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.secret.fastalign.data.FastaData;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.main.MatchResult;

public final class SimHash
{
	private final int kmerSize;
	private final ArrayList<SequenceSimHash> hashes;
	private final int numWords;

	public SimHash(int kmerSize, int numWords)
	{
		this.kmerSize = kmerSize;
		this.hashes = new ArrayList<SequenceSimHash>();
		this.numWords = numWords;
	}
	
	public void addData(FastaData data)
	{
		for (Sequence seq : data.getSequences())
			addSequence(seq);
	}
	
	public void addSequence(Sequence seq)
	{
		SequenceSimHash hash = new SequenceSimHash(seq, this.kmerSize, this.numWords);
		
		this.hashes.add(hash);
	}
	
	public List<MatchResult> findMatches(Sequence seq, double acceptScore)
	{
		SequenceSimHash simHash = new SequenceSimHash(seq, this.kmerSize, this.numWords);
		
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (SequenceSimHash hash : this.hashes)
		{
			if (seq.getId().equals(hash.getSequenceId()))
				continue;
			
			double score = simHash.score(hash);
			
			if (score>acceptScore)
				results.add(new MatchResult(seq.getId(), hash.getSequenceId(), score, 0));
		}
		
		Collections.sort(results);
		
		return results;
	}
}
