package com.secret.fastalign.hash;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public final class MinHash
{
	private final ArrayList<HashMap<Long, ArrayList<SequenceId>>> hashes;
	private final HashMap<SequenceId, Integer> seqLengths;
	private final int numHashes;
	
	private static long hashBasic(String kmer, long salt) 
  {
	  long key = 5381 + kmer.hashCode();
		   
	  key = key * 37;
	  key = key + salt;
	  key ^= key >> 33;
	  key *= 0xff51afd7ed558ccdL;
	  key ^= key >> 33;
	  key *= 0xc4ceb9fe1a85ec53L;
	  key ^= key >> 33;
		 
	  return(key);
 }
	
	public MinHash(int numHashes)
	{
		this.numHashes = numHashes;
	  this.hashes = new ArrayList<HashMap<Long, ArrayList<SequenceId>>>();
	  this.seqLengths = new HashMap<SequenceId, Integer>();
	  
	  //allocate the hashtables
	  for (int iter=0; iter<this.numHashes; iter++)
	  	this.hashes.add(new HashMap<Long, ArrayList<SequenceId>>());
	}
	
	public void addSequence(Sequence seq)
	{
		//record the length
		this.seqLengths.put(seq.getId(), seq.length());
		
		long[] hashValues = getHashes(seq);
		
		int count = 0;
		for (HashMap<Long, ArrayList<SequenceId>> hash : this.hashes)
		{
			ArrayList<SequenceId> currList = hash.get(hashValues[count]);
			
			if (currList==null)
				currList = new ArrayList<SequenceId>();
			
			currList.add(seq.getId());
			
			count++;
		}
	}
	
	public List<MatchResult> findMatches(Sequence seq, double threshold)
	{
		HashMap<SequenceId, Integer> counts = new HashMap<SequenceId, Integer>();
		
		//go through each hash and collect hits
		long[] hashes = getHashes(seq);
			
		for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
		{
			ArrayList<SequenceId> matchedList = this.hashes.get(hashIndex).get(hashes[hashIndex]);
			
			//if some matches exist add them
			if (matchedList!=null)
			{
				for (SequenceId matchedId : matchedList)
				{
					//get current count in the list
					Integer currentCount = counts.get(matchedId);
					
					//increment the count
					if (currentCount==null)
						counts.put(matchedId, 1);
					else
						counts.put(matchedId, currentCount+1);
				}
			}
		}
				
		//compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>();
		
		for (SequenceId id : counts.keySet())
		{
			//get the count
			int count = counts.get(id);
			
			//get the sequence length ratios
			double ratio = (double)this.seqLengths.get(id)/(double)seq.length();
			if (ratio<1.0)
				ratio = 1.0/ratio;
			
			double score = (double)count*ratio;
			
			if (score>threshold)
				matches.add(new MatchResult(seq.getId(), id, score));
		}
		
		//sort the result
		Collections.sort(matches);
		
		return matches;
	}
	
	public void addData(ReadData data)
	{
		for (Sequence seq : data.getSequences())
		{
			addSequence(seq);
		}
	}
	
	private long[] getHashes(Sequence seq)
	{
		long[] hashes = new long[this.numHashes];
		
		for (int iter=0; iter<seq.numKmers(); iter++)
		{
			String kmer = seq.getKmer(iter);
			
			for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
			{
				long salt = hashIndex*10;
				
				//compute the hash
				long hash = hashBasic(kmer, salt);
				
				//store minimum hash
				if (hash < hashes[hashIndex])
					hashes[hashIndex] = hash;
			}
		}
		
		return hashes;
	}
}
