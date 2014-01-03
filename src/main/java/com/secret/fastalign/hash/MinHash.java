package com.secret.fastalign.hash;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import com.secret.fastalign.utils.CyclicHash;

public final class MinHash
{
	private final class HitInfo
	{
		public int count;
		public int uniqueHits;
		public int posSumFrom;
		public int posSumTo;
		
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

		public void addUniqueHit(int posFrom, int posTo)
		{
			this.count++;
			this.posSumFrom+=posFrom;
			this.posSumTo+=posTo;
			this.uniqueHits++;
		}
		
		public void addNonUniqueHit()
		{
			this.count++;
		}
		
		public int fromShift()
		{
			double averageFrom = (double)this.posSumFrom/(double)this.uniqueHits;
			double averageTo = (double)this.posSumTo/(double)this.uniqueHits;
			
			return (int)Math.round(averageTo-averageFrom);
		}

	}
	
	private final ArrayList<HashMap<Integer, ArrayList<KmerInfo>>> hashes;
	private final HashMap<SequenceId, Integer> seqLengths;
	private final int numHashes;
	private final int kmerSize;
	
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
	
	public MinHash(int numHashes, int kmerSize)
	{
		this.kmerSize = kmerSize;
		this.numHashes = numHashes;
	  this.hashes = new ArrayList<HashMap<Integer, ArrayList<KmerInfo>>>();
	  this.seqLengths = new HashMap<SequenceId, Integer>();
	  
	  //allocate the hashtables
	  for (int iter=0; iter<this.numHashes; iter++)
	  	this.hashes.add(new HashMap<Integer, ArrayList<KmerInfo>>());
	}
	
	public void addSequence(Sequence seq)
	{
		//record the length
		this.seqLengths.put(seq.getId(), seq.length());
		
		int[][] hashValues = getHashesStandard(seq);
		
		int count = 0;
		for (HashMap<Integer, ArrayList<KmerInfo>> hash : this.hashes)
		{
			ArrayList<KmerInfo> currList = hash.get(hashValues[0][count]);
			
			if (currList==null)
			{
				currList = new ArrayList<KmerInfo>();
				hash.put(hashValues[0][count], currList);
			}
			
			currList.add(new KmerInfo(seq.getId(), hashValues[1][count]));
			
			count++;
		}
	}
	
	public List<MatchResult> findMatches(Sequence seq, double readError)
	{
		HashMap<SequenceId, HitInfo> matchHitMap = new HashMap<SequenceId, HitInfo>();
		
		//go through each hash and collect hits
		int[][] hashes = getHashesStandard(seq);
		//int[][] hashes = getHashesCyclic(seq);
			
		for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
		{
			ArrayList<KmerInfo> currentHashMatchList = this.hashes.get(hashIndex).get(hashes[0][hashIndex]);
			
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
						if (hashes[1][hashIndex]>=0 && matchedId.getPosition()>=0)
							matchHitMap.put(matchedId.getId(), new HitInfo(hashes[1][hashIndex], matchedId.getPosition()));
						else
							matchHitMap.put(matchedId.getId(), new HitInfo());
					}
					else
					{
						//if this kmer is unique then use it for positioning
						if (hashes[1][hashIndex]>=0 && matchedId.getPosition()>=0)
							currentHitInfo.addUniqueHit(hashes[1][hashIndex], matchedId.getPosition());
						else
							currentHitInfo.addNonUniqueHit();
					}
				}
			}
		}
				
		//compute the proper counts for all sets and remove below threshold
		ArrayList<MatchResult> matches = new ArrayList<MatchResult>();
		
		//compute the probability that a matching kmer would match under error
		double probKmerEqual = Math.pow(1.0-readError, (double)this.kmerSize);
		
		for (SequenceId id : matchHitMap.keySet())
		{
			//get the hit info
			HitInfo hit = matchHitMap.get(id);
			
			//get the count
			double countPercent = (double)hit.count/(double)this.numHashes;
			
			//get the sequence length ratios
			double len1 = (double)this.seqLengths.get(id);
			double len2 = (double)seq.length();
			
			//compute the maximum number of kmers that is possible to have in common
			double maxSharedKmers = Math.min(len1-this.kmerSize, len2-this.kmerSize);
			
			//compute the scale to adjust the score based on sequence sizes and how many they have in common
			double adjScale = (len1+len2-maxSharedKmers*probKmerEqual)/maxSharedKmers;
			
			//adjust the score based on lengths
			double score = countPercent*adjScale;
			
			if (score>probKmerEqual*.9)
				matches.add(new MatchResult(seq.getId(), id, score, hit.fromShift()));
		}
		
		//sort the result
		Collections.sort(matches);
		
		return matches;
	}
	
	public void addData(FastaData data)
	{
		for (Sequence seq : data.getSequences())
		{
			addSequence(seq);
		}
	}
	
	/*
	private int[][] getHashesCyclic(Sequence seq)
	{
		//allocate the new hashes
		int[][] hashes = new int[2][this.numHashes];
		Arrays.fill(hashes[0], Integer.MAX_VALUE);
		
		//get the string
		String s = seq.getString();
		
		CyclicHash ch = new CyclicHash(this.kmerSize);
		
		//put initial values on
		int k = 0;
		for(; k<this.kmerSize;k++) 
		{
			ch.eat(s.charAt(k));
		}
				
		//get the first value
		int rollinghash = ch.eat(s.charAt(k)); // the first or last 32-(n-1) bits are 
		// do something with the hash value
						
		//start rolling
		for(;k<s.length();k++) 
		{
			//System.out.println(k-this.kmerSize);
			rollinghash = ch.update(s.charAt(k-this.kmerSize), s.charAt(k))<<this.kmerSize;
			
			char salt = ' ';
			ch.eat(s.charAt(k));
			
			for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
			{
				salt = (char)hashIndex;
			
				if (rollinghash<hashes[0][hashIndex])
				{
					hashes[0][hashIndex] = rollinghash;
				}
				if (rollinghash==hashes[0][hashIndex])
				{
					hashes[1][hashIndex] = -1;
				}
			}
		}
		
		return hashes;
	}
	*/
	
	private int[][] getHashesStandard(Sequence seq)
	{
		//allocate the new hashes
		int[][] hashes = new int[2][this.numHashes];
		Arrays.fill(hashes[0], Integer.MAX_VALUE);
		
		for (int kmerIndex=0; kmerIndex<seq.numKmers(this.kmerSize); kmerIndex++)
		{
			String kmer = seq.getKmer(kmerIndex, this.kmerSize);
			
			for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
			{
				long salt = hashIndex*10;
				
				//compute the hash
				int hash = (int)hashBasic(kmer, salt);
				
				//store minimum hash
				if (hash < hashes[0][hashIndex])
				{
					hashes[0][hashIndex] = hash;
					hashes[1][hashIndex] = kmerIndex;
				}
				else
				if (hash == hashes[0][hashIndex])
					hashes[1][hashIndex] = -1;
			}
		}
		
		return hashes;
	}
}
