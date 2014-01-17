package com.secret.fastalign.minhash;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;
import com.secret.fastalign.data.FastaData;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.data.SequenceId;
import com.secret.fastalign.main.MatchResult;
import com.secret.fastalign.utils.RabinKarpHash;

public final class MinHash
{
	private final class HitInfo
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
	
	public final static int SUPER_SHINGLE_SIZE = 2;
	
	private final ArrayList<HashMap<Integer, ArrayList<KmerInfo>>> hashes;
	private final int kmerSize;
	private final int numHashes;
	private final ArrayList<Integer> seeds;
	private final HashMap<SequenceId, Integer> seqLengths;
	
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
	
	private static final long hashBasic(final String kmer, final long salt) 
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
	
	public MinHash(int numHashes, int kmerSize)
	{
		this.kmerSize = kmerSize;
		this.numHashes = numHashes;
	  this.hashes = new ArrayList<HashMap<Integer, ArrayList<KmerInfo>>>();
	  this.seqLengths = new HashMap<SequenceId, Integer>();
	  
	  this.seeds = new ArrayList<Integer>();
	  
	  Random rand = new Random(1);
	  
	  //allocate the hashtables
	  for (int iter=0; iter<this.numHashes*SUPER_SHINGLE_SIZE; iter++)
	  {
	  	this.seeds.add(rand.nextInt());
	  	this.hashes.add(new HashMap<Integer, ArrayList<KmerInfo>>());
	  }
	}
	
	public void addData(FastaData data)
	{
		for (Sequence seq : data.getSequences())
		{
			addSequence(seq);
		}
	}
	
	public void addSequence(Sequence seq)
	{
		//record the length
		this.seqLengths.put(seq.getId(), seq.length());
		
		int[][] hashValues = getHashes(seq);
		
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
			
			if (count>=this.numHashes)
				break;
		}
	}
	
	public List<MatchResult> findMatches(Sequence seq, double probEqualKmerMatch)
	{
		HashMap<SequenceId, HitInfo> matchHitMap = new HashMap<SequenceId, HitInfo>();
		
		//get seq length
		int len1 = seq.length();

		//go through each hash and collect hits
		int[][] hashes = getHashes(seq);
			
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
		
		for (SequenceId id : matchHitMap.keySet())
		{
			//do not store matches to yourself
			if (id.equals(seq.getId()))
				continue;
			
			//get the hit info
			HitInfo hit = matchHitMap.get(id);
			
			double countPercent = (double)hit.count/(double)this.numHashes;
			
			//get the sequence length ratios
			int len2 = this.seqLengths.get(id);
						
			double jaccardAccept = jaccardRate(probEqualKmerMatch, len1, len2, 1.0);
			
			//System.out.println(""+countPercent +" "+jaccardAccept);
			
			//if (countPercent>=jaccardAccept)
			{
				//double jaccardRatio = countPercent/jaccardRate(probEqualKmerMatch, len1, len2, 1.0);
				//matches.add(new MatchResult(seq.getId(), id, jaccardRatio, hit.fromShift()));
			  matches.add(new MatchResult(seq.getId(), id, countPercent*(double)(Math.max(len1,len2)/Math.min(len1,len2)), hit.fromShift()));
			}
		}
		
		//sort the result
		Collections.sort(matches);
		
		return matches;
	}
	
	public int[][] getHashes(Sequence seq)
	{
		//return getHashesCyclic(seq);
		return getHashesStandard(seq);
	}
	
	public int[][] getHashesCyclic(Sequence seq)
	{
		//allocate the new hashes
		int[][] hashes = new int[2][this.numHashes*SUPER_SHINGLE_SIZE];
		int[][] hashesFinal = new int[2][this.numHashes];

		Arrays.fill(hashes[0], Integer.MAX_VALUE);
		
		//get the string
		String seqString = seq.getString();
		char[] seqArray = seqString.toCharArray();
		
		RabinKarpHash ch = new RabinKarpHash(this.kmerSize);
		
		int kmerIndex;
		for(kmerIndex = 0; kmerIndex<this.kmerSize-1; kmerIndex++) 
		{
			ch.eat(seqArray[kmerIndex]);
		}
		
		int hash = ch.eat(seqArray[kmerIndex++]);
		
		//set the value
		for (int hashIndex=0; hashIndex<hashes[0].length; hashIndex++)
		{
			int seed = this.seeds.get(hashIndex);

			hashes[0][hashIndex] = hash*22695477+seed;
			hashes[1][hashIndex] = kmerIndex;
		}
		
		//start rolling
		for(; kmerIndex<seqArray.length; kmerIndex++) 
		{
			//System.out.println(k-this.kmerSize);
			hash = ch.update(seqArray[kmerIndex-this.kmerSize], seqArray[kmerIndex]);
			
			for (int hashIndex=0; hashIndex<hashes[0].length; hashIndex++)
			{
				int seed = this.seeds.get(hashIndex);

				//int currHash = this.hf.newHasher().putInt(hash).putInt(seed).hashCode();
				int currHash = hash*22695477+seed;
				
				if (currHash <= hashes[0][hashIndex])
				{
					hashes[0][hashIndex] = currHash;
					
					//record the position of the kmer
					if (currHash == hashes[0][hashIndex])
						hashes[1][hashIndex] = -1;
					else
						hashes[1][hashIndex] = kmerIndex;
				}
			}
			
		}
		
		for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
		{
			hashesFinal[0][hashIndex] = hashes[0][SUPER_SHINGLE_SIZE*hashIndex];
			for (int superSize=1; superSize<SUPER_SHINGLE_SIZE; superSize++)
				hashesFinal[0][hashIndex] = hashesFinal[0][hashIndex]*22695477+hashes[0][SUPER_SHINGLE_SIZE*hashIndex+superSize];
			hashesFinal[1][hashIndex] = hashes[1][SUPER_SHINGLE_SIZE*hashIndex+0];
		}
		
		return hashesFinal;
	}
	
	public int[][] getHashesStandard(Sequence seq)
	{
		//allocate the new hashes
		int[][] hashes = new int[2][this.numHashes*SUPER_SHINGLE_SIZE];
		int[][] hashesFinal = new int[2][this.numHashes];
		
		Arrays.fill(hashes[0], Integer.MAX_VALUE);
		
		HashFunction hasher = Hashing.murmur3_32(0);
		
		//get the string
		String seqString = seq.getString();
		char[] seqArray = seqString.toCharArray();
			
		//start rolling
		for(int kmerIndex = 0; kmerIndex<seqArray.length-kmerIndex; kmerIndex++) 
		{
			//System.out.println(k-this.kmerSize);
			for (int hashIndex=0; hashIndex<hashes[0].length; hashIndex++)
			{
				int currHash = Integer.MAX_VALUE;
				for (int subKmerIndex = 0; subKmerIndex<this.kmerSize; subKmerIndex++)
				{
					int hash = hasher.newHasher().putChar(seqArray[kmerIndex+subKmerIndex]).putInt(subKmerIndex).putInt(hashIndex).hash().asInt();
					
					if (hash<currHash)
						currHash = hash;
				}

				if (currHash <= hashes[0][hashIndex])
				{
					hashes[0][hashIndex] = currHash;
					
					//record the position of the kmer
					if (currHash == hashes[0][hashIndex])
						hashes[1][hashIndex] = -1;
					else
						hashes[1][hashIndex] = kmerIndex;
				}
			}			
	}
		
		for (int hashIndex=0; hashIndex<this.numHashes; hashIndex++)
		{
			hashesFinal[0][hashIndex] = hashes[0][SUPER_SHINGLE_SIZE*hashIndex];
			for (int superSize=1; superSize<SUPER_SHINGLE_SIZE; superSize++)
				hashesFinal[0][hashIndex] = hashesFinal[0][hashIndex]*22695477+hashes[0][SUPER_SHINGLE_SIZE*hashIndex+superSize];
			hashesFinal[1][hashIndex] = hashes[1][SUPER_SHINGLE_SIZE*hashIndex+0];
		}
		
		return hashesFinal;
	}
	
	public double maxPercentInCommon(int len1, int len2)
	{
		double numKmers1 = (double)(len1-this.kmerSize);
		double numKmers2 = (double)(len2-this.kmerSize);
		
		double numCommon = (double)Math.min(numKmers1,numKmers2);
		
		return numCommon/(numKmers1+numKmers2);
		
	}
	
	public double jaccardRate(double probKmerMatches, int len1, int len2, double acceptRatio)
	{
		double numKmers1 = (double)(len1-this.kmerSize);
		double numKmers2 = (double)(len2-this.kmerSize);
		
		double numCommon = (double)Math.min(numKmers1,numKmers2)*probKmerMatches;
		
		double jaccard = numCommon*acceptRatio/(numKmers1+numKmers2-numCommon*acceptRatio);
		
		return jaccard;
	}
}
