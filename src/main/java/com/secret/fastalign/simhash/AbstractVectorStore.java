package com.secret.fastalign.simhash;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import com.secret.fastalign.data.FastaData;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.data.SequenceId;
import com.secret.fastalign.main.MatchResult;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Pair;

public abstract class AbstractVectorStore<T extends VectorHash<T>>
{
	protected static final int SUB_KMER_SIZE = 4;
	protected static final int SUB_WORD_SIZE = 2;
	protected static final int SUB_STRING_SIZE = 100;
	
	protected final int kmerSize;
	protected final int numWords;
	protected final ConcurrentHashMap<SequenceId, SequenceVectorHashes<T>> sequenceVectorsHash;
	
	public AbstractVectorStore(int kmerSize, int numWords) 
	{
		this.kmerSize = kmerSize;
		this.numWords = numWords;

		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, SequenceVectorHashes<T>>();
	}
	

	public void addData(final FastaData data)
	{
  	//figure out number of cores
  	final int numThreads = Runtime.getRuntime().availableProcessors()*2;
  	ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);
  	
    for (int iter=0; iter<numThreads; iter++)
  	{
    	final int currThread = iter;
  		Runnable task = new Runnable()
			{					
				@Override
				public void run()
				{
			    for (int currIter=currThread; currIter<data.size(); currIter+=numThreads)
			    {
			    	boolean success = addSequence(data.getSequence(currIter));
			    	
			    	if (!success)
			    		throw new FastAlignRuntimeException("Error: Unable to add duplicate sequence "+data.getSequence(currIter).getId()+".");
			    }
				}
			};
  	
    	//enqueue the task
			execSvc.execute(task);					
  	}
  	
		//shutdown the service
    execSvc.shutdown();
    try
		{
			execSvc.awaitTermination(365L, TimeUnit.DAYS);
		} 
    catch (InterruptedException e) 
    {
    	execSvc.shutdownNow();
    	throw new FastAlignRuntimeException("Unable to finish all tasks.");
    }
	}
	
	public boolean addDirectionalSequence(Sequence seq)
	{
		//put the result into the hashmap
		SequenceVectorHashes<T> simHash = this.sequenceVectorsHash.put(seq.getId(), getVectorHash(seq, this.kmerSize, this.numWords));
		
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
	
	public boolean addSequence(Sequence seq)
	{
		//add forward sequence
		boolean success = addDirectionalSequence(seq);
		
		//add reverse sequence
		if (success)
		{
			Sequence reverse = seq.getReverseCompliment();		
			success = addDirectionalSequence(reverse);
		}

		return success;
	}
	
	public abstract SequenceVectorHashes<T> getVectorHash(Sequence seq, int kmerSize, int numWords);
	
	public List<MatchResult> findMatches(Sequence seq, double acceptScore)
	{
		//see if already have the sequence stored
		SequenceVectorHashes<T> simHash = this.sequenceVectorsHash.get(seq.getId());
		
		//if not get the sequence
		if (simHash==null)
			simHash = getVectorHash(seq, this.kmerSize, this.numWords);
		
		return findMatches(simHash, acceptScore);
	}
	
	public ArrayList<MatchResult> findMatches(final double acceptScore)
	{
  	//figure out number of cores
  	final int numThreads = Runtime.getRuntime().availableProcessors()*2;
  	ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);
  	
  	//allocate the storage and get the list of valeus
  	final Collection<SequenceVectorHashes<T>> storedHashes = this.sequenceVectorsHash.values();
  	final ArrayList<MatchResult> combinedList = new ArrayList<MatchResult>();

  	//for each thread create a task
    for (int iter=0; iter<numThreads; iter++)
  	{
    	final int currThread = iter;
  		Runnable task = new Runnable()
			{ 			
				@Override
				public void run()
				{
	  			Iterator<SequenceVectorHashes<T>> iterator = storedHashes.iterator();
	  			List<MatchResult> localMatches = new ArrayList<MatchResult>();
	  			
	  			//skip to the initial value
		    	for (int skipIter=0; skipIter<currThread; skipIter++)
		    	{
		    		if (iterator.hasNext())
		    			iterator.next();
		    		else
		    			break;
		    	}

		    	//if there are values left
		    	while (iterator.hasNext())
			    {
		    		SequenceVectorHashes<T> nextSequence = iterator.next();
		    		
		    		//only search the forward sequences
		    		if (nextSequence.getSequenceId().isForward())
		      		localMatches.addAll(findMatches(nextSequence, acceptScore));

	      		//skip the required number of iterations
			    	for (int skipIter=0; skipIter<numThreads-1; skipIter++)
			    	{
			    		if (iterator.hasNext())
			    			iterator.next();
			    		else
			    			break;
			    	}
			    }
		    	
	    		//combine the results
	    		synchronized (combinedList)
					{
						combinedList.addAll(localMatches);
					}

				}
			};
  	
    	//enqueue the task
			execSvc.execute(task);					
  	}
  	
		//shutdown the service
    execSvc.shutdown();
    try
		{
			execSvc.awaitTermination(365L, TimeUnit.DAYS);
		} 
    catch (InterruptedException e) 
    {
    	execSvc.shutdownNow();
    	throw new FastAlignRuntimeException("Unable to finish all tasks.");
    }
		
		return combinedList;
	}

	public int size()
	{
		return this.sequenceVectorsHash.size();
	}
	
	public List<MatchResult> findMatches(SequenceVectorHashes<T> seqHash, double acceptScore)
	{		
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (SequenceVectorHashes<T> hash : this.sequenceVectorsHash.values())
		{
			if (seqHash.getSequenceId().getHeaderId()==hash.getSequenceId().getHeaderId())
				continue;
			
			//compute the initial score
			double score = seqHash.jaccord(hash);
			
			if (score>=acceptScore)
			{
				Pair<Double,Integer> result = seqHash.orderedScore(hash);
				score = result.x;
				int shift = result.y*SUB_STRING_SIZE;
				
				results.add(new MatchResult(seqHash.getSequenceId(), hash.getSequenceId(), score, shift));
			}
		}
		
		//not decided if should sort
		Collections.sort(results);
		
		return results;
	}
}
