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

public abstract class AbstractVectorStore<T extends VectorHash<T>>
{
	protected final int kmerSize;
	protected final ConcurrentHashMap<SequenceId, T> sequenceVectorsHash;
	protected final int numWords;
	
	public AbstractVectorStore(int kmerSize, int numWords) 
	{
		this.kmerSize = kmerSize;
		this.sequenceVectorsHash = new ConcurrentHashMap<SequenceId, T>();
		this.numWords = numWords;
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
			    	T simHashPrev = addSequence(data.getSequence(currIter));
			    	
			    	if (simHashPrev!=null)
			    		System.err.println("Warning: Duplicate sequence id "+simHashPrev.getSequenceId()+".");
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
	
	public T addSequence(Sequence seq)
	{
		Sequence reverse = seq.getReverseCompliment();

		//store the complement
		T simHash = this.sequenceVectorsHash.get(reverse.getId());
		
		if (simHash==null)
			simHash = this.sequenceVectorsHash.put(reverse.getId(), getVectorHash(reverse));
		
		simHash = this.sequenceVectorsHash.get(seq.getId());
		
		if (simHash==null)
			simHash = this.sequenceVectorsHash.put(seq.getId(), getVectorHash(seq));
					
		return simHash;
	}
	
	public abstract T getVectorHash(Sequence seq);
	
	public List<MatchResult> findMatches(Sequence seq, double acceptScore)
	{
		//see if already have the sequence stored
		T simHash = this.sequenceVectorsHash.get(seq.getId());
		
		//if not get the sequence
		if (simHash==null)
			simHash = getVectorHash(seq);
		
		return findMatches(simHash, acceptScore);
	}
	
	public ArrayList<MatchResult> findMatches(final double acceptScore)
	{
  	//figure out number of cores
  	final int numThreads = Runtime.getRuntime().availableProcessors()*2;
  	ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);
  	
  	//allocate the storage and get the list of valeus
  	final Collection<T> storedHashes = this.sequenceVectorsHash.values();
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
	  			Iterator<T> iterator = storedHashes.iterator();
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
		    		T nextSequence = iterator.next();
		    		
		    		//only search the forward sequences
		    		if (nextSequence.getSequenceId().isForward())
		    		{
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

	
	public List<MatchResult> findMatches(T seqHash, double acceptScore)
	{		
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (T hash : this.sequenceVectorsHash.values())
		{
			if (seqHash.getSequenceId().getLongId()==(hash.getSequenceId().getLongId()))
				continue;
			
			double score = seqHash.correlation(hash);
			
			if (score>acceptScore)
				results.add(new MatchResult(seqHash.getSequenceId(), hash.getSequenceId(), score, 0));
		}
		
		//not decided if should sort
		Collections.sort(results);
		
		return results;
	}
}
