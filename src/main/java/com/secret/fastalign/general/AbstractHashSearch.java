package com.secret.fastalign.general;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import com.secret.fastalign.utils.FastAlignRuntimeException;

public abstract class AbstractHashSearch<H extends AbstractSequenceHashes<H>, T extends AbstractReducedSequence<H,T>>
{

	protected final int kmerSize;
	protected final int numWords;

	public AbstractHashSearch(int kmerSize, int numWords)
	{
		this.kmerSize = kmerSize;
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
	
	public abstract Collection<SequenceId> getStoredForwardSequenceIds();
	
	public abstract T getStoredSequenceHash(SequenceId id);

	public abstract boolean addDirectionalSequence(Sequence seq);

	public List<MatchResult> findMatches(Sequence seq, double acceptScore)
	{
		//see if already have the sequence stored
		T hashes = getStoredSequenceHash(seq.getId());
		
		//if not get the sequence
		if (hashes==null)
			hashes = getSequenceHash(seq);
		
		return findMatches(hashes, acceptScore, false);
	}

	public abstract List<MatchResult> findMatches(T hashes, double acceptScore, boolean allToAll);
	
	public abstract T getSequenceHash(Sequence seq);

	public ArrayList<MatchResult> findMatches(final double acceptScore)
	{
		//figure out number of cores
		final int numThreads = Runtime.getRuntime().availableProcessors()*2;
		ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);
		
		//allocate the storage and get the list of valeus
		final Collection<SequenceId> storedHashes = getStoredForwardSequenceIds();
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
	  			Iterator<SequenceId> iterator = storedHashes.iterator();
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
		    		//get the sequence hashes
		    		SequenceId nextSequence = iterator.next();
		    		T sequenceHashes = getStoredSequenceHash(nextSequence);
		    		
		    		//only search the forward sequences
	      		localMatches.addAll(findMatches(sequenceHashes, acceptScore, true));
	
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

	public abstract int size();

}