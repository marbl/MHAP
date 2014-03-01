package com.secret.fastalign.general;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.minhash.SequenceMinHashes;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public abstract class AbstractHashSearch
{
	private final AtomicLong matchesProcessed;
	private final AtomicLong sequencesSearched;

	protected final int numThreads;
	protected final int minStoreLength;
	private final boolean storeResults;
	protected static BufferedWriter outWriter = new BufferedWriter(new OutputStreamWriter(System.out), 8*1024*1024);

	public AbstractHashSearch(int numThreads, int minStoreLength, boolean storeResults)
	{
		this.numThreads = numThreads;
		this.minStoreLength = minStoreLength;
		this.storeResults = storeResults;
		this.matchesProcessed = new AtomicLong();
		this.sequencesSearched = new AtomicLong();
	}
	
	protected void addData(final SequenceMinHashStreamer data)
	{
		//figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);
		
		final AtomicInteger counter = new AtomicInteger();
	  for (int iter=0; iter<this.numThreads; iter++)
		{
			Runnable task = new Runnable()
			{					
				@Override
				public void run()
				{
					try
					{
			    	SequenceMinHashes seqHashes = data.dequeue(false);
				    while(seqHashes != null)
				    {
			    		addSequence(seqHashes);

			    		int currCount = counter.incrementAndGet();
				    	if (currCount%5000==0)
				    		System.err.println("Current # sequences stored: "+currCount+"...");
				    	
				    	seqHashes = data.dequeue(false);
				    }
			    }
					catch (IOException e)
					{
						throw new FastAlignRuntimeException(e);
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

	protected abstract boolean addSequence(SequenceMinHashes seqHashes);
	
	public ArrayList<MatchResult> findMatches(final double acceptScore)
	{
		//figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);
		
		//allocate the storage and get the list of valeus
		final ArrayList<MatchResult> combinedList = new ArrayList<MatchResult>();
		final ConcurrentLinkedQueue<SequenceId> seqList = new ConcurrentLinkedQueue<SequenceId>(getStoredForwardSequenceIds());
	
		//for each thread create a task
	  for (int iter=0; iter<this.numThreads; iter++)
		{
			Runnable task = new Runnable()
			{ 			
				@Override
				public void run()
				{
	  			List<MatchResult> localMatches = new ArrayList<MatchResult>();

      		//get next sequence
	  			SequenceId nextSequence = seqList.poll();

	  			while (nextSequence!=null)
			    {		    		
		    		SequenceMinHashes sequenceHashes = getStoredSequenceHash(nextSequence);
		    		
		    		//only search the forward sequences
	      		localMatches.addAll(findMatches(sequenceHashes, acceptScore, true));
	      		
	      		//record search
	      		AbstractHashSearch.this.sequencesSearched.getAndIncrement();
		    		
	      		//get next sequence
		    		nextSequence = seqList.poll();

		    		//output stored results
		    		if (nextSequence==null || localMatches.size()>20000)
		    		{
			    		//count the number of matches
			  			AbstractHashSearch.this.matchesProcessed.getAndAdd(localMatches.size());
				    	
			  			if (AbstractHashSearch.this.storeResults)
			  			{	  			
				    		//combine the results
				    		synchronized (combinedList)
								{
									combinedList.addAll(localMatches);
								}
			  			}
			  			else
			  				outputResults(localMatches);
			  			
			  			localMatches.clear();
		    		}
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
	
	public ArrayList<MatchResult> findMatches(final SequenceMinHashStreamer data, final double acceptScore)
	{
		//figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);
		
		//allocate the storage and get the list of valeus
		final ArrayList<MatchResult> combinedList = new ArrayList<MatchResult>();
	
		//for each thread create a task
	  for (int iter=0; iter<this.numThreads; iter++)
		{
			Runnable task = new Runnable()
			{ 			
				@Override
				public void run()
				{
	  			List<MatchResult> localMatches = new ArrayList<MatchResult>();
	  			
	  			try
	  			{
		  			SequenceMinHashes sequenceHashes = data.dequeue(true);
	
		  			while (sequenceHashes!=null)
				    {		    		
			    		//only search the forward sequences
		      		localMatches.addAll(findMatches(sequenceHashes, acceptScore, false));
	
		      		//record search
		      		AbstractHashSearch.this.sequencesSearched.getAndIncrement();
		      		
		      		//get the sequence hashes
		      		sequenceHashes = data.dequeue(true);			    		

			    		//output stored results
			    		if (sequenceHashes==null || localMatches.size()>20000)
			    		{
				    		//count the number of matches
				  			AbstractHashSearch.this.matchesProcessed.getAndAdd(localMatches.size());
					    	
				  			if (AbstractHashSearch.this.storeResults)
				  			{	  			
					    		//combine the results
					    		synchronized (combinedList)
									{
										combinedList.addAll(localMatches);
									}
				  			}
				  			else
				  				outputResults(localMatches);
				  			
				  			localMatches.clear();
			    		}
				    }
		  		}
	  			catch (IOException e)
	  			{
	  				throw new FastAlignRuntimeException(e);
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

	protected abstract List<MatchResult> findMatches(SequenceMinHashes hashes, double acceptScore, boolean allToAll);

	public long getMatchesProcessed()
	{
		return this.matchesProcessed.get();
	}
	
	public abstract Collection<SequenceId> getStoredForwardSequenceIds();
	
	public abstract SequenceMinHashes getStoredSequenceHash(SequenceId id);

	protected void outputResults(List<MatchResult> matches)
	{
		if (this.storeResults || matches.isEmpty())
			return;
		
		try
		{
			synchronized (outWriter)
			{
				for (MatchResult currResult : matches)
				{
					outWriter.write(currResult.toString());
					outWriter.newLine();
				}

				outWriter.flush();
			}
		}
		catch (IOException e)
		{
			throw new FastAlignRuntimeException(e);
		}
	}

	public abstract int size();

	/**
	 * @return the sequencesSearched
	 */
	public long getNumberSequencesSearched()
	{
		return this.sequencesSearched.get();
	}

	/**
	 * @param outWriter the outWriter to set
	 */
	public static void setOutWriter(BufferedWriter outWriter)
	{
		AbstractHashSearch.outWriter = outWriter;
	}

}