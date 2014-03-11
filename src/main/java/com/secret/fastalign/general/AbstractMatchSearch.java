package com.secret.fastalign.general;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.ReadBuffer;

public abstract class AbstractMatchSearch<H extends SequenceHashes>
{
	private final AtomicLong matchesProcessed;
	private final AtomicLong sequencesSearched;

	protected final int numThreads;
	private final boolean storeResults;
	protected static BufferedWriter outWriter = new BufferedWriter(new OutputStreamWriter(System.out), 8*1024*1024);

	public AbstractMatchSearch(int numThreads, boolean storeResults)
	{
		this.numThreads = numThreads;
		this.storeResults = storeResults;
		this.matchesProcessed = new AtomicLong();
		this.sequencesSearched = new AtomicLong();
	}
	
	protected void addData(final AbstractSequenceHashStreamer<H> data)
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
						ReadBuffer buf = new ReadBuffer();
			    	H seqHashes = data.dequeue(false, buf);
				    while(seqHashes != null)
				    {
			    		addSequence(seqHashes);

			    		int currCount = counter.incrementAndGet();
				    	if (currCount%5000==0)
				    		System.err.println("Current # sequences stored: "+currCount+"...");
				    	
				    	seqHashes = data.dequeue(false, buf);
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

	protected abstract boolean addSequence(H seqHashes);
	
	public ArrayList<MatchResult> findMatches()
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
		    		H sequenceHashes = getStoredSequenceHash(nextSequence);
		    		
		    		//only search the forward sequences
	      		localMatches.addAll(findMatches(sequenceHashes, true));
	      		
	      		//record search
	      		AbstractMatchSearch.this.sequencesSearched.getAndIncrement();
		    		
	      		//get next sequence
		    		nextSequence = seqList.poll();

		    		//output stored results
		    		if (nextSequence==null || localMatches.size()>20000)
		    		{
			    		//count the number of matches
			  			AbstractMatchSearch.this.matchesProcessed.getAndAdd(localMatches.size());
				    	
			  			if (AbstractMatchSearch.this.storeResults)
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
	  
	  flushOutput();

	  return combinedList;
	}
	
	public ArrayList<MatchResult> findMatches(final AbstractSequenceHashStreamer<H> data) throws IOException
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
	  				ReadBuffer buf = new ReadBuffer();
	  				
		  			H sequenceHashes = data.dequeue(true, buf);
	
		  			while (sequenceHashes!=null)
				    {		    		
			    		//only search the forward sequences
		      		localMatches.addAll(findMatches(sequenceHashes, false));
	
		      		//record search
		      		AbstractMatchSearch.this.sequencesSearched.getAndIncrement();
		      		
		      		//get the sequence hashes
		      		sequenceHashes = data.dequeue(true, buf);			    		

			    		//output stored results
			    		if (sequenceHashes==null || localMatches.size()>20000)
			    		{
				    		//count the number of matches
				  			AbstractMatchSearch.this.matchesProcessed.getAndAdd(localMatches.size());
					    	
				  			if (AbstractMatchSearch.this.storeResults)
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
		
	  flushOutput();
	  
		return combinedList;	
	}

	protected abstract List<MatchResult> findMatches(H hashes, boolean toSelf);

	public long getMatchesProcessed()
	{
		return this.matchesProcessed.get();
	}
	
	public abstract List<SequenceId> getStoredForwardSequenceIds();
	
	public abstract H getStoredSequenceHash(SequenceId id);

	protected void flushOutput()
	{
		try
		{
			outWriter.flush();
		}
		catch (IOException e)
		{
			throw new FastAlignRuntimeException(e);
		}
	}

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
		AbstractMatchSearch.outWriter = outWriter;
	}

}