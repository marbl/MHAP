package com.secret.fastalign.general;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.minhash.SequenceMinHashes;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Utils;

public final class SequenceMinHashStreamer
{
	private final DataInputStream buffInput;
	private final FastaData fastaData;
	private final HashSet<Integer> filter;
	private final int kmerSize;
	private final AtomicLong numberProcessed;
	private final AtomicLong numberSubSequencesProcessed;
	private final int numHashes;
	private final boolean readingFasta;
	private final ConcurrentLinkedQueue<SequenceMinHashes> sequenceHashList;
	private final int subKmerSize;
	private final int subSequenceSize;
	private boolean readClosed;

	public SequenceMinHashStreamer(String file) throws FileNotFoundException
	{
  	this.readingFasta = false;
		this.sequenceHashList = new ConcurrentLinkedQueue<SequenceMinHashes>();
		this.fastaData = null;
		this.kmerSize = 0;
		this.numHashes = 0;
		this.subSequenceSize = 0;
		this.subKmerSize = 0;
		this.filter = null;
		this.numberProcessed = new AtomicLong();
		this.numberSubSequencesProcessed = new AtomicLong();
		this.readClosed = false;

		this.buffInput = new DataInputStream(new BufferedInputStream(new FileInputStream(file), Utils.BUFFER_BYTE_SIZE));  	
	}
	
	public SequenceMinHashStreamer(String file, int kmerSize, int numHashes, int subSequenceSize, int subKmerSize, HashSet<Integer> filter) throws IOException
	{	
		this.sequenceHashList = new ConcurrentLinkedQueue<SequenceMinHashes>();
		this.fastaData = new FastaData(file);
		this.kmerSize = kmerSize;
		this.numHashes = numHashes;
		this.subSequenceSize = subSequenceSize;
		this.subKmerSize = subKmerSize;
		this.filter = filter;
		this.numberProcessed = new AtomicLong();
		this.numberSubSequencesProcessed = new AtomicLong();
		this.readingFasta = true;
		this.buffInput = null;
		this.readClosed = false;
	}
	
	public SequenceMinHashes dequeue(boolean fwdOnly) throws IOException
	{
		enqueue(fwdOnly);

		return this.sequenceHashList.poll();
	}
	
	private boolean enqueue(boolean fwdOnly) throws IOException
	{
		SequenceMinHashes seqHashes;
		if (this.readingFasta)
		{
			Sequence seq = this.fastaData.dequeue();
			while (seq!=null && fwdOnly && !seq.getId().isForward())
					seq = this.fastaData.dequeue();
			
			if (seq==null)
				return false;
			
			//compute the hashes
			seqHashes = new SequenceMinHashes(seq, this.kmerSize, this.numHashes, this.subSequenceSize, this.subKmerSize, false, this.filter);
			
			this.sequenceHashList.add(seqHashes);

			//increment counter
			this.numberProcessed.getAndIncrement();
			this.numberSubSequencesProcessed.getAndAdd(seqHashes.getMainHashes().numSubSequences());

			int numProcessed = getNumberProcessed();
	  	if (numProcessed%5000==0)
	  		System.err.println("Current # sequences loaded and processed from file: "+numProcessed+"...");

	  	if (!fwdOnly)
			{
				//compute the hashes
				seqHashes = new SequenceMinHashes(seq.getReverseCompliment(), this.kmerSize, this.numHashes, this.subSequenceSize, this.subKmerSize, false, this.filter);
				
				this.sequenceHashList.add(seqHashes);
	
				//increment counter
				this.numberProcessed.getAndIncrement();
				this.numberSubSequencesProcessed.getAndAdd(seqHashes.getMainHashes().numSubSequences());

				numProcessed = getNumberProcessed();
				if (numProcessed%5000==0)
		  		System.err.println("Current # sequences loaded and processed from file: "+numProcessed+"...");
			}
		}
		else
		{			
			//read the binary file
			seqHashes = readFromBinary();
			while (seqHashes!=null && fwdOnly && !seqHashes.getSequenceId().isForward())
			{
				seqHashes = readFromBinary();
			}
			
			//do nothing and return
			if (seqHashes==null)
				return false;
			
			this.sequenceHashList.add(seqHashes);
	
			//increment counter
			this.numberProcessed.getAndIncrement();
			this.numberSubSequencesProcessed.getAndAdd(seqHashes.getMainHashes().numSubSequences());

			int numProcessed = getNumberProcessed();
	  	if (numProcessed%5000==0)
	  		System.err.println("Current # sequences loaded and processed from file: "+numProcessed+"...");
		}
		
		return true;
	}
	
	public synchronized void enqueueFullFile(final boolean fwdOnly, int numThreads) throws IOException
	{
		//figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);
		
		//for each thread create a task
	  for (int iter=0; iter<numThreads; iter++)
		{
			Runnable task = new Runnable()
			{ 			
				@Override
				public void run()
				{
	        try
					{
						while (enqueue(fwdOnly))
						{
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
	
	public long getFastaProcessed()
	{
		if (this.fastaData==null)
			return 0;
		
		return this.fastaData.getNumberProcessed();
	}
	
	public int getNumberProcessed()
	{
		return this.numberProcessed.intValue();
	}
	
	public int getNumberSubSequencesProcessed()
	{
		return this.numberSubSequencesProcessed.intValue();
	}
	
	private SequenceMinHashes readFromBinary() throws IOException
	{
		byte[] byteArray = null;
		synchronized (this.buffInput)
		{
			if (this.readClosed)
				return null;
			
			try
			{
				//get the size in bytes
				int byteSize = this.buffInput.readInt();
				
				//allocate the array
				byteArray = new byte[byteSize];
				
				//read that many bytes
				this.buffInput.read(byteArray);
			}
			catch(EOFException e)
			{
	  		this.buffInput.close();
	  		this.readClosed = true;			
	  		
	  		return null;
			}
		}
			
		//get as byte array stream
  	SequenceMinHashes seqHashes = SequenceMinHashes.fromByteStream(new DataInputStream(new ByteArrayInputStream(byteArray)));

  	return seqHashes;
  }
	
	public void writeToBinary(String file, final boolean fwdOnly, int numThreads) throws IOException
	{
    OutputStream output = null;
    try 
    {
      output = new BufferedOutputStream(new FileOutputStream(file), Utils.BUFFER_BYTE_SIZE);
      final OutputStream finalOutput = output;
      
  		//figure out number of cores
  		ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);
  		
  		//for each thread create a task
  	  for (int iter=0; iter<numThreads; iter++)
  		{
  			Runnable task = new Runnable()
  			{ 			
  				@Override
  				public void run()
  				{
  	        SequenceMinHashes seqHashes;

  	        try
						{
							seqHashes = dequeue(fwdOnly);
	  	        while (seqHashes!=null)
	  	        {
  	        		byte[] byteArray = seqHashes.getAsByteArray();
  	        		int arraySize = byteArray.length;
  	        		
  	        		//store the size as byte array
  	        		byte[] byteSize = new byte[] {
  	                (byte)(arraySize >>> 24),
  	                (byte)(arraySize >>> 16),
  	                (byte)(arraySize >>> 8),
  	                (byte)arraySize};
  	        		
	  	        	synchronized (finalOutput)
								{	  	        		
	  	        		finalOutput.write(byteSize);
		  	        	finalOutput.write(byteArray);									
								}
	  	        	
	  	        	seqHashes = dequeue(fwdOnly);
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
  	  
  	  finalOutput.flush();
    }
    finally 
    {
      output.close();
    }
  }
}
