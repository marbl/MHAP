package com.secret.fastalign.general;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.util.Iterator;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.ReadBuffer;
import com.secret.fastalign.utils.Utils;

public abstract class AbstractSequenceHashStreamer<H extends SequenceHashes>
{
	private final FastaData fastaData;
	private final AtomicLong numberProcessed;
	private final boolean readingFasta;
	private final ConcurrentLinkedQueue<H> sequenceHashList;
	
	public AbstractSequenceHashStreamer(FastaData data, boolean readingFasta)
	{
		this.fastaData = data;
		this.readingFasta = readingFasta;
		this.sequenceHashList = new ConcurrentLinkedQueue<H>();
		this.numberProcessed = new AtomicLong();
	}
	
	public H dequeue(boolean fwdOnly, ReadBuffer buf) throws IOException
	{
		enqueue(fwdOnly, buf);

		return this.sequenceHashList.poll();
	}
	
	private boolean enqueue(boolean fwdOnly, ReadBuffer buf) throws IOException
	{
		H seqHashes;
		if (this.readingFasta)
		{
			Sequence seq = this.fastaData.dequeue();
						
			//compute the hashes
			seqHashes = null;
			if (seq!=null)
				seqHashes = getHashes(seq);		

			if (seqHashes==null)
				return false;
			processAddition(seqHashes);
			
			this.sequenceHashList.add(seqHashes);
	
			//fasta files are all fwd
			if (!fwdOnly)
			{
				//compute the hashes
				seqHashes = getHashes(seq.getReverseCompliment());
				
				this.sequenceHashList.add(seqHashes);
				processAddition(seqHashes);
			}
		}
		else
		{
			//read the binary file
			seqHashes = readFromBinary(buf, fwdOnly);
			while (seqHashes!=null && fwdOnly && !seqHashes.getSequenceId().isForward())
			{
				seqHashes = readFromBinary(buf, fwdOnly);
			}
			
			//do nothing and return
			//record
			if (seqHashes==null)
				return false;

			processAddition(seqHashes);
			
			this.sequenceHashList.add(seqHashes);

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
					ReadBuffer buf = new ReadBuffer();
					
	        try
					{
						while (enqueue(fwdOnly, buf))
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
	
	public abstract H getHashes(Sequence seq);
	
	protected void processAddition(H seqHashes)
	{
		//increment counter
		this.numberProcessed.getAndIncrement();

		int numProcessed = getNumberProcessed();
  	if (numProcessed%5000==0)
  		System.err.println("Current # sequences loaded and processed from file: "+numProcessed+"...");
	}
	
	public Iterator<H> getDataIterator()
	{
		return this.sequenceHashList.iterator();
	}
	
	public int getNumberProcessed()
	{
		return this.numberProcessed.intValue();
	}
	
	public int getNumberSubSequencesProcessed()
	{
		return this.getNumberProcessed();
	}
	
	protected abstract H readFromBinary(ReadBuffer buf, boolean fwdOnly) throws IOException;

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
  	        H seqHashes;
  	        ReadBuffer buf = new ReadBuffer();
  	        

  	        try
						{
							seqHashes = dequeue(fwdOnly, buf);
	  	        while (seqHashes!=null)
	  	        {
  	        		byte[] byteArray = seqHashes.getAsByteArray();
  	        		int arraySize = byteArray.length;
  	        		byte isFwd = seqHashes.getSequenceId().isForward() ? (byte)1 :(byte)0;
  	        		
  	        		//store the size as byte array
  	        		byte[] byteSize = ByteBuffer.allocate(5).put(isFwd).putInt(arraySize).array();
  	        		
	  	        	synchronized (finalOutput)
								{
	  	        		finalOutput.write(byteSize);
		  	        	finalOutput.write(byteArray);									
								}
	  	        	
	  	        	seqHashes = dequeue(fwdOnly, buf);
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
