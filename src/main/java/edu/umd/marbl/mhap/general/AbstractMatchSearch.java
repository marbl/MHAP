/* 
 * MHAP package
 * 
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2014 by Konstantin Berlin and Sergey Koren
 * University Of Maryland
 * 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */
package edu.umd.marbl.mhap.general;

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

import edu.umd.marbl.mhap.sketch.SequenceSketch;
import edu.umd.marbl.mhap.sketch.SequenceSketchStreamer;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.ReadBuffer;
import edu.umd.marbl.mhap.utils.Utils;

public abstract class AbstractMatchSearch
{
	private final AtomicLong matchesProcessed;
	protected final int numThreads;

	private final AtomicLong sequencesSearched;
	private final boolean storeResults;

	public final static int NUM_ELEMENTS_PER_OUTPUT = 20000;
	protected final static BufferedWriter STD_OUT_BUFFER = new BufferedWriter(new OutputStreamWriter(System.out),
			Utils.BUFFER_BYTE_SIZE);

	public AbstractMatchSearch(int numThreads, boolean storeResults)
	{
		this.numThreads = numThreads;
		this.storeResults = storeResults;
		this.matchesProcessed = new AtomicLong();
		this.sequencesSearched = new AtomicLong();
	}

	protected void addData(final SequenceSketchStreamer data)
	{
		// figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);

		final AtomicInteger counter = new AtomicInteger();
		for (int iter = 0; iter < this.numThreads; iter++)
		{
			Runnable task = new Runnable()
			{
				@Override
				public void run()
				{
					try
					{
						ReadBuffer buf = new ReadBuffer();
						SequenceSketch seqHashes = data.dequeue(false, buf);
						while (seqHashes != null)
						{
							addSequence(seqHashes);

							int currCount = counter.incrementAndGet();
							if (currCount % 5000 == 0)
								System.err.println("Current # sequences stored: " + currCount + "...");

							seqHashes = data.dequeue(false, buf);
						}
					}
					catch (IOException e)
					{
						throw new MhapRuntimeException(e);
					}
				}
			};

			// enqueue the task
			execSvc.execute(task);
		}

		// shutdown the service
		execSvc.shutdown();
		try
		{
			execSvc.awaitTermination(365L, TimeUnit.DAYS);
		}
		catch (InterruptedException e)
		{
			execSvc.shutdownNow();
			throw new MhapRuntimeException("Unable to finish all tasks.");
		}
	}

	protected abstract boolean addSequence(SequenceSketch seqHashes);

	public ArrayList<MatchResult> findMatches()
	{
		// figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);

		// allocate the storage and get the list of valeus
		final ArrayList<MatchResult> combinedList = new ArrayList<MatchResult>();
		final ConcurrentLinkedQueue<SequenceId> seqList = new ConcurrentLinkedQueue<SequenceId>(
				getStoredForwardSequenceIds());

		// for each thread create a task
		for (int iter = 0; iter < this.numThreads; iter++)
		{
			Runnable task = new Runnable()
			{
				@Override
				public void run()
				{
					List<MatchResult> localMatches = new ArrayList<MatchResult>();

					// get next sequence
					SequenceId nextSequence = seqList.poll();

					while (nextSequence != null)
					{
						SequenceSketch sequenceHashes = getStoredSequenceHash(nextSequence);

						// only search the forward sequences
						localMatches.addAll(findMatches(sequenceHashes, true));

						// record search
						AbstractMatchSearch.this.sequencesSearched.getAndIncrement();

						// get next sequence
						nextSequence = seqList.poll();

						// output stored results
						if (nextSequence == null || localMatches.size() >= NUM_ELEMENTS_PER_OUTPUT)
						{
							// count the number of matches
							AbstractMatchSearch.this.matchesProcessed.getAndAdd(localMatches.size());

							if (AbstractMatchSearch.this.storeResults)
							{
								// combine the results
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

			// enqueue the task
			execSvc.execute(task);
		}

		// shutdown the service
		execSvc.shutdown();
		try
		{
			execSvc.awaitTermination(365L, TimeUnit.DAYS);
		}
		catch (InterruptedException e)
		{
			execSvc.shutdownNow();
			throw new MhapRuntimeException("Unable to finish all tasks.");
		}

		flushOutput();

		return combinedList;
	}

	protected abstract List<MatchResult> findMatches(SequenceSketch hashes, boolean toSelf);

	public ArrayList<MatchResult> findMatches(final SequenceSketchStreamer data) throws IOException
	{
		// figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);

		// allocate the storage and get the list of valeus
		final ArrayList<MatchResult> combinedList = new ArrayList<MatchResult>();

		// for each thread create a task
		for (int iter = 0; iter < this.numThreads; iter++)
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

						SequenceSketch sequenceHashes = data.dequeue(true, buf);

						while (sequenceHashes != null)
						{
							// only search the forward sequences
							localMatches.addAll(findMatches(sequenceHashes, false));

							// record search
							AbstractMatchSearch.this.sequencesSearched.getAndIncrement();

							// get the sequence hashes
							sequenceHashes = data.dequeue(true, buf);

							// output stored results
							if (sequenceHashes == null || localMatches.size() >= NUM_ELEMENTS_PER_OUTPUT)
							{
								// count the number of matches
								AbstractMatchSearch.this.matchesProcessed.getAndAdd(localMatches.size());

								if (AbstractMatchSearch.this.storeResults)
								{
									// combine the results
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
						throw new MhapRuntimeException(e);
					}
				}
			};

			// enqueue the task
			execSvc.execute(task);
		}

		// shutdown the service
		execSvc.shutdown();
		try
		{
			execSvc.awaitTermination(365L, TimeUnit.DAYS);
		}
		catch (InterruptedException e)
		{
			execSvc.shutdownNow();
			throw new MhapRuntimeException("Unable to finish all tasks.");
		}

		flushOutput();

		return combinedList;
	}

	protected void flushOutput()
	{
		try
		{
			STD_OUT_BUFFER.flush();
		}
		catch (IOException e)
		{
			throw new MhapRuntimeException(e);
		}
	}

	public long getMatchesProcessed()
	{
		return this.matchesProcessed.get();
	}

	/**
	 * @return the sequencesSearched
	 */
	public long getNumberSequencesSearched()
	{
		return this.sequencesSearched.get();
	}

	public abstract List<SequenceId> getStoredForwardSequenceIds();

	public abstract SequenceSketch getStoredSequenceHash(SequenceId id);

	protected void outputResults(List<MatchResult> matches)
	{
		if (this.storeResults || matches.isEmpty())
			return;

		try
		{
			synchronized (STD_OUT_BUFFER)
			{
				for (MatchResult currResult : matches)
				{
					STD_OUT_BUFFER.write(currResult.toString());
					STD_OUT_BUFFER.newLine();
				}

				STD_OUT_BUFFER.flush();
			}
		}
		catch (IOException e)
		{
			throw new MhapRuntimeException(e);
		}
	}

	public abstract int size();

}