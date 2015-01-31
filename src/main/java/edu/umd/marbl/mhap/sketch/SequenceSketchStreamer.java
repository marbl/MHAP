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
package edu.umd.marbl.mhap.sketch;

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
import java.nio.ByteBuffer;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.ReadBuffer;
import edu.umd.marbl.mhap.utils.Utils;

public class SequenceSketchStreamer
{
	private final DataInputStream buffInput;
	private final FastaData fastaData;
	private final HashSet<Long> filter;
	private final KmerCounts kmerCounter;
	private final int kmerSize;
	private final AtomicLong numberProcessed;
	private final AtomicLong numberSubSequencesProcessed;
	private final int numHashes;
	private final int offset;
	private final boolean weighted;

	private final int orderedKmerSize;
	private boolean readClosed;
	private final boolean readingFasta;
	private final ConcurrentLinkedQueue<SequenceSketch> sequenceHashList;

	public SequenceSketchStreamer(String file, int offset) throws FileNotFoundException
	{
		this.fastaData = null;
		this.readingFasta = false;
		this.sequenceHashList = new ConcurrentLinkedQueue<SequenceSketch>();
		this.numberProcessed = new AtomicLong();
		this.kmerCounter = null;
		this.weighted = true;

		this.kmerSize = 0;
		this.numHashes = 0;
		this.orderedKmerSize = 0;
		this.filter = null;
		this.numberSubSequencesProcessed = new AtomicLong();
		this.readClosed = false;
		this.offset = offset;

		this.buffInput = new DataInputStream(new BufferedInputStream(new FileInputStream(file), Utils.BUFFER_BYTE_SIZE));
	}

	public SequenceSketchStreamer(String file, int kmerSize, int numHashes, int subSequenceSize, int orderedKmerSize,
			HashSet<Long> filter, KmerCounts kmerCounter, boolean weighted, int offset) throws IOException
	{
		this.fastaData = new FastaData(file, offset);
		this.readingFasta = true;
		this.sequenceHashList = new ConcurrentLinkedQueue<SequenceSketch>();
		this.numberProcessed = new AtomicLong();

		this.weighted = weighted;
		this.kmerCounter = kmerCounter;
		this.kmerSize = kmerSize;
		this.numHashes = numHashes;
		this.orderedKmerSize = orderedKmerSize;
		this.filter = filter;
		this.numberSubSequencesProcessed = new AtomicLong();
		this.buffInput = null;
		this.readClosed = false;
		this.offset = offset;
	}

	public SequenceSketch dequeue(boolean fwdOnly, ReadBuffer buf) throws IOException
	{
		enqueue(fwdOnly, buf);

		return this.sequenceHashList.poll();
	}

	private boolean enqueue(boolean fwdOnly, ReadBuffer buf) throws IOException
	{
		SequenceSketch seqHashes;
		if (this.readingFasta)
		{
			Sequence seq = this.fastaData.dequeue();

			// compute the hashes
			seqHashes = null;
			if (seq != null)
				seqHashes = getSketch(seq);

			if (seqHashes == null)
				return false;
			processAddition(seqHashes);

			this.sequenceHashList.add(seqHashes);

			// fasta files are all fwd
			if (!fwdOnly)
			{
				// compute the hashes
				seqHashes = getSketch(seq.getReverseCompliment());

				this.sequenceHashList.add(seqHashes);
				processAddition(seqHashes);
			}
		}
		else
		{
			// read the binary file
			seqHashes = readFromBinary(buf, fwdOnly);
			while (seqHashes != null && fwdOnly && !seqHashes.getSequenceId().isForward())
			{
				seqHashes = readFromBinary(buf, fwdOnly);
			}

			// do nothing and return
			// record
			if (seqHashes == null)
				return false;

			processAddition(seqHashes);

			this.sequenceHashList.add(seqHashes);

		}

		return true;
	}

	public synchronized void enqueueFullFile(final boolean fwdOnly, int numThreads) throws IOException
	{
		// figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);

		// for each thread create a task
		for (int iter = 0; iter < numThreads; iter++)
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

	public Iterator<SequenceSketch> getDataIterator()
	{
		return this.sequenceHashList.iterator();
	}

	public int getFastaProcessed()
	{
		if (this.fastaData == null)
			return 0;

		return this.fastaData.getNumberProcessed();
	}

	public SequenceSketch getSketch(Sequence seq)
	{
		// compute the hashes
		return new SequenceSketch(seq, this.kmerSize, this.numHashes, this.orderedKmerSize, false, this.filter, this.kmerCounter, this.weighted);
	}

	public int getNumberProcessed()
	{
		return this.numberProcessed.intValue();
	}

	public int getNumberSubSequencesProcessed()
	{
		return this.numberSubSequencesProcessed.intValue();
	}

	protected void processAddition(SequenceSketch seqHashes)
	{
		// increment counter
		this.numberProcessed.getAndIncrement();

		int numProcessed = getNumberProcessed();
		if (numProcessed % 5000 == 0)
			System.err.println("Current # sequences loaded and processed from file: " + numProcessed + "...");

		if (seqHashes != null)
			this.numberSubSequencesProcessed.getAndAdd(1);
	}

	protected SequenceSketch readFromBinary(ReadBuffer buf, boolean fwdOnly) throws IOException
	{
		byte[] byteArray = null;
		synchronized (this.buffInput)
		{
			if (this.readClosed)
				return null;

			try
			{
				boolean keepReading = true;
				while (keepReading)
				{
					byte isFwd = this.buffInput.readByte();

					if (!fwdOnly || isFwd == 1)
						keepReading = false;

					// get the size in bytes
					int byteSize = this.buffInput.readInt();

					// allocate the array
					byteArray = buf.getBuffer(byteSize);
					// byteArray = new byte[byteSize];

					// read that many bytes
					this.buffInput.read(byteArray, 0, byteSize);
				}
			}
			catch (EOFException e)
			{
				this.buffInput.close();
				this.readClosed = true;

				return null;
			}
		}

		// get as byte array stream
		SequenceSketch seqHashes = SequenceSketch.fromByteStream(new DataInputStream(
				new ByteArrayInputStream(byteArray)), this.offset);

		return seqHashes;
	}

	public void writeToBinary(String file, final boolean fwdOnly, int numThreads) throws IOException
	{
		OutputStream output = null;
		try
		{
			output = new BufferedOutputStream(new FileOutputStream(file), Utils.BUFFER_BYTE_SIZE);
			final OutputStream finalOutput = output;

			// figure out number of cores
			ExecutorService execSvc = Executors.newFixedThreadPool(numThreads);

			// for each thread create a task
			for (int iter = 0; iter < numThreads; iter++)
			{
				Runnable task = new Runnable()
				{
					@Override
					public void run()
					{
						SequenceSketch seqHashes;
						ReadBuffer buf = new ReadBuffer();

						try
						{
							seqHashes = dequeue(fwdOnly, buf);
							while (seqHashes != null)
							{
								byte[] byteArray = seqHashes.getAsByteArray();
								int arraySize = byteArray.length;
								byte isFwd = seqHashes.getSequenceId().isForward() ? (byte) 1 : (byte) 0;

								// store the size as byte array
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

			finalOutput.flush();
		}
		finally
		{
			if (output != null)
				output.close();
		}
	}
}
