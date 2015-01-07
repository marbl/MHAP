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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Locale;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;

import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.Utils;

public class FastaData implements Cloneable
{
	private final BufferedReader fileReader;
	private final int offset;
	private String lastLine;
	private AtomicLong numberProcessed;
	private boolean readFullFile;
	// length of sequences loaded
	private final ConcurrentLinkedQueue<Sequence> sequenceList;

	private static final String[] fastaSuffix = { "fna", "contigs", "contig", "final", "fasta", "fa" };

	private FastaData(ConcurrentLinkedQueue<Sequence> seqList)
	{
		this.sequenceList = new ConcurrentLinkedQueue<Sequence>(seqList);
		this.fileReader = null;
		this.lastLine = null;
		this.readFullFile = true;
		this.numberProcessed = new AtomicLong(this.sequenceList.size());
		this.offset = 0;
	}

	public FastaData(String file, int offset) throws IOException
	{
		try
		{
			this.fileReader = Utils.getFile(file, fastaSuffix);
		}
		catch (Exception e)
		{
			throw new MhapRuntimeException(e);
		}

		this.offset = offset;
		this.lastLine = null;
		this.readFullFile = false;
		this.numberProcessed = new AtomicLong(0);
		this.sequenceList = new ConcurrentLinkedQueue<Sequence>();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public synchronized FastaData clone()
	{
		// enqueue all the data
		try
		{
			enqueueFullFile();
		}
		catch (IOException e)
		{
			throw new MhapRuntimeException(e);
		}

		return new FastaData(this.sequenceList);
	}

	public Sequence dequeue() throws IOException
	{
		Sequence seq;
		synchronized (this.sequenceList)
		{
			if (this.sequenceList.isEmpty())
			{
				enqueueNextSequenceInFile();
			}
	
			// get the sequence
			seq = this.sequenceList.poll();		
		}

		return seq;
	}

	public void enqueueFullFile() throws IOException
	{
		while (enqueueNextSequenceInFile())
		{
		}
	}

	private boolean enqueueNextSequenceInFile() throws IOException
	{
		synchronized (this.fileReader)
		{
			if (this.readFullFile)
				return false;

			// try to read the next line
			if (this.lastLine == null)
			{
				this.lastLine = this.fileReader.readLine();

				// there is no next line
				if (this.lastLine == null)
				{
					this.fileReader.close();
					this.readFullFile = true;
					return false;
				}
			}

			// process the header
			if (!this.lastLine.startsWith(">"))
				throw new MhapRuntimeException("Next sequence does not start with >. Invalid format.");

			// process the current header
			String header = null;
			if (SequenceId.STORE_FULL_ID)
				header = this.lastLine.substring(1).split("[\\s,]+", 2)[0];
			
			//read the first line of the sequence
			this.lastLine = this.fileReader.readLine();

			StringBuilder fastaSeq = new StringBuilder();
			while (true)
			{
				if (this.lastLine == null || this.lastLine.startsWith(">"))
				{
					//generate sequence id
					SequenceId id;
					if (SequenceId.STORE_FULL_ID)
						id = new SequenceId(this.numberProcessed.intValue() + this.offset + 1, true, header);
					else
						id = new SequenceId(this.numberProcessed.intValue() + this.offset + 1);

					Sequence seq = new Sequence(fastaSeq.toString().toUpperCase(Locale.ENGLISH), id);

					// enqueue sequence
					this.sequenceList.add(seq);
					this.numberProcessed.getAndIncrement();

					if (this.lastLine == null)
					{
						this.fileReader.close();
						this.readFullFile = true;
					}

					return true;
				}

				// append the last line
				fastaSeq.append(this.lastLine);
				this.lastLine = this.fileReader.readLine();
			}
		}

	}

	public int getNumberProcessed()
	{
		return this.numberProcessed.intValue();
	}

	public Sequence getSequence(SequenceId id)
	{
		if (id.isForward())
		{
			for (Sequence seq : this.sequenceList)
				if (seq.getId().equals(id))
					return seq;
		}

		id = id.complimentId();
		for (Sequence seq : this.sequenceList)
			if (seq.getId().equals(id))
				return seq.getReverseCompliment();

		return null;
	}

	public boolean isEmpty()
	{
		return this.sequenceList.isEmpty() && this.readFullFile;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable
	{
		super.finalize();
		this.fileReader.close();
	}
}
