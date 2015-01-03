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
package edu.umd.marbl.mhap.minhash;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLong;

import edu.umd.marbl.mhap.general.AbstractSequenceHashStreamer;
import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.ReadBuffer;
import edu.umd.marbl.mhap.utils.Utils;

public class SequenceMinHashStreamer extends AbstractSequenceHashStreamer<SequenceMinHashes>
{
	private final DataInputStream buffInput;
	private final HashSet<Integer> filter;
	private final int kmerSize;
	private final AtomicLong numberSubSequencesProcessed;
	private final int numHashes;
	private final int offset;
	private boolean readClosed;
	private final int orderedKmerSize;

	public SequenceMinHashStreamer(String file, int offset) throws FileNotFoundException
	{
		super(null, false);
		
		this.kmerSize = 0;
		this.numHashes = 0;
		this.orderedKmerSize = 0;
		this.filter = null;
		this.numberSubSequencesProcessed = new AtomicLong();
		this.readClosed = false;
		this.offset = offset;

		this.buffInput = new DataInputStream(new BufferedInputStream(new FileInputStream(file), Utils.BUFFER_BYTE_SIZE));  	
	}
	
	public SequenceMinHashStreamer(String file, int kmerSize, int numHashes, int subSequenceSize, int orderedKmerSize, 
			HashSet<Integer> filter, int offset) throws IOException
	{	
		super(new FastaData(file, offset), true);

		this.kmerSize = kmerSize;
		this.numHashes = numHashes;
		this.orderedKmerSize = orderedKmerSize;
		this.filter = filter;
		this.numberSubSequencesProcessed = new AtomicLong();
		this.buffInput = null;
		this.readClosed = false;
		this.offset = offset;
	}

	@Override
	public SequenceMinHashes getHashes(Sequence seq)
	{
		//compute the hashes
		return new SequenceMinHashes(seq, this.kmerSize, this.numHashes, this.orderedKmerSize, false, this.filter);
	}
	
	@Override
	public int getNumberSubSequencesProcessed()
	{
		return this.numberSubSequencesProcessed.intValue();
	}

	@Override
	protected void processAddition(SequenceMinHashes seqHashes)
	{
		super.processAddition(seqHashes);
		if (seqHashes!=null)
			this.numberSubSequencesProcessed.getAndAdd(1);
	}

	@Override
	protected SequenceMinHashes readFromBinary(ReadBuffer buf, boolean fwdOnly) throws IOException
	{
		byte[] byteArray = null;
		synchronized (this.buffInput)
		{
			if (this.readClosed)
				return null;
			
			try
			{
				boolean keepReading = true;
				while(keepReading)
				{
					byte isFwd = this.buffInput.readByte();
					
					if (!fwdOnly || isFwd==1)
						keepReading = false;
					
					//get the size in bytes
					int byteSize = this.buffInput.readInt();
					
					//allocate the array
					byteArray = buf.getBuffer(byteSize);
					//byteArray = new byte[byteSize];				
					
					//read that many bytes
					this.buffInput.read(byteArray, 0, byteSize);
				}
			}
			catch(EOFException e)
			{
	  		this.buffInput.close();
	  		this.readClosed = true;			
	  		
	  		return null;
			}
		}
			
		//get as byte array stream
  	SequenceMinHashes seqHashes = SequenceMinHashes.fromByteStream(new DataInputStream(new ByteArrayInputStream(byteArray)), this.offset);

  	return seqHashes;
  }
}
