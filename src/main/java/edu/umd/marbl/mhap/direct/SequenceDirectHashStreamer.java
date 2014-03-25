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
package edu.umd.marbl.mhap.direct;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;

import edu.umd.marbl.mhap.general.AbstractSequenceHashStreamer;
import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.ReadBuffer;
import edu.umd.marbl.mhap.utils.Utils;

public class SequenceDirectHashStreamer extends AbstractSequenceHashStreamer<SequenceDirectHashes>
{
	private final DataInputStream buffInput;
	private final HashSet<Integer> filter;
	private final int kmerSize;
	private final int offset;
	private boolean readClosed;
	private final int orderedKmerSize;

	public SequenceDirectHashStreamer(String file, int offset) throws FileNotFoundException
	{
		super(null, false);
		
		this.kmerSize = 0;
		this.orderedKmerSize = 0;
		this.filter = null;
		this.readClosed = false;
		this.offset = offset;

		this.buffInput = new DataInputStream(new BufferedInputStream(new FileInputStream(file), Utils.BUFFER_BYTE_SIZE));  	
	}
	
	public SequenceDirectHashStreamer(String file, int kmerSize, int orderedKmerSize,
			HashSet<Integer> filter, int offset) throws IOException
	{	
		super(new FastaData(file, offset), true);

		this.kmerSize = kmerSize;
		this.orderedKmerSize = orderedKmerSize;
		this.filter = filter;
		this.buffInput = null;
		this.readClosed = false;
		this.offset = offset;
	}

	@Override
	public SequenceDirectHashes getHashes(Sequence seq)
	{
		//compute the hashes
		return new SequenceDirectHashes(seq, this.kmerSize, this.orderedKmerSize, this.filter);
	}
	
	@Override
	protected SequenceDirectHashes readFromBinary(ReadBuffer buf, boolean fwdOnly) throws IOException
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
		SequenceDirectHashes seqHashes = SequenceDirectHashes.fromByteStream(new DataInputStream(new ByteArrayInputStream(byteArray)), this.offset);

  	return seqHashes;
  }
}
