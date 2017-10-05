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
package edu.umd.marbl.mhap.impl;

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.Serializable;

import edu.umd.marbl.mhap.sketch.FrequencyCounts;
import edu.umd.marbl.mhap.sketch.MinHashSketch;
import edu.umd.marbl.mhap.sketch.BottomOverlapSketch;
import edu.umd.marbl.mhap.sketch.ZeroNGramsFoundException;

public final class SequenceSketch implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3155689614837922443L;

	private final SequenceId id;
	private final MinHashSketch mainHashes;
	private final BottomOverlapSketch orderedHashes;
	//private final MinHashBitSequenceSubSketches alignmentSketches;
	private final int sequenceLength;

	public final static int BIT_SKETCH_SIZE = 20;
	public final static int SUBSEQUENCE_SIZE = 50;
	public final static int BIT_KMER_SIZE = 7;

	public static SequenceSketch fromByteStream(DataInputStream input, int offset) throws IOException
	{
		try
		{
			// input.

			// dos.writeBoolean(this.id.isForward());
			boolean isFwd = input.readBoolean();

			// dos.writeLong(this.id.getHeaderId());
			// dos.writeUTF(this.id.getHeader());
			
			SequenceId id = new SequenceId(input.readLong() + offset, isFwd, input.readUTF());
			
			//dos.writeInt(this.sequenceLength);
			int sequenceLength = input.readInt();

			// dos.write(this.mainHashes.getAsByteArray());
			MinHashSketch mainHashes = MinHashSketch.fromByteStream(input);

			if (mainHashes == null)
				throw new MhapRuntimeException("Unexpected data read error.");

			BottomOverlapSketch orderedHashes = null;
			orderedHashes = BottomOverlapSketch.fromByteStream(input);
			if (orderedHashes == null)
				throw new MhapRuntimeException("Unexpected data read error when reading ordered k-mers.");

			return new SequenceSketch(id, sequenceLength, mainHashes, orderedHashes);

		}
		catch (EOFException e)
		{
			return null;
		}
	}

	public SequenceSketch(SequenceId id, int sequenceLength, MinHashSketch mainHashes, BottomOverlapSketch orderedHashes)
	{
		this.sequenceLength = sequenceLength;
		this.id = id;
		this.mainHashes = mainHashes;
		this.orderedHashes = orderedHashes;
	}

	public SequenceSketch(Sequence seq, int kmerSize, int numHashes, int orderedKmerSize, int orderedSketchSize, FrequencyCounts kmerFilter, boolean doReverseCompliment, double repeatWeight) throws ZeroNGramsFoundException
	{
		this.sequenceLength = seq.length();
		this.id = seq.getId();
		
		//do not do reverse compliment for minhash, since unordered
		this.mainHashes = new MinHashSketch(seq.getSquenceString(), kmerSize, numHashes, kmerFilter, false, repeatWeight);
		
		//do not do reverse compliment
		this.orderedHashes = new BottomOverlapSketch(seq.getSquenceString(), orderedKmerSize, orderedSketchSize, false);
	}

	public SequenceSketch createOffset(int offset)
	{
		return new SequenceSketch(this.id.createOffset(offset), this.sequenceLength, this.mainHashes, this.orderedHashes);
	}

	public byte[] getAsByteArray()
	{
		byte[] mainHashesBytes = this.mainHashes.getAsByteArray();		
		byte[] orderedHashesBytes = this.orderedHashes.getAsByteArray();
		
		//get size		
		ByteArrayOutputStream bos = new ByteArrayOutputStream(mainHashesBytes.length+orderedHashesBytes.length);
		DataOutputStream dos = new DataOutputStream(bos);

		try
		{
			dos.writeBoolean(this.id.isForward());
			dos.writeLong(this.id.getHeaderId());
			dos.writeUTF(this.id.getHeader());
			dos.writeInt(this.sequenceLength);
			dos.write(mainHashesBytes);
			dos.write(orderedHashesBytes);

			dos.flush();
			return bos.toByteArray();
		}
		catch (IOException e)
		{
			throw new MhapRuntimeException("Unexpected IO error.");
		}
	}
	
	public MinHashSketch getMinHashes()
	{
		return this.mainHashes;
	}

	public BottomOverlapSketch getOrderedHashes()
	{
		return this.orderedHashes;
	}
	
	public SequenceId getSequenceId()
	{
		return this.id;
	}

	public int getSequenceLength()
	{
		return this.sequenceLength;
	}
}
