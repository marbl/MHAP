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
import java.util.HashSet;

import edu.umd.marbl.mhap.sketch.MinHashSketch;
import edu.umd.marbl.mhap.sketch.MinHashSketchSequence;
import edu.umd.marbl.mhap.sketch.NGramCounts;
import edu.umd.marbl.mhap.sketch.OrderedNGramHashes;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;

public final class SequenceSketch implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3155689614837922443L;

	private final SequenceId id;
	private final MinHashSketch mainHashes;
	private final OrderedNGramHashes orderedHashes;
	private final MinHashSketchSequence bitSequence;
	private final int sequenceLength;

	public final static double SHIFT_CONSENSUS_PERCENTAGE = 0.75;
	public final static int BIT_SKETCH_SIZE = 6;
	public final static int SUBSEQUENCE_SIZE = 200;
	public final static int BIT_KMER_SIZE = 7;
	public final static int MIN_OVERLAP_SIZE = 1000;


	public static SequenceSketch fromByteStream(DataInputStream input, int offset) throws IOException
	{
		try
		{
			// input.

			// dos.writeBoolean(this.id.isForward());
			boolean isFwd = input.readBoolean();

			// dos.writeInt(this.id.getHeaderId());
			SequenceId id = new SequenceId(input.readInt() + offset, isFwd);
			
			//dos.writeInt(this.sequenceLength);
			int sequenceLength = input.readInt();

			// dos.write(this.mainHashes.getAsByteArray());
			MinHashSketch mainHashes = MinHashSketch.fromByteStream(input);

			if (mainHashes == null)
				throw new MhapRuntimeException("Unexpected data read error.");

			// dos.write(this.orderedHashes.getAsByteArray());
			OrderedNGramHashes orderedHashes = OrderedNGramHashes.fromByteStream(input);
			if (orderedHashes == null)
				throw new MhapRuntimeException("Unexpected data read error.");

			return new SequenceSketch(id, sequenceLength, mainHashes, orderedHashes);

		}
		catch (EOFException e)
		{
			return null;
		}
	}

	public SequenceSketch(SequenceId id, int sequenceLength, MinHashSketch mainHashes, OrderedNGramHashes orderedHashes)
	{
		this.sequenceLength = sequenceLength;
		this.id = id;
		this.mainHashes = mainHashes;
		this.orderedHashes = orderedHashes;
		this.bitSequence = null;
	}

	public SequenceSketch(Sequence seq, int kmerSize, int numHashes, int orderedKmerSize, boolean storeHashes,
			HashSet<Long> filter, NGramCounts kmerCount, boolean weighted)
	{
		this.sequenceLength = seq.length();
		this.id = seq.getId();
		this.mainHashes = new MinHashSketch(seq.toString(), kmerSize, numHashes, filter, kmerCount, weighted);
		this.orderedHashes = new OrderedNGramHashes(seq.toString(), orderedKmerSize);
		
		//bit sequence
		this.bitSequence = new MinHashSketchSequence(seq.toString(), BIT_KMER_SIZE, SUBSEQUENCE_SIZE, MIN_OVERLAP_SIZE, -0.15, BIT_SKETCH_SIZE);
	}

	public SequenceSketch createOffset(int offset)
	{
		return new SequenceSketch(this.id.createOffset(offset), this.sequenceLength, this.mainHashes, this.orderedHashes);
	}

	public byte[] getAsByteArray()
	{
		ByteArrayOutputStream bos = new ByteArrayOutputStream(this.mainHashes.numHashes() * 4
				+ this.orderedHashes.size() * 2);
		DataOutputStream dos = new DataOutputStream(bos);

		try
		{
			dos.writeBoolean(this.id.isForward());
			dos.writeInt(this.id.getHeaderId());
			dos.writeInt(this.sequenceLength);
			dos.write(this.mainHashes.getAsByteArray());
			dos.write(this.orderedHashes.getAsByteArray());

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

	public OrderedNGramHashes getOrderedHashes()
	{
		return this.orderedHashes;
	}
	
	public MinHashSketchSequence getBitSequence()
	{
		return this.bitSequence;
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
