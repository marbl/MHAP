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

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashSet;

import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.FastAlignRuntimeException;
import edu.umd.marbl.mhap.utils.Utils;

public final class MinHash implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 8846482698636860862L;
	private final int[] minHashes;
	private final int seqLength;
	
	public static MinHash fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			//store the size
			//bb.putInt(this.seqLength);
			int seqLength = input.readInt();
			
			//bb.putInt(this.minHashes.length);
			int hashNum = input.readInt();
			
			//store the array
			int[] minHashes = new int[hashNum];
			for (int hash=0; hash<hashNum; hash++)
			{
				//bb.putInt(this.minHashes[seq][hash]);
				minHashes[hash] = input.readInt();
			}
			
			return new MinHash(seqLength, minHashes);
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	private MinHash(int seqLength, int[] minHashes)
	{
		this.seqLength = seqLength;
		this.minHashes = minHashes;
	}
	
	public MinHash(Sequence seq, int kmerSize, int numHashes, HashSet<Integer> filter)
	{
		this.seqLength = seq.length();

		this.minHashes = Utils.computeKmerMinHashes(seq.getString(), kmerSize, numHashes, filter);
	}

	public byte[] getAsByteArray()
	{
		ByteBuffer bb = ByteBuffer.allocate(4*(2+this.minHashes.length));
		
		//store the size
		bb.putInt(this.seqLength);
		bb.putInt(this.minHashes.length);
		
		//store the array
		for (int hash=0; hash<this.minHashes.length; hash++)
			bb.putInt(this.minHashes[hash]); 
    
    return bb.array();
	}
	
	public final int getSequenceLength()
	{
		return this.seqLength;
	}
	
	public final double jaccard(MinHash h)
	{
		int count = 0;
		int size = this.minHashes.length;
		
		if (h.minHashes.length!=size)
			throw new FastAlignRuntimeException("MinHashes must be of same length in order to be comapred.");
		
		for (int iter=0; iter<size; iter++)
		{
			if (this.minHashes[iter]==h.minHashes[iter])
				count++;
		}
		
		return (double)count/(double)size;
	}

	/**
	 * @return the minHashes
	 */
	public final int[] getMinHashArray()
	{
		return this.minHashes;
	}
	
	public final int numHashes()
	{
		return this.minHashes.length;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "MinHash "+Arrays.toString(this.minHashes) + "";
	}
}
