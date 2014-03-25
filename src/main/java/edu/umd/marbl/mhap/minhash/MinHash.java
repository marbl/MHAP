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
import edu.umd.marbl.mhap.utils.Utils;

public final class MinHash implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 8846482698636860862L;
	private final int[][] minHashes;
	private final int seqLength;
	
	public static MinHash fromByteStream(DataInputStream input) throws IOException
	{
		try
		{
			//store the size
			//bb.putInt(this.seqLength);
			int seqLength = input.readInt();
			
			//bb.putInt(this.minHashes.length);
			int seqNum = input.readInt();
			
			//bb.putInt(this.minHashes[0].length);
			int hashNum = input.readInt();
			
			//store the array
			int[][] minHashes = new int[seqNum][];
			for (int seq=0; seq<seqNum; seq++)
			{
				minHashes[seq] = new int[hashNum];
				for (int hash=0; hash<hashNum; hash++)
				{
					//bb.putInt(this.minHashes[seq][hash]);
					minHashes[seq][hash] = input.readInt();
				}
			}
			
			return new MinHash(seqLength, minHashes);
		}
		catch (EOFException e)
		{
			return null;
		}
	}
	
	private MinHash(int seqLength, int[][] minHashes)
	{
		this.seqLength = seqLength;
		this.minHashes = minHashes;
	}
	
	public MinHash(Sequence seq, int kmerSize, int numHashes, int subSequenceSize, HashSet<Integer> filter)
	{
		this.seqLength = seq.length();

		int numberSubSeq = seq.length()/subSequenceSize+1;
		if (seq.length()%subSequenceSize<kmerSize)
			numberSubSeq--;
		
		//adjust for more equal distribution
		subSequenceSize = seq.length()/numberSubSeq+1;
		
		this.minHashes = new int[numberSubSeq][];
		for (int iter=0; iter<numberSubSeq; iter++)
		{
			String subString = seq.getString().substring(iter*subSequenceSize, Math.min(seq.length(), (iter+1)*subSequenceSize));

			//get the hashes
			if (subString.length()>=kmerSize)
				this.minHashes[iter] = Utils.computeKmerMinHashes(subString, kmerSize, numHashes, filter);
		}	
	}

	public byte[] getAsByteArray()
	{
		ByteBuffer bb = ByteBuffer.allocate(4*(3+this.minHashes.length*this.minHashes[0].length));
		
		//store the size
		bb.putInt(this.seqLength);
		bb.putInt(this.minHashes.length);
		bb.putInt(this.minHashes[0].length);
		
		//store the array
		for (int seq=0; seq<this.minHashes.length; seq++)
			for (int hash=0; hash<this.minHashes[0].length; hash++)
				bb.putInt(this.minHashes[seq][hash]); 
    
    return bb.array();
	}
	
	public final int getSequenceLength()
	{
		return this.seqLength;
	}

	/**
	 * @return the minHashes
	 */
	public final int[][] getSubSeqMinHashArray()
	{
		return this.minHashes;
	}
	
	public final int numSubSequences()
	{
		return this.minHashes.length;
	}
	
	public final int numHashes()
	{
		return this.minHashes[0].length;
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
