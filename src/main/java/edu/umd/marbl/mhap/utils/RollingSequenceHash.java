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
package edu.umd.marbl.mhap.utils;

public final class RollingSequenceHash
{
	private final int kmerSize;

	public RollingSequenceHash(int kmerSize)
	{
		this.kmerSize = kmerSize;
	}
	
	public final int[] hashInt(final String seq)
	{
		//get the string
		final char[] seqArray = seq.toCharArray();
		
		//RabinKarpHash ch = new RabinKarpHash(this.kmerSize);
		CyclicHash ch = new CyclicHash(this.kmerSize);
		
		//allocate hashes
		int[] hashes = new int[seqArray.length-this.kmerSize+1];
		
		int kmerIndex;
		for(kmerIndex = 0; kmerIndex<this.kmerSize-1; kmerIndex++) 
		{
			ch.eat(seqArray[kmerIndex]);
		}
		
		int count = 0;
		int hash = ch.eat(seqArray[kmerIndex++]);		
		hashes[count++] = hash;
		
		//start rolling
		for(; kmerIndex<seqArray.length; kmerIndex++) 
		{
			//System.out.println(k-this.kmerSize);
			hash = ch.update(seqArray[kmerIndex-this.kmerSize], seqArray[kmerIndex]);			
			hashes[count++] = hash;
		}
		
		return hashes;
	}
}
