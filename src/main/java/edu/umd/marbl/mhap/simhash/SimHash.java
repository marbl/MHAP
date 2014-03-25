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
package edu.umd.marbl.mhap.simhash;

import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.Utils;

public final class SimHash
{
	private long[] bits;
	
	public SimHash(Sequence seq, int kmerSize, int numberWords)
	{
		//compute the hashes
		long[][] hashes = Utils.computeKmerHashes(seq, kmerSize, numberWords);
		
		recordHashes(hashes, kmerSize, numberWords);
	}
	
	private final void recordHashes(final long[][] hashes, final int kmerSize, final int numWords)
	{
		final int[] counts = new int[numWords*64];

		//perform count for each kmer
		for (int kmerIndex=0; kmerIndex<hashes.length; kmerIndex++)
		{
			long[] kmerHashes = hashes[kmerIndex];
			
		  for (int longIndex=0; longIndex<numWords; longIndex++)
		  {	      
		  	final long val = kmerHashes[longIndex];
		  	final int offset = longIndex*64;

		  	long mask = 0b1;
		  	
	      for (int bit=0; bit<64; bit++)
	      {
	        /* if not different then increase count */
	        if ((val&mask)==0b0)
	          counts[offset+bit]--;
	        else
	        	counts[offset+bit]++;
	        	
	        mask = mask << 1;
	      }
		  }		  
		}
		
		this.bits = new long[numWords];
	  for (int longIndex=0; longIndex<numWords; longIndex++)
	  {	      
	  	final int offset = longIndex*64;
	  	long val = 0b0;
	  	long mask = 0b1;
	  	
      for (int bit=0; bit<64; bit++)
      {
      	if (counts[offset+bit]>0)
      		val = val | mask;
      	
      	//adjust the mask
        mask = mask << 1;
      }
      
      this.bits[longIndex] = val;
	  }	  
	}
}
