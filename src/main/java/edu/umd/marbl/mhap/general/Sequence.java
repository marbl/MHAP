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

import edu.umd.marbl.mhap.utils.Utils;

public final class Sequence
{
	private final String sequence;
	private final SequenceId id;
	
	public Sequence(int[] sequence, SequenceId id)
	{
		this.id = id;
		
		StringBuilder s = new StringBuilder();
		for (int iter=0; iter<sequence.length; iter++)
		{
			switch(sequence[iter])
			{
				case 0 : s.append("U"); break;
				case 1 : s.append("C"); break;
				case 2 : s.append("G"); break;
				case 3 : s.append("T"); break;
				default : throw new RuntimeException("Uknown integer value.");
			}
		}
		
		this.sequence = s.toString();
	}
	
	public Sequence(String sequence, SequenceId id)
	{
		this.sequence = sequence;
		this.id = id;
	}
	
	public String getString()
	{
		return this.sequence;
	}
	
	public SequenceId getId()
	{
		return this.id;
	}
	
	public Sequence getReverseCompliment()
	{
		return new Sequence(Utils.rc(this.sequence), this.id.complimentId());
	}
	
	public String getKmer(int index, int kmerSize)
	{
		return this.sequence.substring(index, index+kmerSize);
	}
	
	public int numKmers(int kmerSize)
	{
		return this.sequence.length()-kmerSize+1;
	}
	
	public int length()
	{
		return this.sequence.length();
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		StringBuilder str = new StringBuilder();
		
		str.append(">"+this.id+"\n");
		str.append(this.sequence);
		
		return str.toString();
	}
}
