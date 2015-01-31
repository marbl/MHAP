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

import edu.umd.marbl.mhap.utils.MhapRuntimeException;

public abstract class AbstractBitSketch<T extends AbstractBitSketch<T>> implements VectorHash<T>,
		Comparable<T>
{
	protected final long[] bits;

	protected AbstractBitSketch(long[] bits)
	{
		this.bits = bits;
	}

	public final double adjScore(T sh)
	{
		double score = jaccard(sh);

		return score;
	}

	public long[] getBits()
	{
		return this.bits;
	}

	@Override
	public int compareTo(final T sim)
	{
		for (int bitIndex = 0; bitIndex < this.bits.length; bitIndex++)
		{
			if (this.bits[bitIndex] < sim.bits[bitIndex])
				return -1;
			if (this.bits[bitIndex] > sim.bits[bitIndex])
				return 1;
		}

		return 0;
	}

	public final int getIntersectionCount(final T sh)
	{
		if (this.bits.length != sh.bits.length)
			throw new MhapRuntimeException("Size of bits in tables must match.");

		int count = 0;
		for (int longIndex = 0; longIndex < this.bits.length; longIndex++)
		{
			final long xor = this.bits[longIndex] ^ sh.bits[longIndex];

			count += Long.bitCount(xor);
		}

		return this.bits.length * 64 - count;
	}

	@Override
	public final double jaccard(final T sh)
	{
		int count = getIntersectionCount(sh);
		
		double jaccard = ((double) count / (double) (this.bits.length * 64) - 0.5) * 2.0;

		return Math.max(0.0, jaccard);
	}

	@Override
	public String toString()
	{
		StringBuilder s = new StringBuilder();
		for (int longIndex = 0; longIndex < this.bits.length; longIndex++)
		{
			long mask = 1L << 63;

			for (int bit = 63; bit >= 0; bit--)
			{
				if ((this.bits[longIndex] & mask) == 0)
					s.append("0");
				else
					s.append("1");

				mask = mask >>> 1;
			}
		}
		return s.toString();
	}
}