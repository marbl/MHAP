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
 * Copyright (c) 2015 by Konstantin Berlin and Sergey Koren
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
package edu.umd.marbl.mhap.sketch;

import java.util.concurrent.atomic.LongAdder;

public final class CountMin<T extends Object> implements Counter<T>
{
	private final LongAdder[][] countTable;
	private final int depth;
	private final int seed;

	private final LongAdder totalAdded;
	private final int width;

	public CountMin(double eps, double confidence, int seed)
	{
		// 2/w = eps ; w = 2/eps
		// 1/2^depth <= 1-confidence ; depth >= -log2 (1-confidence)

		// estimate the table size
		// this.width = (int) Math.ceil((double)2 / eps);
		// this.depth = (int) Math.ceil(-Math.log(1.0 - confidence) /
		// Math.log(2));
		// this.seed = seed;

		this((int) Math.ceil(-Math.log(1.0 - confidence) / Math.log(2)), (int) Math.ceil((double) 2 / eps), seed);
	}

	public CountMin(int depth, int width, int seed)
	{
		this.depth = depth;
		this.width = width;
		this.seed = seed;

		this.countTable = new LongAdder[depth][width];
		this.totalAdded = new LongAdder();

		// zero all the elements
		for (int iter1 = 0; iter1 < depth; iter1++)
			for (int iter2 = 0; iter2 < width; iter2++)
				this.countTable[iter1][iter2] = new LongAdder();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.invincea.labs.pace.hash.Counter#add(java.lang.Object)
	 */
	@Override
	public void add(T obj)
	{
		add(obj, 1);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.invincea.labs.pace.hash.Counter#add(java.lang.Object, int)
	 */
	@Override
	public void add(T obj, long increment)
	{
		if (increment <= 0)
			throw new SketchRuntimeException("Positive value expected for increment.");

		// compute the hash
		int[] hashes = HashUtils.computeHashesInt(obj, depth, seed);

		for (int iter = 0; iter < depth; iter++)
		{
			this.countTable[iter][((hashes[iter] << 1) >>> 1) % width].add(increment);
		}
		
		//store the total
		this.totalAdded.add(increment);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.invincea.labs.pace.hash.Counter#getCount(java.lang.Object)
	 */
	@Override
	public long getCount(Object obj)
	{
		// compute the hash
		int[] hashes = HashUtils.computeHashesInt(obj, depth, seed);

		long mincount = Long.MAX_VALUE;

		for (int iter = 0; iter < depth; iter++)
		{
			long value = this.countTable[iter][((hashes[iter] << 1) >>> 1) % width].longValue();
			if (mincount > value)
				mincount = value;
		}

		return mincount;
	}

	public int getDepth()
	{
		return this.depth;
	}

	public int getWidth()
	{
		return this.width;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see com.invincea.labs.pace.hash.Counter#maxCount()
	 */
	@Override
	public long maxCount()
	{
		throw new SketchRuntimeException("Method not implemented.");
	}

	public long totalAdded()
	{
		return this.totalAdded.longValue();
	}

}
