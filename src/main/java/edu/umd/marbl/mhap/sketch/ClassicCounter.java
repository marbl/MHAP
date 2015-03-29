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

import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.LongAdder;

public final class ClassicCounter<T extends Object> implements Counter<T>
{
	private final HashMap<Object, LongAdder> map;
	private final LongAdder numAdditions;
	private final AtomicLong maxCount;
	
	public ClassicCounter(int size)
	{
		this.map = new HashMap<>(size);
		this.maxCount = new AtomicLong();
		this.numAdditions = new LongAdder();
	}

	@Override
	public long getCount(Object obj)
	{
		LongAdder adder = map.get(obj);
		if (adder==null)
			return 0;
		
		return map.get(obj).longValue();
	}

	@Override
	public void add(Object obj)
	{
		add(obj, 1);
	}

	@Override
	public long maxCount()
	{
		return this.maxCount.longValue();
	}

	@Override
	public void add(Object obj, long count)
	{
		LongAdder adder = null;
		synchronized (this.map)
		{
			adder = this.map.get(obj);
			if (adder==null)
			{
				adder = new LongAdder();
				this.map.put(obj, adder);
			}
		}
		
		adder.add(count);
		
		// assumes value always increasing
		if (maxCount.longValue() < count)
		{
			synchronized (maxCount)
			{
				//TODO fix
				if (maxCount.longValue() < count)
					maxCount.set(count);
			}
		}
		
		this.numAdditions.add(count);
	}

}
