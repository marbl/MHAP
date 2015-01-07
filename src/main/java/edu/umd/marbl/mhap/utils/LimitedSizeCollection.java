/* 
 *  
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2013 by Konstantin Berlin 
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

import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.PriorityQueue;

public final class LimitedSizeCollection<T extends Comparable<T>> implements Collection<T>
{
	public enum Priority
	{
		MAX_VALUES, MIN_VALUES;
	}

	private T best;
	private int maxSize;
	private final PriorityQueue<T> queue;

	public LimitedSizeCollection(int maxSize)
	{
		this(maxSize, Priority.MIN_VALUES);
	}

	public LimitedSizeCollection(int maxSize, Priority priority)
	{
		// initiate with reverse queue
		if (priority == Priority.MIN_VALUES)
		{
			this.queue = new PriorityQueue<T>(maxSize, new Comparator<T>()
			{
				@Override
				public final int compare(T s1, T s2)
				{
					return s2.compareTo(s1);
				}
			});
		}
		else
		{
			this.queue = new PriorityQueue<T>(maxSize, new Comparator<T>()
			{
				@Override
				public final int compare(T s1, T s2)
				{
					return s1.compareTo(s2);
				}
			});
		}

		this.maxSize = maxSize;
		this.best = null;
	}

	@Override
	public boolean add(T o)
	{
		if (o == null)
			return false;

		if (this.maxSize <= 0)
			return false;

		// if can fit just add
		if (this.queue.size() < this.maxSize)
		{
			this.queue.add(o);
		}
		else if (this.queue.comparator().compare(o, this.queue.peek()) > 0)
		{
			this.queue.add(o);
			this.queue.poll();
		}
		else
			return false;

		if (this.best == null || this.queue.comparator().compare(o, this.best) > 0)
		{
			this.best = o;
		}

		return true;
	}

	@Override
	public boolean addAll(Collection<? extends T> c)
	{
		for (T elem : c)
			add(elem);

		return true;
	}

	@Override
	public void clear()
	{
		this.best = null;

		this.queue.clear();
	}

	@Override
	public boolean contains(Object o)
	{
		return this.queue.contains(o);
	}

	@Override
	public boolean containsAll(Collection<?> c)
	{
		return this.queue.containsAll(c);
	}

	public T getBest()
	{
		return this.best;
	}

	public Collection<T> getCollection()
	{
		return this.queue;
	}

	public T getWorst()
	{
		return this.queue.peek();
	}

	@Override
	public boolean isEmpty()
	{
		return this.queue.isEmpty();
	}

	public boolean isFull()
	{
		return this.queue.size() >= this.maxSize;
	}

	@Override
	public Iterator<T> iterator()
	{
		return this.queue.iterator();
	}

	@Override
	public boolean remove(Object o)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean removeAll(Collection<?> c)
	{
		this.best = null;

		return this.queue.removeAll(c);
	}

	public T removeWorst()
	{
		return this.queue.poll();
	}

	@Override
	public boolean retainAll(Collection<?> c)
	{
		throw new UnsupportedOperationException();
	}

	public void setSize(int maxSize)
	{
		this.maxSize = maxSize;
		while (this.queue.size() > this.maxSize)
			this.queue.poll();
	}

	@Override
	public int size()
	{
		return this.queue.size();
	}

	@Override
	public Object[] toArray()
	{
		return this.queue.toArray();
	}

	@Override
	public <Y> Y[] toArray(Y[] a)
	{
		return this.queue.toArray(a);
	}
}