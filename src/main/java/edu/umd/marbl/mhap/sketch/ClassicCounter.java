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
