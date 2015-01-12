package edu.umd.marbl.mhap.sketch;

import java.util.concurrent.atomic.LongAdder;

import edu.umd.marbl.mhap.general.Counter;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.Utils;

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
			throw new MhapRuntimeException("Positive value expected for increment.");

		// compute the hash
		int[] hashes = Utils.computeHashesInt(obj, depth, seed);

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
		int[] hashes = Utils.computeHashesInt(obj, depth, seed);

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
		throw new MhapRuntimeException("Method not implemented.");
	}

	public long totalAdded()
	{
		return this.totalAdded.longValue();
	}

}
