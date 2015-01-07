package edu.umd.marbl.mhap.utils;

public final class HitCounter
{
	public int count;

	public HitCounter()
	{
		this.count = 0;
	}
	
	public HitCounter(int count)
	{
		this.count = count;
	}

	public void addHit()
	{
		this.count++;
	}
	
	public void addHits(int counts)
	{
		this.count+=counts;
	}
}