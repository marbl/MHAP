package edu.umd.marbl.mhap.sketch;

public interface Counter<T extends Object> extends Filter
{
	public long getCount(T obj);

	public void add(T obj);

	public long maxCount();

	public void add(T obj, long count);

}