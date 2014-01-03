package com.secret.fastalign.utils;
public class SortablePair<ComparableType extends Comparable<ComparableType>,AnyType> extends Pair<ComparableType,AnyType> implements Comparable<SortablePair<ComparableType,AnyType>>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2817516347839329908L;

	/**
	 * Instantiates a new sortable pair.
	 * 
	 * @param x
	 *          the x
	 * @param y
	 *          the y
	 */
	public SortablePair(ComparableType x, AnyType y)
	{
		super(x,y);
	}
		
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public final int compareTo(SortablePair<ComparableType, AnyType> pair)
	{
		return this.x.compareTo(pair.x);
	}
}
