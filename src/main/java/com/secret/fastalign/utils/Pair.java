package com.secret.fastalign.utils;

import java.io.Serializable;

public class Pair<A,B> implements Serializable
{
	public final A x;
	
	public final B y;
	/**
	 * 
	 */
	private static final long serialVersionUID = -5782450990742961765L;

	public Pair(A x, B y)
	{
		this.x = x;
		this.y = y;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		
		Pair<?, ?> other = (Pair<?,?>) obj;
		
		if (this.x == null)
		{
			if (other.x != null)
				return false;
		}
		else if (!this.x.equals(other.x))
			return false;
		if (this.y == null)
		{
			if (other.y != null)
				return false;
		}
		else if (!this.y.equals(other.y))
			return false;
		
		return true;
	}

  @Override
  public int hashCode() {
      return 31 * hashcode(this.x) + hashcode(this.y);
  }

  // todo move this to a helper class.
  private static int hashcode(Object o) {
      return o == null ? 0 : o.hashCode();
  }

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "[x=" + this.x + ", y=" + this.y + "]";
	}
}