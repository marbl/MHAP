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
package edu.umd.marbl.mhap.utils;

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