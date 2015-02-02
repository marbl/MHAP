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
package edu.umd.marbl.mhap.impl;

import java.io.Serializable;

public final class SequenceId implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2181572437818064822L;
	private final int id;
	private final boolean isFwd;
	private final String strId;
	
	public static boolean STORE_FULL_ID = false; 
	
	public SequenceId(int id)
	{
		this(id, true);
	}
	
	public SequenceId(int id, boolean isFwd)
	{
		this.id = id;
		this.isFwd = isFwd;
		this.strId = null;
	}
	
	public SequenceId(int id, boolean isFwd, String strId)
	{
		this.id = id;
		this.isFwd = isFwd;
		this.strId = strId;
	}
	
	public SequenceId createOffset(int offset)
	{
		return new SequenceId(this.id+offset, this.isFwd, this.strId);
	}
	
	public SequenceId complimentId()
	{
		return new SequenceId(this.id, !this.isFwd, this.strId);
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
		SequenceId other = (SequenceId) obj;
		
		return (this.id==other.id) && (this.isFwd == other.isFwd);
	}
	
	public boolean isForward()
	{
		return this.isFwd;
	}
	
	public int getHeaderId()
	{
		return this.id;
	}

	public String getHeader()
	{
		if (this.strId!=null)
			return this.strId;
		
		return String.valueOf(this.id);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode()
	{
		return this.isFwd? this.id : -this.id;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return ""+getHeader()+(this.isFwd ? "(fwd)" : "(rev)");
	}
	
	public String toStringInt()
	{
		return ""+getHeader()+(this.isFwd ? " 1" : " 0");
	}
}
