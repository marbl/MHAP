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
package edu.umd.marbl.mhap.align;

public class AlignElementString implements AlignElement<AlignElementString>
{
	private final double EXACT_MATCH_SCORE = 1.0;
	
	private final double MISMATCH_SCORE = -1.0;
	private final char[] s;
	
	public AlignElementString(String s)
	{
		this.s = s.toCharArray();
	}
	
	@Override
	public int length()
	{
		return s.length;
	}

	@Override
	public double similarityScore(AlignElementString e, int i, int j)
	{
		if (this.s[i]==e.s[j])
			return EXACT_MATCH_SCORE;
		else
			return MISMATCH_SCORE;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return new String(s);
	}

	@Override
	public String toString(AlignElementString match, int i, int j)
	{
		return ""+s[i];
	}

	@Override
	public String toString(int i)
	{
		return ""+s[i];
	}
}
