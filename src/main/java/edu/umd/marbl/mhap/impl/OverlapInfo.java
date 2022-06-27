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
package edu.umd.marbl.mhap.impl;

import java.lang.Math;

public final class OverlapInfo
{
	public final int a1;
	public final int a2;
	public final int b1;
	public final int b2;
	public final double rawScore;
	public final double score;
	
	public static OverlapInfo EMPTY = new OverlapInfo(0.0, 0.0, 0, 0, 0, 0);
	
	public OverlapInfo(double score, double rawScore, int a1, int a2, int b1, int b2)
	{
		this.score = score;
		this.rawScore = rawScore;
		this.a1 = a1;
		this.a2 = a2;
		this.b1 = b1;
		this.b2 = b2;
	}
	
	/**
	 * @return the score
	 */
	public double getScore()
	{
		return this.score;
	}

	/**
	 * @return the length of the overlap
	 */
	public double getLength()
	{
		return (Math.abs(this.a1 - this.a2) + Math.abs(this.b1 - this.b2)) / 2;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString()
	{
		return "[score="+this.score+", a1="+this.a1+" a2="+this.a2+", b1="+this.b1+" b2="+this.b2+"]";
	}

}
