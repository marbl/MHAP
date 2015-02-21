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


public final class MatchResult implements Comparable<MatchResult>
{
	private final SequenceId fromId;
	private final SequenceId toId;
	private final int a1;
	private final int a2;
	private final int b1;
	private final int b2;
	private final double score;
	private final double rawScore;
	private final int fromLength;
	private final int toLength;
	
	protected MatchResult(SequenceId fromId, SequenceId toId, OverlapInfo overlap, int fromLength, int toLength)
	{
		this.fromId = fromId;
		this.toId = toId;
		
		this.fromLength = fromLength;
		this.toLength = toLength;
		
		this.a1 = getFromId().isForward() ? overlap.a1 : fromLength-overlap.a2-1;
		this.a2 = getFromId().isForward() ? overlap.a2 : fromLength-overlap.a1-1;
		this.b1 = getToId().isForward() ? overlap.b1 : toLength-overlap.b2-1;
		this.b2 = getToId().isForward() ? overlap.b2 : toLength-overlap.b1-1;
		
		this.rawScore = overlap.rawScore;
		
		if (overlap.score>1.0)
			this.score = 1.0;
		else
			this.score = overlap.score;
	}

	/**
	 * @return the fromId
	 */
	public SequenceId getFromId()
	{
		return this.fromId;
	}

	/**
	 * @return the toId
	 */
	public SequenceId getToId()
	{
		return this.toId;
	}

	/**
	 * @return the score
	 */
	public double getScore()
	{
		return this.score;
	}

	@Override
	public int compareTo(MatchResult o)
	{
		return -Double.compare(this.score, o.score);
	}
	
	@Override
	public String toString()
	{
		return String.format("%s %s %.6f %.6f %d %d %d %d %d %d %d %d",
				getFromId().getHeader(), 
				getToId().getHeader(),
				(1.0-getScore())*100.0,
				this.rawScore,
				getFromId().isForward() ? 0 : 1,
				this.a1,
				this.a2,
				this.fromLength,
				getToId().isForward() ? 0 : 1,
				this.b1,
				this.b2,
				this.toLength);
	}


}
