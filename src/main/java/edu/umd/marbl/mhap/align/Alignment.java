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

import java.util.Iterator;
import java.util.List;

public final class Alignment<S extends AlignElement<S>>
{
	public enum Operation
	{
		MATCH,
		INSERT,
		DELETE;
	}
	
	private final double score;
	//private final double gapOpen;
	private final List<Operation> operations;
	private final S a;
	private final S b;
	private int a1;
	private int a2;
	private int b1;
	private int b2;
	
	protected Alignment(S a, S b, int a1, int a2, int b1, int b2, double score, double gapOpen, List<Operation> operations)
	{
		this.score = score;
		this.operations = operations;
		this.a = a;
		this.b = b;
		this.a1 = a1;
		this.a2 = a2;
		this.b1 = b1;
		this.b2 = b2;
		//this.gapOpen = gapOpen;
	}
	
	public double getOverlapScore(int minMatches)
	{
		int i = 0;
		int j = 0;
		
		Iterator<Operation> iter = this.operations.iterator();
		
		if (!iter.hasNext())
			return 0.0;
		
		//remove the start
		Operation o = iter.next();			
		while (o==Operation.DELETE)
		{
			i++;
			
			if (iter.hasNext())
				o = iter.next();
			else
				return 0.0;
		}
		if (i==0)
		{
			while (o==Operation.INSERT)
			{
				if (iter.hasNext())
					o = iter.next();
				else
					return 0.0;
			}
		}
		
		double score = 0.0;
		int count = 0;
		while (true)
		{
			if (o == Operation.DELETE)
			{	
				i++;
			}
			else
			if (o == Operation.INSERT)
			{
				//count++;
				j++;
			}
			else
			{
				score = score + a.similarityScore(b, i, j);
				count++;
				i++;
				j++;
			}
			
			if (iter.hasNext())
				o = iter.next();
			else
				break;
		}
		
		//System.err.println(this.operations);
		//System.err.println("HI="+count+" "+minMatches);
		
		if (count<minMatches)
			return 0.0;
		
		if (score<=0.0)
			return 0.0;
		
		return score/(double)(count);
	}
	
	protected List<Operation> getOperations()
	{
		return this.operations;
	}
	
	public double getScore()
	{
		return this.score;
	}
	
	public int getA1()
	{				
		return this.a1;
	}
	
	public int getA2()
	{				
		return this.a2;
	}
	
	public int getB1()
	{
		return this.b1;
	}
	
	public int getB2()
	{
		return this.b2;
	}
	
	public String outputAlignmentSelf()
	{
		StringBuilder str = new StringBuilder();

		int i = 0;
		int j = 0;
		
		int count = 0;
		while(i<a.length() || j<b.length())
		{
			Operation o;
			if (count<this.operations.size())
				o = this.operations.get(count);
			else
			if (i<a.length())
				o = Operation.DELETE;
			else
				o = Operation.INSERT;
						
			if (o == Operation.DELETE)
			{	
				String aStr;
				if (j>=b.length())
					aStr = a.toString(i);
				else
					aStr = a.toString(b, i, j);
				
				str.append(aStr);
				i++;
			}
			else
			if (o == Operation.INSERT)
			{
				String bStr;
				if (i>=a.length())
					bStr = b.toString(j);
				else
					bStr = b.toString(a, j, i);
				for (int space=0; space<bStr.length(); space++)
					str.append("-");
				j++;
			}
			else
			{
				str.append(a.toString(b, i, j));
				i++;
				j++;
			}
			
			count++;
		}
		
		return str.toString();		
	}
	
	public String outputAlignmentOther()
	{
		StringBuilder str = new StringBuilder();

		int i = 0;
		int j = 0;
		int count = 0;
		while(i<a.length() || j<b.length())
		{
			Operation o;
			if (count<this.operations.size())
				o = this.operations.get(count);
			else
			if (i<a.length())
				o = Operation.DELETE;
			else
				o = Operation.INSERT;
			
			if (o == Operation.INSERT)
			{
				String bStr;
				if (i>=a.length())
					bStr = b.toString(j);
				else
					bStr = b.toString(a, j, i);

				str.append(bStr);
				j++;
			}
			else
			if (o == Operation.DELETE)
			{
				String aStr;
				if (j>=b.length())
					aStr = a.toString(i);
				else
					aStr = a.toString(b, i, j);
				
				for (int space=0; space<aStr.length(); space++)
					str.append("-");
				i++;
			}
			else
			{
				str.append(b.toString(a, j, i));
				i++;
				j++;
			}
			
			count++;
		}
		
		return str.toString();
	}
	
	public String outputAlignment()
	{
		StringBuilder str = new StringBuilder();
		
		str.append("Sequence 1: ");
		str.append(outputAlignmentSelf()+"\n");
		str.append("Sequence 2: ");
		str.append(outputAlignmentOther()+"\n");
		
		return str.toString();
	}
}
