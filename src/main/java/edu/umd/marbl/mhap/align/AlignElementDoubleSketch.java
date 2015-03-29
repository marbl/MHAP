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

import edu.umd.marbl.mhap.impl.OverlapInfo;
import edu.umd.marbl.mhap.sketch.Sketch;

public final class AlignElementDoubleSketch<T extends Sketch<T>> implements AlignElement<AlignElementDoubleSketch<T>>
{
	private final T[] elements;
	private final int seqLength;
	private final int stepSize;
	
	public AlignElementDoubleSketch(T[] sketchArray, int stepSize, int seqLength)
	{
		this.elements = sketchArray;
		this.stepSize = stepSize;
		this.seqLength = seqLength;
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementDoubleSketch<T>> aligner, AlignElementDoubleSketch<T> b)
	{
		Alignment<AlignElementDoubleSketch<T>> alignment = localAlignOneSkip(aligner, b);
		
		int a1 = alignment.getA1()*2;
		int a2 = alignment.getA2()*2;
		int b1 = alignment.getB1()*2;
		int b2 = alignment.getB2()*2;
		
		if (alignment.getScore()<0.0)
			return new OverlapInfo(0.0, 0.0, 0, 0, 0, 0);
		
		int offsetStart = similarityOffset(b, alignment.getA1(), alignment.getB1());
		int offsetEnd = similarityOffset(b, alignment.getA2(), alignment.getB2());
		
		if (offsetStart>0)
			a1++;
		else
		if (offsetStart<0)
			b1++;
		if (offsetEnd>0)
			a2++;
		else
		if (offsetEnd<0)
			b2++;
		
		a1 = a1*this.stepSize;
		a2 = Math.min(getSequenceLength()-1, (a2*this.stepSize+this.stepSize-1));

		b1 = b1*b.stepSize;		
		b2 = Math.min(b.getSequenceLength()-1, (b2*b.stepSize+b.stepSize-1));
		
		double score = alignment.getScore();
		
		//int overlapSize = Math.max(a2-a1, b2-b1);
		//double relOverlapSize = (double)overlapSize/(double)this.stepSize;
		//score = score/relOverlapSize;
		
		return new OverlapInfo(score/100000.0, score, a1, a2, b1, b2);
	}

	public int getSequenceLength()
	{
		return this.seqLength;
	}

	public T getSketch(int index)
	{
		return this.elements[index];
	}
	
	public int getStepSize()
	{
		return this.stepSize;
	}
	
	@Override
	public int length()
	{
		int val = this.elements.length/2;
		if (this.elements.length%2!=0)
			val++;
		
		return val;
	}
	
	public Alignment<AlignElementDoubleSketch<T>> localAlignOneSkip(Aligner<AlignElementDoubleSketch<T>> aligner, AlignElementDoubleSketch<T> b)
	{
		return aligner.localAlignOneSkip(this, b);
	}
	
	@Override
	public double similarityScore(AlignElementDoubleSketch<T> e, int i, int j)
	{
		double max = this.elements[2*i].similarity(e.elements[2*j]);
		
		if ((2*i+1)<this.elements.length)
			max = Math.max(max, this.elements[2*i+1].similarity(e.elements[2*j]));
		if ((2*j+1)<e.elements.length)
			max = Math.max(max, this.elements[2*i].similarity(e.elements[2*j+1]));
		
		return max;
	}
	
	private int similarityOffset(AlignElementDoubleSketch<T> e, int i, int j)
	{
		double max = this.elements[2*i].similarity(e.elements[2*j]);
		int diff = 0;
		
		if ((2*i+1)<this.elements.length)
		{
			double val = this.elements[2*i+1].similarity(e.elements[2*j]);
			if (max<val)
			{
				max = val;
				diff = 1;
			}
		}
		if ((2*j+1)<e.elements.length)
		{
			double val = this.elements[2*i].similarity(e.elements[2*j+1]);
			if (max<val)
			{
				max = val;
				diff = -1;
			}
		}
		
		return diff;
	}
	
	@Override
	public String toString(AlignElementDoubleSketch<T> match, int i, int j)
	{
		return toString();
	}
	
	@Override
	public String toString(int i)
	{
		return this.elements[i].toString();
	}
}
