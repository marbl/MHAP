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

public final class AlignElementSketch<T extends Sketch<T>> implements AlignElement<AlignElementSketch<T>>
{
	private final T[] elements;
	private final int seqLength;
	private final int stepSize;
	
	public AlignElementSketch(T[] sketchArray, int stepSize, int seqLength)
	{
		this.elements = sketchArray;
		this.stepSize = stepSize;
		this.seqLength = seqLength;
	}
	
	public OverlapInfo getOverlapInfo(Aligner<AlignElementSketch<T>> aligner, AlignElementSketch<T> b)
	{
		Alignment<AlignElementSketch<T>> alignment = localAlignOneSkip(aligner, b);
		
		int a1 = alignment.getA1();
		int a2 = alignment.getA2();
		int b1 = alignment.getB1();
		int b2 = alignment.getB2();
		
		a1 = alignment.getA1()*this.stepSize;		
		a2 = Math.min(getSequenceLength()-1, alignment.getA2()*this.stepSize+this.stepSize-1);
			
		b1 = alignment.getB1()*b.stepSize;		
		b2 = Math.min(b.getSequenceLength()-1, alignment.getB2()*b.stepSize+b.stepSize-1);
		
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
		return this.elements.length;
	}
	
	public Alignment<AlignElementSketch<T>> localAlignOneSkip(Aligner<AlignElementSketch<T>> aligner, AlignElementSketch<T> b)
	{
		return aligner.localAlignOneSkip(this, b);
	}
	
	@Override
	public double similarityScore(AlignElementSketch<T> e, int i, int j)
	{
		return this.elements[i].similarity(e.elements[j]);
	}
	
	@Override
	public String toString(AlignElementSketch<T> match, int i, int j)
	{
		return toString();
	}
	
	@Override
	public String toString(int i)
	{
		return this.elements[i].toString();
	}
}
