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
package edu.umd.marbl.mhap.sketch;

import java.util.HashMap;
import java.util.Map;


public final class FrequencyCounts
{
	private final double filterCutoff;
	private final Map<Long,Double> fractionCounts;
	private final double maxIdfValue;
	private final double maxValue;
	private final double minIdfValue;
	private final double minValue;
	
	public FrequencyCounts(Map<Long,Double> fractionCounts, double filterCutoff)
	{
		this.fractionCounts = new HashMap<>(fractionCounts);
		this.filterCutoff = filterCutoff;
		
		double maxValue = Double.NEGATIVE_INFINITY;
		for (double val : this.fractionCounts.values())
			maxValue = Math.max(maxValue, val);
		
		this.maxValue = maxValue;
		this.minValue = this.filterCutoff;
		
		this.minIdfValue = idf(this.maxValue);
		this.maxIdfValue = idf(this.minValue);
	}
	
	public boolean contains(long hash)
	{
		return this.fractionCounts.containsKey(hash);
	}
	
	public double documentFrequencyRatio(long hash)
	{
		Double val = this.fractionCounts.get(hash);
		if (val == null)
			val = this.minValue;
		
		return val;
	}
	
	public double getFilterCutoff()
	{
		return this.filterCutoff;
	}
	
	public double idf(double freq)
	{
		return Math.log(this.maxValue/freq);
		//return Math.log1p(this.maxValue/freq);
	}
	
	public double idf(long hash)
	{
		double freq = documentFrequencyRatio(hash);
		return idf(freq); 
	}
	
	public double idfDiscrete(long hash, int maxValue)
	{
		Double val = this.fractionCounts.get(hash);
		if (val == null)
			return maxValue;
		
		//get the true value
		double idf = idf(val);
		
		//scale it to match max
		double scale = (maxIdf()-minIdf())/(double)(maxValue-1.0);
		
		return 1.0+(idf-minIdf())/scale;
	}
	
	public double inverseDocumentFrequency(long hash)
	{
		return 1.0/documentFrequencyRatio(hash);
	}
	
	public double maxIdf()
	{
		return this.maxIdfValue;
	}
	
	public double minIdf()
	{
		return this.minIdfValue;
	}
}
