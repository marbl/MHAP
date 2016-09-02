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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import com.google.common.hash.BloomFilter;

import edu.umd.marbl.mhap.impl.MhapRuntimeException;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.longs.Long2DoubleOpenHashMap;

public final class FrequencyCounts
{
	private final double filterCutoff;
	private final Map<Long,Double> fractionCounts;
	private final Set<Integer> kmerSizes;
	private final double maxIdfValue;
	private final double maxValue;
	private final double minIdfValue;
	private final double minValue;
	private final boolean noTf;
	private final double offset;
	private final int removeUnique;
	private final BloomFilter<Long> validMers;
	
	public static final double REPEAT_SCALE = 3.0;
	
	public FrequencyCounts(BufferedReader bf, double filterCutoff, double offset, int removeUnique, boolean noTf, int numThreads) throws IOException
	{
		//removeUnique = 0: do nothing extra to k-mers not specified in the file
		//removeUnique = 1: remove k-mers not specified in the file from the sketch
		//removeUnique = 2: supress k-mers not specified in the file the same as max supression
		
		if (removeUnique<0 || removeUnique>2)
			throw new MhapRuntimeException("Unknown removeUnique option "+removeUnique+".");
		
		if (offset<0.0 || offset>=1.0)
			throw new MhapRuntimeException("Offset can only be between 0 and 1.0.");

		this.kmerSizes = new IntOpenHashSet();
		this.removeUnique = removeUnique;
		this.noTf = noTf;
		
		// generate hashset
		Long2DoubleOpenHashMap validMap = new Long2DoubleOpenHashMap();
		BloomFilter<Long> validMers;

		//the max value observed in the list
		AtomicReference<Double> maxValue = new AtomicReference<Double>(Double.NEGATIVE_INFINITY);

		//read in the first line to generate the bloom filter
		String line = bf.readLine();
		try
		{
			long size;
			if (line==null)
			{
				System.err.println("Warning, k-mer filter file is empty. Assuming zero entries.");
				size = 1L;
			}
			else
			{
				size = Long.parseLong(line);
			
				if (size<0L)
					throw new MhapRuntimeException("K-mer filter file size line must have positive long value.");
				else
				if (size==0L)
				{
					System.err.println("Warning, k-mer filter file has zero elements.");
					size = 1L;
				}
			}
			
			//if no nothing, no need to store the while list
			if (removeUnique>0)
				validMers = BloomFilter.create((value, sink) -> sink.putLong(value), size, 1.0e-5);
			else
				validMers = null;
		}
		catch (Exception e)
		{
			throw new MhapRuntimeException("K-mer filter file first line must contain estimated number of k-mers in the file (long).");
		}
		
		final ThreadPoolExecutor executor = new ThreadPoolExecutor(numThreads, numThreads, 100L, TimeUnit.MILLISECONDS,
				new LinkedBlockingQueue<Runnable>(10000), new ThreadPoolExecutor.CallerRunsPolicy());
		
		line = bf.readLine();			
		while (line != null)
		{
			String currLine = line;
			
			executor.submit(() -> 
			{
				try
				{
					String[] str = currLine.split("\\s+", 3);
					
					if (str.length < 1)
						throw new MhapRuntimeException("K-mer filter file must have at least one column [k-mer]. Line="+currLine);

					//store the kmer sizes in the list
					synchronized (this.kmerSizes)
					{
						this.kmerSizes.add(str[0].length());
					}					
					
					long[] hash = HashUtils.computeSequenceHashesLong(str[0], str[0].length(), 0);
					
					if (str.length >= 2)
					{
						double percent = Double.parseDouble(str[1]);
						
						// if greater, add to hashset
						if (percent >= filterCutoff)
						{
							maxValue.getAndUpdate(v -> Math.max(v, percent));
							
							//store the max percent
							synchronized (validMap)
							{
								validMap.put(hash[0], percent);								
							}
						}
					}
		
					//store in the bloom filter
					if (removeUnique>0)
						synchronized (validMers)
						{
							validMers.put(hash[0]);							
						}
				}
				catch (Exception e)
				{
					System.err.println(e);
				}
	
			});
			
			// read the next line
			line = bf.readLine();									
		}
		
		executor.shutdown();
		try
		{
			executor.awaitTermination(5L, TimeUnit.DAYS);
		}
		catch (InterruptedException e)
		{
			executor.shutdownNow();
			throw new RuntimeException("Unable to finish all tasks.");
		}
		
		//trim the hashtable to the right size
		validMap.trim();
	
		this.validMers = validMers;
		this.fractionCounts = validMap;
		this.filterCutoff = filterCutoff;
		this.offset = offset;
		this.maxValue = maxValue.get();
		this.minValue = this.filterCutoff;
		
		this.minIdfValue = idf(this.maxValue);
		this.maxIdfValue = idf(this.minValue);
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
	
	public List<Integer> getKmerSizes()
	{
		return new ArrayList<>(this.kmerSizes);
	}
	
	public double idf(double freq)
	{
		return Math.log(this.maxValue/freq-this.offset);
		//return Math.log1p(this.maxValue/freq);
	}
	
	public double idf(long hash)
	{
		double freq = documentFrequencyRatio(hash);
		return idf(freq); 
	}
	
	public double inverseDocumentFrequency(long hash)
	{
		return 1.0/documentFrequencyRatio(hash);
	}
	
	public boolean isPopular(long hash)
	{
		return this.fractionCounts.containsKey(hash);
	}

	public boolean keepKmer(long hash)
	{
		if (this.removeUnique==1)
			return this.validMers.mightContain(hash);
			
		return true;
	}
	
	public double maxIdf()
	{
		return this.maxIdfValue;
	}
	
	public double minIdf()
	{
		return this.minIdfValue;
	}
	
	public double scaledIdf(long hash)
	{
		return scaledIdf(hash, REPEAT_SCALE);
	}
	
	public double scaledIdf(long hash, double maxValue)
	{
		if (this.removeUnique==2 && this.validMers!=null && !this.validMers.mightContain(hash))
			return 1.0;			
		
		Double val = this.fractionCounts.get(hash);
		if (val == null)
			return maxValue;
		
		//get the true value
		double idf = idf(val);
		
		//scale it to match max
		double scale = (maxIdf()-minIdf())/(double)(maxValue-1.0);
		
		return 1.0+(idf-minIdf())/scale;
	}

	public double tfWeight(int weight)
	{
		if (this.noTf)
			return 1.0;
		
		return (double)weight;
	}
}
