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
package edu.umd.marbl.mhap.main;

import java.io.BufferedReader;
import java.util.TreeMap;

import edu.umd.marbl.mhap.utils.Utils;

public class GetHistogramStats {
	private static final int NUM_SD = 7;
	private TreeMap<Integer, Long> histogram = new TreeMap<Integer, Long>();
	private double percent = 0.99;
	private double mean = 0;
	private double stdev = 0;
	private long cut = 0;

	public GetHistogramStats(String fileName, double p) {
		try {
			BufferedReader bf = Utils.getFile(fileName, "hist");
			String line = null;

			while ((line = bf.readLine()) != null) {
				String[] split = line.trim().split("\\s+");
				int val = Integer.parseInt(split[0]);
				long count = Long.parseLong(split[1]);
				this.histogram.put(val, count);
			}
			bf.close();
			this.percent = p;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void process() throws NumberFormatException {
		double variance = 0;
		double sum = 0;
		long total = 0;

		for (int val : this.histogram.keySet()) {
			long count = this.histogram.get(val);
			for (long i = 0; i < count; i++) {
				total++;
				double delta = (val - this.mean);
				this.mean += (delta / total);
				variance += delta * (val - this.mean);
				sum += val;
			}
		}
		variance /= total;
		this.stdev = Math.sqrt(variance);

		double runningSum = 0;
		for (int val : this.histogram.keySet()) {
			long count = this.histogram.get(val);
			runningSum += (double) val * count;
			if ((runningSum / sum) > this.percent) {
				this.cut = val;
				break;
			}
		}
	}

	@Override
	public String toString() {
		return Utils.DECIMAL_FORMAT.format(this.mean) + "\t" + Utils.DECIMAL_FORMAT.format(this.stdev) + "\t" + "\t" + this.cut
				+ "\t" + Utils.DECIMAL_FORMAT.format(this.mean + NUM_SD * this.stdev);
	}
	
	public static void main(String[] args) throws NumberFormatException {
		GetHistogramStats s = new GetHistogramStats(args[0], Double.parseDouble(args[1]));
		s.process();
		System.out.println(s.toString());
	}
}
