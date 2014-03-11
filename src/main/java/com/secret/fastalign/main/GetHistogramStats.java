package com.secret.fastalign.main;

import java.io.BufferedReader;
import java.util.TreeMap;

import com.secret.fastalign.utils.Utils;

public class GetHistogramStats {
	private static final int NUM_SD = 7;
	private TreeMap<Integer, Long> histogram = new TreeMap<Integer, Long>();
	private double percent = 0.99;
	private double mean = 0;
	private double stdev = 0;
	private int cut = 0;

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
		int total = 0;

		for (int val : this.histogram.keySet()) {
			long count = this.histogram.get(val);
			for (int i = 0; i < count; i++) {
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
