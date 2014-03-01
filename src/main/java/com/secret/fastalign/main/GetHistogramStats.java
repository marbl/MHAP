package com.secret.fastalign.main;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.TreeMap;

import com.secret.fastalign.utils.Utils;

public class GetHistogramStats {
	private static final NumberFormat nf = new DecimalFormat(
			"############.####");
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
				histogram.put(val, count);
			}
			bf.close();
			percent = p;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void process() throws NumberFormatException, IOException {
		double variance = 0;
		double sum = 0;
		int total = 0;

		for (int val : histogram.keySet()) {
			long count = histogram.get(val);
			for (int i = 0; i < count; i++) {
				total++;
				double delta = (val - mean);
				mean += (delta / total);
				variance += delta * (val - mean);
				sum += val;
			}
		}
		variance /= total;
		stdev = Math.sqrt(variance);

		double runningSum = 0;
		for (int val : histogram.keySet()) {
			long count = histogram.get(val);
			runningSum += (double) val * count;
			if ((runningSum / sum) > percent) {
				cut = val;
				break;
			}
		}
	}

	@Override
	public String toString() {
		return nf.format(mean) + "\t" + nf.format(stdev) + "\t" + "\t" + cut
				+ "\t" + nf.format(mean + NUM_SD * stdev);
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		GetHistogramStats s = new GetHistogramStats(args[0], Double.parseDouble(args[1]));
		s.process();
		System.out.println(s.toString());
	}
}
