package com.secret.fastalign.main;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.GregorianCalendar;
import java.util.Calendar;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.HashSet;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import com.secret.fastalign.utils.IntervalTree;
import com.secret.fastalign.utils.Utils;

public class EstimateROC {
	private static class Pair {
		public int first;
		public int second;
		public String identifier;

		public Pair(int first, int second) {
			this.first = first;
			this.second = second;
		}

		public Pair(int first, int second, String third) {
			this.first = first;
			this.second = second;
			this.identifier = third;
		}

		public int size() {
			return (Math.max(first, (int) second)
					- Math.min(first, (int) second) + 1);
		}
	}

	private static Random generator = null;
	public static int seed = 0;

	private HashMap<String, IntervalTree<Integer>> clusters = new HashMap<String, IntervalTree<Integer>>();
	private HashMap<String, String> seqToChr = new HashMap<String, String>();
	private HashSet<String> chrs = new HashSet<String>();
	private TreeMap<String, Pair> seqToPosition = new TreeMap<String, Pair>();
	private HashMap<Integer, String> seqToName = new HashMap<Integer, String>();
	private HashMap<String, HashSet<String>> overlaps = new HashMap<String, HashSet<String>>();
	private HashSet<String> ovlNames = new HashSet<String>();

	private int minOvlLen = 500;
	private long tp = 0;
	private long fn = 0;
	private long tn = 0;
	private long fp = 0;
	private double ppv = 0;
	private long numMatches = 0;
	private long numCompared = 0;

	public static void printUsage() {
		System.err
				.println("This program uses random sampling to estimate PPV/Sensitivity/Specificity");
		System.err.println("The program requires 3 arguments:");
		System.err
				.println("\t1. A blasr M4 file mapping sequences to a reference (or reference subset)");
		System.err
				.println("\t2. All-vs-all mappings of same sequences in CA ovl format");
		System.err.println("\t3. Minimum overlap length");
	}

	public static void main(String[] args) throws Exception {
		int NUM_TP_TRIALS = 5000;
		int NUM_FP_TRIALS = NUM_TP_TRIALS;

		if (args.length < 3) {
			printUsage();
			System.exit(1);
		}
		EstimateROC g = new EstimateROC(Integer.parseInt(args[2]));

		// load and cluster reference
		g.processReference(args[0]);

		// load matches
		g.processOverlaps(args[1]);

		g.estimateSensitivity(NUM_TP_TRIALS);

		// now estimate FP/TN by picking random match and checking reference
		// mapping
		g.estimateSpecificity(NUM_FP_TRIALS);

		System.out.println("Estimated sensitivity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tp / (g.tp + g.fn)));
		System.out.println("Estimated specificity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tn / (g.fp + g.tn)));
		System.out.println("Estimated PPV:\t "
				+ Utils.DECIMAL_FORMAT.format(g.ppv));
	}

	@SuppressWarnings("unused")
	public EstimateROC(int minOvlLen) {
		this.minOvlLen = minOvlLen;
		if (false) {
			GregorianCalendar t = new GregorianCalendar();
			int t1 = t.get(Calendar.SECOND);
			int t2 = t.get(Calendar.MINUTE);
			int t3 = t.get(Calendar.HOUR_OF_DAY);
			int t4 = t.get(Calendar.DAY_OF_MONTH);
			int t5 = t.get(Calendar.MONTH);
			int t6 = t.get(Calendar.YEAR);
			seed = t6 + 65 * (t5 + 12 * (t4 + 31 * (t3 + 24 * (t2 + 60 * t1))));
		}

		generator = new Random(seed);
	}

	private String pickRandomSequence() {
		int val = generator.nextInt(seqToName.size());
		return seqToName.get(val);
	}

	private HashSet<String> getSequenceMatches(String id, int min) {
		String chr = seqToChr.get(id);
		Pair p1 = seqToPosition.get(id);
		List<Integer> intersect = clusters.get(chr).get(p1.first,
				(long) p1.second);
		HashSet<String> result = new HashSet<String>();

		Iterator<Integer> it = intersect.iterator();
		while (it.hasNext()) {
			String id2 = seqToName.get(it.next());
			Pair p2 = seqToPosition.get(id2);
			String chr2 = seqToChr.get(id2);
			int overlap = Utils.getRangeOverlap(p1.first, (int) p1.second,
					p2.first, (int) p2.second);
			if (overlap >= min && !id.equalsIgnoreCase(id2)) {
				result.add(id2);
			}
		}

		return result;
	}

	@SuppressWarnings("unused")
	private void processOverlaps(String file) throws Exception {
		BufferedReader bf = new BufferedReader(new InputStreamReader(
				new FileInputStream(file)));

		String line = null;
		while ((line = bf.readLine()) != null) {
			String[] splitLine = line.trim().split("\\s+");

			String id = null;
			String id2 = null;
			if (splitLine.length == 7) {
				id = splitLine[0];
				id2 = splitLine[1];
				double score = Double.parseDouble(splitLine[5]) * 5;
				int aoffset = Integer.parseInt(splitLine[3]);
				int boffset = Integer.parseInt(splitLine[4]);
				boolean isFwd = ("N".equals(splitLine[2]));
			} else if (splitLine.length == 13) {
				id = splitLine[0];
				if (id.indexOf("/") != -1) {
					id = id.substring(0, splitLine[0].indexOf("/"));
				}
				id2 = splitLine[1];
			}

			if (id == null || id2 == null) {
				continue;
			}
			if (id.equalsIgnoreCase(id2)) {
				continue;
			}
			if (seqToChr.get(id) == null || seqToChr.get(id2) == null) {
				continue;
			}
			String ovlName = (id.compareTo(id2) <= 0 ? id + "_" + id2 : id2
					+ "_" + id);
			if (ovlNames.contains(ovlName)) {
				continue;
			}
			ovlNames.add(ovlName);
			numMatches++;
			if (overlaps.get(id) == null) {
				overlaps.put(id, new HashSet<String>());
			}
			overlaps.get(id).add(id2);
			if (overlaps.get(id2) == null) {
				overlaps.put(id2, new HashSet<String>());
			}
			overlaps.get(id2).add(id);
		}
	}

	/**
	 * We are parsing file of the format 18903/0_100 ref000001|lambda_NEB3011
	 * -462 96.9697 0 0 99 100 0 2 101 48502 254 21589/0_100
	 * ref000001|lambda_NEB3011 -500 100 0 0 100 100 1 4 104 48502 254
	 * 15630/0_100 ref000001|lambda_NEB3011 -478 98 0 0 100 100 0 5 105 48502
	 * 254
	 **/
	@SuppressWarnings("unused")
	private void processReference(String file) throws Exception {
		BufferedReader bf = new BufferedReader(new InputStreamReader(
				new FileInputStream(file)));
		String line = null;
		int counter = 0;
		while ((line = bf.readLine()) != null) {
			String[] splitLine = line.trim().split("\\s+");

			String id = splitLine[0];
			if (id.indexOf("/") != -1) {
				id = id.substring(0, splitLine[0].indexOf("/"));
			}
			int start = Integer.parseInt(splitLine[5]);
			int end = Integer.parseInt(splitLine[6]);
			int length = Integer.parseInt(splitLine[7]);
			int startInRef = Integer.parseInt(splitLine[9]);
			int endInRef = Integer.parseInt(splitLine[10]);
			String chr = splitLine[1];
			if (!clusters.containsKey(chr)) {
				clusters.put(chr, new IntervalTree<Integer>());
			}
			clusters.get(chr).addInterval((long) startInRef, (long) endInRef,
					counter);
			seqToPosition.put(id, new Pair(startInRef, endInRef));
			seqToChr.put(id, chr);
			chrs.add(chr);
			seqToName.put(counter, id);
			counter++;
		}
		bf.close();
		for (String chr : clusters.keySet()) {
			clusters.get(chr).build();
		}
	}

	private void checkMatches(String id, HashSet<String> matches) {
		for (String m : matches) {
			numCompared++;
			if (overlaps.get(id).contains(m) || overlaps.get(m).contains(id)) {
				tp++;
			} else {
				fn++;
			}
		}
	}

	private void estimateSensitivity(int numTrials) {
		// we estimate TP/FN by randomly picking a sequence, getting its
		// cluster, and checking our matches
		for (int i = 0; i < numTrials; i++) {
			// pick cluster
			String id = pickRandomSequence();
			HashSet<String> matches = getSequenceMatches(id, minOvlLen);
			checkMatches(id, matches);
		}
	}

	private void estimateSpecificity(int numTrials) {
		long numFPCompared = 0;

		// we estimate FP/TN by randomly picking two seqeunces
		for (int i = 0; i < numTrials; i++) {
			// pick cluster
			String id = pickRandomSequence();
			String other = pickRandomSequence();
			while (id.equalsIgnoreCase(other)) {
				other = pickRandomSequence();
			}
			HashSet<String> matches = getSequenceMatches(id, 0);

			if (overlaps.get(id).contains(other)
					|| overlaps.get(other).contains(id)) {
				if (!matches.contains(other)) {
					fp++;
				}
				numFPCompared++;
			} else {
				if (!matches.contains(other)) {
					tn++;
				}
			}
		}

		// now our formula for PPV. We want to estimate the number of true
		// matches of our total matches
		// we know # total matches
		double numEstimatedTP = ((double) tp / numCompared) * numMatches;
		double numEstimatedFP = ((double) fp / numFPCompared) * numMatches;
		ppv = numEstimatedTP / (numEstimatedTP + numEstimatedFP);
	}
}
