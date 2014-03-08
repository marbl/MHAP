package com.secret.fastalign.main;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import com.secret.fastalign.utils.IntervalTree;
import com.secret.fastalign.utils.Utils;

public class EstimateROC {
	private static class Pair {
		public int first;
		public int second;

		public Pair(int first, int second) {
			this.first = first;
			this.second = second;
		}

		@SuppressWarnings("unused")
		public int size() {
			return (Math.max(first, (int) second)
					- Math.min(first, (int) second) + 1);
		}
	}

	private static Random generator = null;
	public static int seed = 0;

	private HashMap<String, IntervalTree<Integer>> clusters = new HashMap<String, IntervalTree<Integer>>();
	private HashMap<String, String> seqToChr = new HashMap<String, String>();
	private HashMap<String, Pair> seqToPosition = new HashMap<String, Pair>();
	//private String[] seqToName = null;
	private HashMap<Integer, String> seqToName = new HashMap<Integer, String>();
	private HashSet<String> ovlNames = new HashSet<String>();
	private String[] ovlToName = null;
	
	private int minOvlLen = 500;
	private int numTrials = 1000;
	private long tp = 0;
	private long fn = 0;
	private long tn = 0;
	private long fp = 0;
	private double ppv = 0;

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
		int numTrials = 5000;

		if (args.length < 3) {
			printUsage();
			System.exit(1);
		}
		if (args.length > 3) {
			numTrials = Integer.parseInt(args[3]);
		}
		EstimateROC g = new EstimateROC(Integer.parseInt(args[2]), numTrials);

		// load and cluster reference
		System.err.print("Loading reference...");
		g.processReference(args[0]);
		System.err.println("done.");

		// load matches
		System.err.print("Loading matches...");
		g.processOverlaps(args[1]);
		System.err.println("done");

		System.err.println("Computing sensitivity...");
		g.estimateSensitivity();

		// now estimate FP/TN by picking random match and checking reference
		// mapping
		System.err.println("Computing specificity...");
		g.estimateSpecificity();

		// last but not least PPV, pick random subset of our matches and see what percentage are true
		System.err.println("Computing PPV...");
		g.estimatePPV();

		System.out.println("Estimated sensitivity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tp / (g.tp + g.fn)));
		System.out.println("Estimated specificity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tn / (g.fp + g.tn)));
		System.out.println("Estimated PPV:\t "
				+ Utils.DECIMAL_FORMAT.format(g.ppv));
	}

	@SuppressWarnings("unused")
	public EstimateROC(int minOvlLen, int numTrials) {
		this.minOvlLen = minOvlLen;
		this.numTrials = numTrials;
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

	private String getOvlName(String id, String id2) {
		return (id.compareTo(id2) <= 0 ? id + "_" + id2 : id2
				+ "_" + id);
	}
	private String pickRandomSequence() {
		int val = generator.nextInt(seqToName.size());
		return seqToName.get(val);
	}
	
	private String pickRandomMatch() {
		int val = generator.nextInt(ovlToName.length);
		return ovlToName[val];
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
			if (!chr.equalsIgnoreCase(chr2)) {
				System.err.println("Error: comparing wrong chromosomes!");
				System.exit(1);
			}
			int overlap = Utils.getRangeOverlap(p1.first, (int) p1.second,
					p2.first, (int) p2.second);
			if (overlap >= min && !id.equalsIgnoreCase(id2)) {
				result.add(id2);
			}
		}

		return result;
	}

	@SuppressWarnings("unused")
	private String[] getOverlapInfo(String line) {
		String[] result = new String[2];
		String[] splitLine = line.trim().split("\\s+");

		if (splitLine.length == 7) {
			result[0]= splitLine[0];
			result[1] = splitLine[1];
			double score = Double.parseDouble(splitLine[5]) * 5;
			int aoffset = Integer.parseInt(splitLine[3]);
			int boffset = Integer.parseInt(splitLine[4]);
			boolean isFwd = ("N".equals(splitLine[2]));
		} else if (splitLine.length == 13) {
			result[0] = splitLine[0];
			if (result[0].indexOf("/") != -1) {
				result[0] = result[0].substring(0, splitLine[0].indexOf("/"));
			}
			result[1] = splitLine[1];
		}
		
		return result;
	}
	private void processOverlaps(String file) throws Exception {
		BufferedReader bf = new BufferedReader(new InputStreamReader(
				new FileInputStream(file)));

		String line = null;
		int counter = 0;
		while ((line = bf.readLine()) != null) {
			String[] result = getOverlapInfo(line);
			String id = result[0];
			String id2 = result[1];

			if (id == null || id2 == null) {
				continue;
			}
			if (id.equalsIgnoreCase(id2)) {
				continue;
			}
			if (seqToChr.get(id) == null || seqToChr.get(id2) == null) {
				continue;
			}
			String ovlName = getOvlName(id, id2);
			if (ovlNames.contains(ovlName)) {
				continue;
			}
			ovlNames.add(ovlName);
			counter++;
		}
		ovlToName = ovlNames.toArray(new String[counter]);
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
			seqToName.put(counter, id);
			counter++;
		}
		bf.close();
		for (String chr : clusters.keySet()) {
			clusters.get(chr).build();
		}
	}

	private boolean overlapExists(String id, String id2) {
		return ovlNames.contains(getOvlName(id, id2));
	}

	private void checkMatches(String id, HashSet<String> matches) {
		for (String m : matches) {
			if (overlapExists(id, m)) {
				tp++;
			} else {
				fn++;
			}
		}
	}

	private void estimateSensitivity() {
		// we estimate TP/FN by randomly picking a sequence, getting its
		// cluster, and checking our matches
		for (int i = 0; i < numTrials; i++) {
			// pick cluster
			String id = pickRandomSequence();
			HashSet<String> matches = getSequenceMatches(id, minOvlLen);
			checkMatches(id, matches);
		}
	}

	private void estimateSpecificity() {
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

			if (overlapExists(id, other)) {
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
	}
	
	private void estimatePPV() {
		int numTP = 0;
		for (int i = 0; i < numTrials; i++) {
			// pick an overlap
			String[] ovl = pickRandomMatch().split("_");
			String id = ovl[0];
			String id2 = ovl[1];
			
			HashSet<String> matches = getSequenceMatches(id, 0);
			if (matches.contains(id2)) {
				numTP++;
			}
		}
		
		// now our formula for PPV. Estimate percent of our matches which are true
		ppv = (double)numTP / numTrials;
	}
}
