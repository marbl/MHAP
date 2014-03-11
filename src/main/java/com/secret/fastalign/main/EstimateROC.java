package com.secret.fastalign.main;

import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.Sequence;
import com.secret.fastalign.utils.IntervalTree;
import com.secret.fastalign.utils.Utils;

public class EstimateROC {
	private static final double MIN_IDENTITY = 0.70;
	private static final int DEFAULT_NUM_TRIALS = 1000;
	private static final int DEFAULT_MIN_OVL = 500;
	
	private static class Pair {
		public int first;
		public int second;

		public Pair(int first, int second) {
			this.first = first;
			this.second = second;
		}

		@SuppressWarnings("unused")
		public int size() {
			return (Math.max(this.first, this.second)
					- Math.min(this.first, this.second) + 1);
		}
	}
	
	private static class Overlap {
		public int afirst;
		public int bfirst;
		public int asecond;
		public int bsecond;
		public boolean isFwd;
		public String id1;
		public String id2;

		public Overlap() {
			// do nothing
		}
	}

	private static Random generator = null;
	public static int seed = 0;

	private HashMap<String, IntervalTree<Integer>> clusters = new HashMap<String, IntervalTree<Integer>>();
	private HashMap<String, String> seqToChr = new HashMap<String, String>(10000000);
	private HashMap<String, Pair> seqToPosition = new HashMap<String, Pair>(10000000);
	private HashMap<Integer, String> seqToName = new HashMap<Integer, String>(10000000);
	private HashSet<String> ovlNames = new HashSet<String>(10000000*100);
	private HashMap<String, Overlap> ovlInfo = new HashMap<String, Overlap>(10000000*100);
	private HashMap<Integer, String> ovlToName = new HashMap<Integer, String>(10000000*100);
	
	private int minOvlLen = DEFAULT_MIN_OVL;
	private int numTrials = DEFAULT_NUM_TRIALS;
	private long tp = 0;
	private long fn = 0;
	private long tn = 0;
	private long fp = 0;
	private double ppv = 0;
	private Sequence[] dataSeq = null;

	public static void printUsage() {
		System.err
				.println("This program uses random sampling to estimate PPV/Sensitivity/Specificity");
		System.err.println("The program requires 2 arguments:");
		System.err
				.println("\t1. A blasr M4 file mapping sequences to a reference (or reference subset)");
		System.err
				.println("\t2. All-vs-all mappings of same sequences in CA ovl format");
		System.err.println("\t3. Minimum overlap length (default: " + DEFAULT_MIN_OVL);
		System.err.println("\t4. Number of random trials, 0 means full compute (default : " + DEFAULT_NUM_TRIALS);
		System.err.println("\t5. Sequences in fasta format.");
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			printUsage();
			System.exit(1);
		}
		EstimateROC g = null;
		if (args.length > 3) {
			g = new EstimateROC(Integer.parseInt(args[2]), Integer.parseInt(args[3]));
		} else if (args.length > 2) {
			g = new EstimateROC(Integer.parseInt(args[2]));
		} else {
			g = new EstimateROC();
		}
		
		System.err.println("Running, reference: " + args[0] + " matches: " + args[1]);
		System.err.println("Number trials:  " + (g.numTrials == 0 ? "all" : g.numTrials));
		System.err.println("Minimum ovl:  " + g.minOvlLen);
		
		// load and cluster reference
		System.err.print("Loading reference...");
		long startTime = System.nanoTime();
		long totalTime = startTime;
		g.processReference(args[0]);
		System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");

		// load matches
		System.err.print("Loading matches...");
		startTime = System.nanoTime();
		g.processOverlaps(args[1]);
		System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");

		if (args.length > 4) {
			// load fasta
			System.err.print("Loading fasta...");
			startTime = System.nanoTime();
			g.loadFasta(args[4]);
			System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
		}
		
		if (g.numTrials == 0) {
			System.err.print("Computing full statistics O(" + g.seqToName.size() + "^2) operations!...");
			startTime = System.nanoTime();
			g.fullEstimate();
			System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
		} else {
			System.err.print("Computing sensitivity...");
			startTime = System.nanoTime();
			g.estimateSensitivity();
			System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
	
			// now estimate FP/TN by picking random match and checking reference
			// mapping
			System.err.print("Computing specificity...");
			startTime = System.nanoTime();
			g.estimateSpecificity();
			System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
	
			// last but not least PPV, pick random subset of our matches and see what percentage are true
			System.err.print("Computing PPV...");
			startTime = System.nanoTime();
			g.estimatePPV();
			System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
		}
		System.err.println("Total time: " + (System.nanoTime() - totalTime) * 1.0e-9 + "s.");
		
		System.out.println("Estimated sensitivity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tp / (g.tp + g.fn)));
		System.out.println("Estimated specificity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tn / (g.fp + g.tn)));
		System.out.println("Estimated PPV:\t "
				+ Utils.DECIMAL_FORMAT.format(g.ppv));
	}

	public EstimateROC() {
		this(DEFAULT_MIN_OVL, DEFAULT_NUM_TRIALS);
	}

	public EstimateROC(int minOvlLen) {
		this(minOvlLen, DEFAULT_NUM_TRIALS);
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
		int val = generator.nextInt(this.seqToName.size());
		return this.seqToName.get(val);
	}
	
	private String pickRandomMatch() {
		int val = generator.nextInt(this.ovlToName.size());
		return this.ovlToName.get(val);
	}
	
	private int getOverlapSize(String id, String id2) {
		String chr = this.seqToChr.get(id);
		String chr2 = this.seqToChr.get(id2);
		Pair p1 = this.seqToPosition.get(id);
		Pair p2 = this.seqToPosition.get(id2);
		if (!chr.equalsIgnoreCase(chr2)) {
			System.err.println("Error: comparing wrong chromosomes!");
			System.exit(1);
		}
		return Utils.getRangeOverlap(p1.first, p1.second,
				p2.first, p2.second);
	}

	private HashSet<String> getSequenceMatches(String id, int min) {
		String chr = this.seqToChr.get(id);
		Pair p1 = this.seqToPosition.get(id);
		List<Integer> intersect = this.clusters.get(chr).get(p1.first,
				p1.second);
		HashSet<String> result = new HashSet<String>();

		Iterator<Integer> it = intersect.iterator();
		while (it.hasNext()) {
			String id2 = this.seqToName.get(it.next());
			Pair p2 = this.seqToPosition.get(id2);
			String chr2 = this.seqToChr.get(id2);
			if (!chr.equalsIgnoreCase(chr2)) {
				System.err.println("Error: comparing wrong chromosomes!");
				System.exit(1);
			}
			int overlap = Utils.getRangeOverlap(p1.first, p1.second,
					p2.first, p2.second);
			if (overlap >= min && !id.equalsIgnoreCase(id2)) {
				result.add(id2);
			}
		}

		return result;
	}

	@SuppressWarnings("unused")
	private Overlap getOverlapInfo(String line) {
		Overlap overlap = new Overlap();
		String[] splitLine = line.trim().split("\\s+");

		try {
			if (splitLine.length == 7 || splitLine.length == 6) {
				overlap.id1 = splitLine[0];
				overlap.id2 = splitLine[1];
				double score = Double.parseDouble(splitLine[5]) * 5;
				int aoffset = Integer.parseInt(splitLine[3]);
				int boffset = Integer.parseInt(splitLine[4]);
				boolean isFwd = ("N".equals(splitLine[2]));
			} else if (splitLine.length == 13) {
				overlap.afirst = Integer.parseInt(splitLine[5]);
				overlap.asecond = Integer.parseInt(splitLine[6]);
				overlap.bfirst = Integer.parseInt(splitLine[9]);
				overlap.bsecond = Integer.parseInt(splitLine[10]);
				overlap.isFwd = (Integer.parseInt(splitLine[8]) == 0);
				if (!overlap.isFwd) {
					overlap.bsecond = Integer.parseInt(splitLine[11]) - Integer.parseInt(splitLine[9]);
					overlap.bfirst = Integer.parseInt(splitLine[11]) - Integer.parseInt(splitLine[10]);
				}
				overlap.id1 = splitLine[0];
				if (overlap.id1.indexOf("/") != -1) {
					overlap.id1 = overlap.id1.substring(0,
							splitLine[0].indexOf("/"));
				}
				if (overlap.id1.indexOf(",") != -1) {
					overlap.id1 = overlap.id1.split(",")[1];
				}
				overlap.id2 = splitLine[1];
				if (overlap.id2.indexOf(",") != -1) {
					overlap.id2 = overlap.id2.split(",")[1];
				}
			}
		} catch (NumberFormatException e) {
			System.err.println("Warning: could not parse input line: " + line
					+ " " + e.getMessage());
		}
		
		return overlap;
	}
	
	private void loadFasta(String file) throws IOException {
		FastaData data = new FastaData(file, 0);
		data.enqueueFullFile();
		this.dataSeq = data.toArray();
	}
	
	private void processOverlaps(String file) throws Exception {
		BufferedReader bf = new BufferedReader(new InputStreamReader(
				new FileInputStream(file)));

		String line = null;
		int counter = 0;
		while ((line = bf.readLine()) != null) {
			Overlap ovl = getOverlapInfo(line);
			String id = ovl.id1;
			String id2 = ovl.id2;

			if (id == null || id2 == null) {
				continue;
			}
			if (id.equalsIgnoreCase(id2)) {
				continue;
			}
			if (this.seqToChr.get(id) == null || this.seqToChr.get(id2) == null) {
				continue;
			}
			String ovlName = getOvlName(id, id2);
			if (this.ovlNames.contains(ovlName)) {
				continue;
			}
			this.ovlNames.add(ovlName);
			this.ovlToName.put(counter, ovlName);
			this.ovlInfo.put(ovlName, ovl);
			counter++;
			
			if (counter % 100000 == 0) {
				System.err.println("Loaded " + counter);
			}
			
		}
		System.err.print("Processed " + this.ovlNames.size() + " overlaps");
		if (this.ovlNames.isEmpty()) {
			System.err
					.println("Error: No sequence matches to reference loaded!");
			System.exit(1);
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
			if (!this.clusters.containsKey(chr)) {
				this.clusters.put(chr, new IntervalTree<Integer>());
			}
			this.clusters.get(chr).addInterval(startInRef, endInRef,
					counter);
			this.seqToPosition.put(id, new Pair(startInRef, endInRef));
			this.seqToChr.put(id, chr);
			this.seqToName.put(counter, id);
			counter++;
		}
		bf.close();
		for (String chr : this.clusters.keySet()) {
			this.clusters.get(chr).build();
		}
		
		System.err.print("Processed " + this.clusters.size() + " chromosomes, "
				+ this.seqToPosition.size() + " sequences matching ref");
		if (this.seqToPosition.isEmpty()) {
			System.err
					.println("Error: No sequence matches to reference loaded!");
			System.exit(1);
		}
	}

	private boolean overlapExists(String id, String id2) {
		return this.ovlNames.contains(getOvlName(id, id2));
	}

	private void checkMatches(String id, HashSet<String> matches) {
		for (String m : matches) {
			if (overlapExists(id, m)) {
				this.tp++;
			} else {
				this.fn++;
			}
		}
	}
	
	private boolean computeDP(String id, String id2) {
		if (this.dataSeq == null) {
			return false;
		}
		Logger logger = Logger.getLogger(SmithWatermanGotoh.class.getName());
		logger.setLevel(Level.OFF);
		logger = Logger.getLogger(MatrixLoader.class.getName());
		logger.setLevel(Level.OFF);
		Overlap ovl = this.ovlInfo.get(getOvlName(id, id2));

		jaligner.Sequence s1 = new jaligner.Sequence(this.dataSeq[Integer.parseInt(ovl.id1)-1].toString().substring(ovl.afirst, ovl.asecond));
		jaligner.Sequence s2 = null;
		if (ovl.isFwd) {
			s2 = new jaligner.Sequence(this.dataSeq[Integer.parseInt(ovl.id2)-1].toString().substring(ovl.bfirst, ovl.bsecond));
		} else {
			s2 = new jaligner.Sequence(this.dataSeq[Integer.parseInt(ovl.id2)-1].getReverseCompliment().toString().substring(ovl.bfirst, ovl.bsecond));
		}
		Alignment alignment;
		try {
			alignment = SmithWatermanGotoh.align(s1, s2, MatrixLoader.load("IDENTITY"), 2f, 0f);
		} catch (MatrixLoaderException e) {
			return false;
		}
		return ((double)alignment.getSimilarity()/s1.length() > MIN_IDENTITY);
	}

	private void estimateSensitivity() {
		// we estimate TP/FN by randomly picking a sequence, getting its
		// cluster, and checking our matches
		for (int i = 0; i < this.numTrials; i++) {
			// pick cluster
			String id = pickRandomSequence();
			HashSet<String> matches = getSequenceMatches(id, this.minOvlLen);
			checkMatches(id, matches);
		}
	}

	private void estimateSpecificity() {
		long numFPCompared = 0;

		// we estimate FP/TN by randomly picking two sequences
		for (int i = 0; i < this.numTrials; i++) {
			// pick cluster
			String id = pickRandomSequence();
			String other = pickRandomSequence();
			while (id.equalsIgnoreCase(other)) {
				other = pickRandomSequence();
			}
			HashSet<String> matches = getSequenceMatches(id, 0);

			if (overlapExists(id, other)) {
				if (!matches.contains(other)) {
					this.fp++;
				}
				numFPCompared++;
			} else {
				if (!matches.contains(other)) {
					this.tn++;
				}
			}
		}
	}
	
	private void estimatePPV() {
		int numTP = 0;
		for (int i = 0; i < this.numTrials; i++) {
			// pick an overlap
			String[] ovl = pickRandomMatch().split("_");
			String id = ovl[0];
			String id2 = ovl[1];
			
			HashSet<String> matches = getSequenceMatches(id, 0);
			if (matches.contains(id2)) {
				numTP++;
			} else {
				if (computeDP(id, id2)) {
					numTP++;
				}
			}
		}
		
		// now our formula for PPV. Estimate percent of our matches which are true
		this.ppv = (double)numTP / this.numTrials;
	}
	
	private void fullEstimate() {
		for (int i = 0; i < this.seqToName.size(); i++) {
			String id = this.seqToName.get(i);
			for (int j = i+1; j < this.seqToName.size(); j++) {
				String id2 = this.seqToName.get(j);
				if (id == null || id2 == null) { continue; }
				HashSet<String> matches = getSequenceMatches(id, 0);

				if (!overlapExists(id, id2)) {
					if (!matches.contains(id2)) {
						this.tn++;
					} else if (getOverlapSize(id, id2) > this.minOvlLen) {
						this.fn++;
					}
				} else {
					if (matches.contains(id2)) {
						this.tp++;
					} else {
						if (computeDP(id, id2)) {
							this.tp++;
						} else {
							this.fp++;
						}
					}
				}
			}
		}
		this.ppv = (double)this.tp / (this.tp+this.fp);
	}
}
