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

import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.NeedlemanWunschGotoh;
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

import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.utils.IntervalTree;
import edu.umd.marbl.mhap.utils.Utils;

public class EstimateROC {
	private static final boolean ALIGN_SW = true;
	private static final double MIN_OVERLAP_DIFFERENCE = 0.8;
	private static final double MIN_IDENTITY = 0.70;
	private static final double MIN_REF_IDENTITY = MIN_IDENTITY + 0.10;
	private static final int DEFAULT_NUM_TRIALS = 10000;
	private static final int DEFAULT_MIN_OVL = 2000;
	private static final boolean DEFAULT_DO_DP = false;
	private static boolean DEBUG = false;
	
	private static class Pair {
		public int first;
		public int second;

		public Pair(int startInRef, int endInRef) {
			this.first = startInRef;
			this.second = endInRef;
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
		
		public int getSize() {
			double first = (double)Math.max(this.asecond, this.afirst) - (double)Math.min(this.asecond, this.afirst);
			first += (double)Math.max(this.bsecond, this.bfirst) - (double)Math.min(this.bsecond, this.bfirst);
			return (int)Math.round(first/2);
		}
		
		@Override
		public String toString() {
			StringBuilder stringBuilder = new StringBuilder();
			stringBuilder.append("Overlap Fwd=" + this.isFwd);
			stringBuilder.append(" Aid=");
			stringBuilder.append(this.id1);
			stringBuilder.append(" (");
			stringBuilder.append(this.afirst);
			stringBuilder.append(", ");
			stringBuilder.append(this.asecond);
			stringBuilder.append("), Bid=");
			stringBuilder.append(this.id2);
			stringBuilder.append(" (");
			stringBuilder.append(this.bfirst);
			stringBuilder.append(", ");
			stringBuilder.append(this.bsecond);
			stringBuilder.append(")");
			return stringBuilder.toString();
		}
	}

	private static Random generator = null;
	public static int seed = 0;

	private HashMap<String, IntervalTree<Integer>> clusters = new HashMap<String, IntervalTree<Integer>>();
	private HashMap<String, String> seqToChr = new HashMap<String, String>(10000000);
	private HashMap<String, Integer> seqToScore = new HashMap<String, Integer>(10000000);
	private HashMap<String, Pair> seqToPosition = new HashMap<String, Pair>(10000000);
	private HashMap<Integer, String> seqToName = new HashMap<Integer, String>(10000000);
	private HashMap<String, Integer> seqNameToIndex = new HashMap<String, Integer>(10000000);
	private HashSet<String> ovlNames = new HashSet<String>(10000000*10);
	private HashMap<String, Overlap> ovlInfo = new HashMap<String, Overlap>(10000000*10);
	private HashMap<Integer, String> ovlToName = new HashMap<Integer, String>(10000000*10);
	
	private int minOvlLen = DEFAULT_MIN_OVL;
	private int numTrials = DEFAULT_NUM_TRIALS;
	private boolean doDP = false;
	private long tp = 0;
	private long fn = 0;
	private long tn = 0;
	private long fp = 0;
	private double ppv = 0;
	private Sequence[] dataSeq = null;

	public static void printUsage() {
		System.err
				.println("This program uses random sampling to estimate PPV/Sensitivity/Specificity");
		System.err.println("The sequences in the fasta file used to generate the truth must be sequentially numbered from 1 to N!");
		System.err
				.println("\t1. A blasr M4 file mapping sequences to a reference (or reference subset)");
		System.err
				.println("\t2. All-vs-all mappings of same sequences in CA ovl format");
		System.err
		.println("\t3. Fasta sequences sequentially numbered from 1 to N.");
		System.err.println("\t4. Minimum overlap length (default: " + DEFAULT_MIN_OVL);
		System.err.println("\t5. Number of random trials, 0 means full compute (default : " + DEFAULT_NUM_TRIALS);
		System.err.println("\t6. Compute DP during PPV true/false");
		System.err.println("\t7. Debug output true/false");
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 3) {
			printUsage();
			System.exit(1);
		}
		EstimateROC g = null;
		if (args.length > 5) {
			g = new EstimateROC(Integer.parseInt(args[3]), Integer.parseInt(args[4]), Boolean.parseBoolean(args[5]));
		} else if (args.length > 4) {
			g = new EstimateROC(Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		} else if (args.length > 3) {
			g = new EstimateROC(Integer.parseInt(args[3]));
		} else {
			g = new EstimateROC();
		}
		if (args.length > 6) {
			DEBUG = Boolean.parseBoolean(args[6]);
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

		// load fasta
		System.err.print("Loading fasta...");
		startTime = System.nanoTime();
		g.loadFasta(args[2]);
		System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
		
		// load matches
		System.err.print("Loading matches...");
		startTime = System.nanoTime();
		g.processOverlaps(args[1]);
		System.err.println("done " + (System.nanoTime() - startTime) * 1.0e-9 + "s.");
		
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
				+ Utils.DECIMAL_FORMAT.format((double) g.tp / (double)(g.tp + g.fn)));
		System.out.println("Estimated specificity:\t"
				+ Utils.DECIMAL_FORMAT.format((double) g.tn / (double)(g.fp + g.tn)));
		System.out.println("Estimated PPV:\t "
				+ Utils.DECIMAL_FORMAT.format(g.ppv));
	}

	public EstimateROC() {
		this(DEFAULT_MIN_OVL, DEFAULT_NUM_TRIALS);
	}

	public EstimateROC(int minOvlLen) {
		this(minOvlLen, DEFAULT_NUM_TRIALS);
	}
	
	public EstimateROC(int minOvlLen, int numTrials) {
		this(minOvlLen, numTrials, DEFAULT_DO_DP);
	}
	
	@SuppressWarnings("unused")
	public EstimateROC(int minOvlLen, int numTrials, boolean doDP) {
		this.minOvlLen = minOvlLen;
		this.numTrials = numTrials;
		this.doDP = doDP;
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

	private static int getSequenceId(String id) {
		return Integer.parseInt(id)-1;
	}
	
	private static String getOvlName(String id, String id2) {
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
			System.err.println("Error: comparing wrong chromosomes betweeen sequences " + id + " and sequence " + id2);
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
				System.err.println("Error: comparing wrong chromosomes betweeen sequences " + id + " and sequence in its cluster " + id2);
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

	private Overlap getOverlapInfo(String line) {
		Overlap overlap = new Overlap();
		String[] splitLine = line.trim().split("\\s+");

		try {
			// CA format
			if (splitLine.length == 7 || splitLine.length == 6) {
				overlap.id1 = splitLine[0];
				overlap.id2 = splitLine[1];
				@SuppressWarnings("unused")
				double score = Double.parseDouble(splitLine[5]) * 5;
				int aoffset = Integer.parseInt(splitLine[3]);
				int boffset = Integer.parseInt(splitLine[4]);
				overlap.isFwd = "N".equalsIgnoreCase(splitLine[2]);
				if (this.dataSeq != null) {
					int alen = this.dataSeq[Integer.parseInt(overlap.id1)-1].length();
					int blen = this.dataSeq[Integer.parseInt(overlap.id2)-1].length();
					overlap.afirst = Math.max(0, aoffset);
					overlap.asecond = Math.min(alen, alen + boffset);
					overlap.bfirst = -1*Math.min(0, aoffset);
					overlap.bsecond = Math.min(blen, blen - boffset);
				}
				//mhap format
			} else if (splitLine.length == 12) {
				overlap.id1 = splitLine[0];
				overlap.id2 = splitLine[1];
				@SuppressWarnings("unused")
				double score = Double.parseDouble(splitLine[2]);
				overlap.isFwd = Integer.parseInt(splitLine[8]) == 0;
				if (this.dataSeq != null) {
					int alen = this.dataSeq[getSequenceId(overlap.id1)].length();
					int blen = this.dataSeq[getSequenceId(overlap.id2)].length();
					overlap.afirst = Integer.parseInt(splitLine[5]);
					overlap.asecond = Integer.parseInt(splitLine[6]);
					overlap.bfirst = Integer.parseInt(splitLine[9]);
					overlap.bsecond = Integer.parseInt(splitLine[10]);
					if (overlap.asecond > alen) {
						overlap.asecond = alen;
					}
					if (overlap.bsecond > blen) {
						overlap.bsecond = blen;
					}
				}
				// blasr format
			} else if (splitLine.length == 13 && !line.contains("[")) {
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
				if (this.dataSeq != null) {
					int alen = this.dataSeq[getSequenceId(overlap.id1)].length();
					int blen = this.dataSeq[getSequenceId(overlap.id2)].length();
					if (overlap.asecond > alen) {
						overlap.asecond = alen;
					}
					if (overlap.bsecond > blen) {
						overlap.bsecond = blen;
					}
				}
				//         1      1,182 n   [ 4,746.. 8,108] x [     0.. 3,896] :   <    982 diffs  ( 34 trace pts)
			} else if (splitLine.length >= 13 && splitLine.length <= 18) {
				overlap.id1 = splitLine[0].replaceAll(",", "");
				overlap.id2 = splitLine[1].replaceAll(",", "");	
				overlap.isFwd = (splitLine[2].equalsIgnoreCase("n"));
				String[] splitTwo = line.split("\\[");
				String aInfo = splitTwo[1].substring(0, splitTwo[1].indexOf("]"));
				String bInfo = splitTwo[2].substring(0, splitTwo[2].indexOf("]"));
				String[] aSplit = aInfo.replaceAll(",", "").split("\\.\\.");
				String[] bSplit = bInfo.replaceAll(",", "").split("\\.\\.");
				overlap.afirst=Integer.parseInt(aSplit[0].trim());
				overlap.asecond=Integer.parseInt(aSplit[1].trim());
				overlap.bfirst=Integer.parseInt(bSplit[0].trim());
				overlap.bsecond=Integer.parseInt(bSplit[1].trim());
				if (!overlap.isFwd) {
					overlap.bsecond = this.dataSeq[getSequenceId(overlap.id2)].length() - Integer.parseInt(bSplit[0].trim());
					overlap.bfirst = this.dataSeq[getSequenceId(overlap.id2)].length() - Integer.parseInt(bSplit[1].trim());
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
		this.dataSeq = new Sequence[data.getNumberProcessed()];
		int i = 0;
		while (!data.isEmpty()) {
			this.dataSeq[i++] = data.dequeue();
		}
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
		
		bf.close();
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
			if (id.indexOf(",") != -1) {
				id = id.split(",")[1];
			}
			double idy = Double.parseDouble(splitLine[3]);
			int start = Integer.parseInt(splitLine[5]);
			int end = Integer.parseInt(splitLine[6]);
			int length = Integer.parseInt(splitLine[7]);
			int seqIsFwd = Integer.parseInt(splitLine[4]);
			if (seqIsFwd != 0) {
				System.err.println("Error: malformed line, first sequences should always be in fwd orientation");
				System.exit(1);
			}
			int startInRef = Integer.parseInt(splitLine[9]);
			int endInRef = Integer.parseInt(splitLine[10]);
			int refLen = Integer.parseInt(splitLine[11]);
			int isRev = Integer.parseInt(splitLine[8]);
			int score = Integer.parseInt(splitLine[2]);
			if (isRev == 1) {
				int tmp = refLen - endInRef;
				endInRef = refLen - startInRef;
				startInRef = tmp;
			}
			if (idy < MIN_REF_IDENTITY*100) {
				continue;
			}
			double diff = ((double)(end - start) / (double)(endInRef-startInRef));
			if (diff < MIN_OVERLAP_DIFFERENCE) {
				continue;
			}
			String chr = splitLine[1];
			if (!this.clusters.containsKey(chr)) {
				this.clusters.put(chr, new IntervalTree<Integer>());
			}
			if (this.seqToPosition.containsKey(id)) {
				if (score < this.seqToScore.get(id)) {
					// replace
					this.seqToPosition.put(id, new Pair(startInRef, endInRef));
					this.seqToChr.put(id, chr);
					this.seqToScore.put(id, score);
				}
			} else {
				this.seqToPosition.put(id, new Pair(startInRef, endInRef));
				this.seqToChr.put(id, chr);
				this.seqToName.put(counter, id);
				this.seqNameToIndex.put(id, counter);
				this.seqToScore.put(id, score);
				counter++;
			}
		}
		bf.close();
		for (String id : this.seqToPosition.keySet()) {
			String chr = this.seqToChr.get(id);
			if (!this.clusters.containsKey(chr)) {
				this.clusters.put(chr, new IntervalTree<Integer>());
			}
			Pair p = this.seqToPosition.get(id);
			this.clusters.get(chr).addInterval(p.first, p.second,
					this.seqNameToIndex.get(id));
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
	
	private boolean overlapMatches(String id, String m) {
		int refOverlap = getOverlapSize(id, m);
		Overlap ovl = this.ovlInfo.get(getOvlName(id, m));
		if (ovl == null) {
			return false;
		}
		int diff = Math.abs(ovl.getSize() - refOverlap);
		double diffPercent = (double)diff / (double)refOverlap;
		if (DEBUG) { System.err.println("Overlap " + ovl + " " + ovl.getSize() + " versus ref " + refOverlap + " " + " diff is " + diff + "(" + diffPercent + ")"); }
		if (diffPercent > 0.3) {
			return false;
		}
		return true;
	}

	private void checkMatches(String id, HashSet<String> matches) {
		for (String m : matches) {
			if (overlapMatches(id, m)) {
				this.tp++;
			} else {
				this.fn++;
				if (DEBUG) { 
					System.err.println("Overlap between sequences: " + id + ", " + m + " is missing.");
					System.err.println(">" + id + " reference location " + this.seqToChr.get(id) + " " + this.seqToPosition.get(id).first + ", " + this.seqToPosition.get(id).second);
					System.err.println(this.dataSeq[Integer.parseInt(id)-1].getString());
					System.err.println(">" + m + " reference location " + this.seqToChr.get(m) + " " + this.seqToPosition.get(m).first + ", " + this.seqToPosition.get(m).second);
					System.err.println(this.dataSeq[Integer.parseInt(m)-1].getString());
				}
			}
		}
	}
	
	private static double getScore(Alignment alignment) {
		char[] sequence1 = alignment.getSequence1();
		char[] sequence2 = alignment.getSequence2();
		int length = Math.max(sequence1.length, sequence2.length);
		int ovlLen = Math.min(sequence1.length, sequence2.length);
		char GAP = '-';
		//int errors = 0;
		int matches = 0;
		for (int i = 0; i <= length; i++)
		{
			char c1 = GAP;
			char c2 = GAP;
			if (i < sequence1.length) {
				c1 = sequence1[i];
			}
			if (i < sequence2.length) {
				c2 = sequence2[i];
			}
			if (c1 != c2 || c1 == GAP || c2 == GAP) {
				//errors++;
			} else {
				matches++;
			}
		}
		return (matches / (double)ovlLen);
	}
	
	private boolean computeDP(String id, String id2) {
		if (this.doDP == false) {
			return false;
		}
		Logger logger = null;
		if (ALIGN_SW) {
			logger = Logger.getLogger(SmithWatermanGotoh.class.getName());
		} else {
			logger = Logger.getLogger(NeedlemanWunschGotoh.class.getName());
		}
		logger.setLevel(Level.OFF);
		logger = Logger.getLogger(MatrixLoader.class.getName());
		logger.setLevel(Level.OFF);
		Overlap ovl = this.ovlInfo.get(getOvlName(id, id2));
		System.err.println("Aligning sequence " + ovl.id1 + " to " + ovl.id2 + " " + ovl.bfirst + " to " + ovl.bsecond + " and " + ovl.isFwd + " and " + ovl.afirst + " " + ovl.asecond);

		jaligner.Sequence s1 = new jaligner.Sequence(this.dataSeq[getSequenceId(ovl.id1)].getString().substring(ovl.afirst, ovl.asecond));
		jaligner.Sequence s2 = null;
		if (ovl.isFwd) {
			s2 = new jaligner.Sequence(this.dataSeq[getSequenceId(ovl.id2)].getString().substring(ovl.bfirst, ovl.bsecond));
		} else {
			s2 = new jaligner.Sequence(Utils.rc(this.dataSeq[getSequenceId(ovl.id2)].getString().substring(ovl.bfirst, ovl.bsecond)));
		}
		Alignment alignment;
		try {
			if (ALIGN_SW) {
				alignment = SmithWatermanGotoh.align(s1, s2, MatrixLoader.load("MATCH"), 2f, 1f);
			} else {
				alignment = NeedlemanWunschGotoh.align(s1, s2, MatrixLoader.load("MATCH"), 2f, 1f);
			}
		} catch (MatrixLoaderException e) {
			return false;
		}
		double score = getScore(alignment); // alignment.getIdentity() / 100;
		if (DEBUG) { 
			System.err.println(alignment.getSummary());
			System.err.println("My score: " + score);
			System.err.println (new jaligner.formats.Pair().format(alignment)); 
		}
		return (score > MIN_IDENTITY && alignment.getLength() > this.minOvlLen);
	}

	private void estimateSensitivity() {
		// we estimate TP/FN by randomly picking a sequence, getting its
		// cluster, and checking our matches
		for (int i = 0; i < this.numTrials; i++) {
			String id = null;
			HashSet<String> matches = null;
			while (matches == null || matches.size() == 0) {
				// pick cluster
				id = pickRandomSequence();
				matches = getSequenceMatches(id, this.minOvlLen);
			}
			
			if (DEBUG) { System.err.println("Estimated sensitivity trial #" + i + " " + id + " matches " + matches); }
			checkMatches(id, matches);
		}
	}

	private void estimateSpecificity() {
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
			int ovlLen = 0;
			String[] ovl = null;
			String ovlName = null;
			while (ovlLen < this.minOvlLen) {
				// pick an overlap
				ovlName = pickRandomMatch();
				Overlap o = this.ovlInfo.get(ovlName);
				ovlLen = Utils.getRangeOverlap(o.afirst, o.asecond, o.bfirst, o.bsecond);
			}
			if (ovlName == null) {
				System.err.println("Could not find any computed overlaps > " + this.minOvlLen);
				System.exit(1);
			} else {
				ovl = ovlName.split("_");
				String id = ovl[0];
				String id2 = ovl[1];
				
				HashSet<String> matches = getSequenceMatches(id, 0);
				if (matches.contains(id2)) {
					numTP++;
				} else {
					if (computeDP(id, id2)) {
						numTP++;
					} else {
						if (DEBUG) { System.err.println("Overlap between sequences: " + id + ", " + id2 + " is not correct."); }
					}
				}
			}
		}
		
		// now our formula for PPV. Estimate percent of our matches which are true
		this.ppv = (double)numTP / (double)this.numTrials;
	}
	
	@SuppressWarnings("cast")
	private void fullEstimate() {
		for (int i = 0; i < this.seqToName.size(); i++) {
			String id = this.seqToName.get(i);
			for (int j = i+1; j < this.seqToName.size(); j++) {
				String id2 = this.seqToName.get(j);
				if (id == null || id2 == null) { continue; }
				HashSet<String> matches = getSequenceMatches(id, 0);

				if (!overlapMatches(id, id2)) {
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
		this.ppv = (double) this.tp / ((double)this.tp+(double)this.fp);
	}
}
