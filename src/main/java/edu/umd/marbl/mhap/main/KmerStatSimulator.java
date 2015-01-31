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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Random;
import java.util.GregorianCalendar;
import java.util.Calendar;
import java.util.HashSet;
import java.io.BufferedReader;
import java.io.PrintStream;

import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.general.SequenceId;
import edu.umd.marbl.mhap.sketch.MinHash;
import edu.umd.marbl.mhap.utils.Utils;

public class KmerStatSimulator {
	private boolean verbose = false;
	private int kmer = -1;
	private int overlap = 100;

	private ArrayList<Double> randomJaccard = new ArrayList<Double>();
	private ArrayList<Double> randomMinHash = new ArrayList<Double>();
	private ArrayList<Double> randomMerCounts = new ArrayList<Double>();
	private String reference = null;
	private double requestedLength = 5000;

	private double sharedCount = 0;
	private ArrayList<Double> sharedJaccard = new ArrayList<Double>();
	private ArrayList<Double> sharedMinHash = new ArrayList<Double>();
	private ArrayList<Double> sharedMerCounts = new ArrayList<Double>();
	private HashMap<String, Integer> skipMers = new HashMap<String, Integer>();

	private int totalTrials = 10000;
	private boolean halfError = false;

	private static Random generator = null;
	public static int seed = 0;

	public static void main(String[] args) throws Exception {
		boolean usage1 = true;
		if (args.length >= 5 && args.length <= 6) {
			usage1=false;
		} else if (args.length >= 7) {
			usage1 = true;
		} else {
			printUsage();
			System.exit(1);
		}

		KmerStatSimulator f = new KmerStatSimulator();
		
		f.totalTrials = Integer.parseInt(args[0]);
		if (usage1) {
			f.requestedLength = Double.parseDouble(args[2]);
			f.kmer = Integer.parseInt(args[1]);
			f.overlap = Integer.parseInt(args[3]);
			if (args.length > 7) {
				f.halfError  = Boolean.parseBoolean(args[7]);
			}
			if (args.length > 8) {
				f.reference = args[8];
			}
			if (f.overlap > f.requestedLength) {
				System.err.println("Cannot have overlap > sequence length");
				System.exit(1);
			}
			if (args.length > 9) {
				f.loadSkipMers(args[9]);
			}

			f.simulate(Double.parseDouble(args[4]), Double.parseDouble(args[5]),
					Double.parseDouble(args[6]));

		} else {
			f.requestedLength = Double.parseDouble(args[1]);
			if (args.length > 5) {
				f.reference = args[5];
			}

			f.simulate(Double.parseDouble(args[2]), Double.parseDouble(args[3]),
					Double.parseDouble(args[4]));			
		}
	}

	public static void printUsage() {
		System.err
				.println("Example usage: simulateSharedKmers <#trials> <kmer size> <seq length> <overlap length> <insertion> <del> <subst> [only one sequence error] [reference genome] [kmers to ignore]");
		System.err
		.println("Usage 2: simulateSharedKmers <#trials> <seq length> <insertion> <del> <subst> [reference genome]");
	}

	@SuppressWarnings("unused")
	public KmerStatSimulator() {
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
	
	private void loadSkipMers(String file) throws Exception {
		BufferedReader bf = Utils.getFile(file, "repeats");
		String line = null;

		while ((line = bf.readLine()) != null) {
			String[] split = line.trim().split("\\s+");
			String mer = split[0].trim();
			int count = Integer.parseInt(split[1]);
			this.skipMers.put(mer, count);
		}
		bf.close();
	}

	private String buildRandomSequence(int length) {
		StringBuilder st = new StringBuilder();

		for (int i = 0; i < length; i++) {
			st.append(getRandomBase(null));
		}
		return st.toString();
	}

	public double compareKmers(String first, String second) {
		HashSet<String> firstSeqs = new HashSet<String>(first.length());
		HashSet<String> totalSeqs = new HashSet<String>(first.length()+second.length());
		HashSet<String> shared = new HashSet<String>(first.length());

		for (int i = 0; i <= first.length() - this.kmer; i++) {
			String fmer = first.substring(i, i + this.kmer);
			if (!this.skipMers.containsKey(fmer)) { 
				firstSeqs.add(fmer);
			}
			totalSeqs.add(fmer);
		}

		for (int i = 0; i <= second.length() - this.kmer; i++) {
			String smer = second.substring(i, i + this.kmer);
			if (firstSeqs.contains(smer)) {
				shared.add(smer);
			} else {
				totalSeqs.add(smer);
			}
		}
		this.sharedCount = shared.size();
		return shared.size() / (double) totalSeqs.size();
	}
	
	public double compareMinHash(String first, String second) {
		MinHash h1 = new MinHash(new Sequence(first, new SequenceId(1)), this.kmer, 1256, null, null, true);
		MinHash h2 = new MinHash(new Sequence(second, new SequenceId(2)), this.kmer, 1256, null, null, true);
		
		return h1.jaccard(h2);
	}

	private char getRandomBase(Character toExclude) {
		Character result = null;

		while (result == null) {
			double base = generator.nextDouble();
			if (base < 0.25) {
				result = 'A';
			} else if (base < 0.5) {
				result = 'C';
			} else if (base < 0.75) {
				result = 'G';
			} else {
				result = 'T';
			}

			if (toExclude != null && toExclude.equals(result)) {
				result = null;
			}
		}

		return result;
	}

	@SuppressWarnings("unused")
	private String getSequence(int firstLen, int firstPos, String sequence,
			double errorRate, StringBuilder profile, StringBuilder realErrorStr) {
		return getSequence(firstLen, firstPos, sequence, errorRate, profile,
				realErrorStr, 0.792, 0.122, 0.086, true);
	}

	private String getSequence(int seqLength, int firstPos, String sequence,
			double errorRate, StringBuilder profile,
			StringBuilder realErrorStr, double insertionRate,
			double deletionRate, double substitutionRate, boolean trimRight) {
				
		StringBuilder firstSeq = new StringBuilder();
		firstSeq.append(sequence.substring(firstPos,
				Math.min(sequence.length(), firstPos + 2 * seqLength)));

		if (firstSeq.length() < 2 * seqLength) {
			firstSeq.append(sequence.substring(
					0,
					Math.min(sequence.length(),
							(2 * seqLength - firstSeq.length()))));
		}
		
		//use a linked list for insertions
		LinkedList<Character> modifiedSequence = new LinkedList<>();
		for (char a : firstSeq.toString().toCharArray())
			modifiedSequence.add(a);

		// now mutate
		int realError = 0;
		ListIterator<Character> iter = modifiedSequence.listIterator();
		while (iter.hasNext()) {
			char i = iter.next();

			if (generator.nextDouble() < errorRate) {
				double errorType = generator.nextDouble();
				if (errorType < substitutionRate) { // mismatch
													// switch base
					
					iter.set(getRandomBase(i));
					
					//firstSeq.setCharAt(i, getRandomBase(firstSeq.charAt(i)));
					//System.err.println("sub");
					realError++;
					i++;
				} else if (errorType < insertionRate + substitutionRate) { // insert
					
					iter.previous();
					iter.add(getRandomBase(null));
					//firstSeq.insert(i, getRandomBase(null));
					// profile.insert(i+1,"X");
					realError++;
					//i += 2;
				} else { // delete
					
					iter.remove();
					// firstSeq.setCharAt(i, 'D');
					// profile.setCharAt(i, '-');
					//System.err.println("delete");
					realError++;
				}
			} else {
				//i++;
			}
		}
		
		firstSeq = new StringBuilder();
		for (char c : modifiedSequence)
			firstSeq.append(c);
		
		realErrorStr.append((double) realError / seqLength);
		
		if (trimRight) {
			return firstSeq.substring(0, seqLength).toString();
		}
		
		return firstSeq.substring(firstSeq.length()-seqLength, firstSeq.length()).toString();
	}

	private void outputStats(ArrayList<Double> values, PrintStream out) {
		double mean = 0.0;
		double variance = 0.0;
		
		int N = 0;		
		for (double d : values) {
			N++;
			mean += d;
		}
		mean = mean/N;
		
		N = 0;
		for (double d : values) {
			N++;
			variance += (d-mean)*(d-mean);
		}
		
		variance /= (N-1);
		
		double stdev = Math.sqrt(variance);
		out.print(mean + "\t" + stdev);
	}

	public void simulate(double insertionRate, double delRate, double subRate)
			throws Exception {
		double errorRate = insertionRate + delRate + subRate;
		double insertionPercentage = insertionRate / errorRate;
		double deletionPercentage = delRate / errorRate;
		double subPercentage = subRate / errorRate;

		if (errorRate < 0 || errorRate > 1) {
			System.err.println("Error rate must be between 0 and 1");
			System.exit(1);
		}
		System.err.println("Started...");

		String[] sequences = null;
		if (this.reference != null) {
			FastaData data = new FastaData(this.reference, 0);
			data.enqueueFullFile();
			sequences = new String[data.getNumberProcessed()];
			int i = 0;
			while (!data.isEmpty())
				sequences[i++] = data.dequeue().getString().toUpperCase().replace("N", "");
		}
		System.err.println("Loaded reference");
		
		for (int i = 0; i < this.totalTrials; i++) {
			if (i % 100 == 0) {
				System.err.println("Done " + i + "/" + this.totalTrials);
			}
			int sequenceLength = (int) this.requestedLength;
			int firstPos = 0;

			String sequence = null;
			int seqID = 0;
			if (this.reference != null) {
				sequence = null;
				while (sequence == null
						|| sequence.length() < 4 * sequenceLength) {
					// pick a sequence from our reference
					seqID = generator.nextInt(sequences.length);
					sequence = sequences[seqID];
				}

				// now pick a position
				firstPos = generator.nextInt(sequence.length());
			} else {
				sequence = buildRandomSequence(sequenceLength * 4);
			}

			// simulate sequence with error
			StringBuilder firstAdj = new StringBuilder();
			StringBuilder errors = new StringBuilder();
			String firstSeq = getSequence(sequenceLength, firstPos, sequence,
					errorRate, firstAdj, errors, insertionPercentage,
					deletionPercentage, subPercentage, false);

			if (this.kmer < 0) { // we were only asked to simulate sequences not compare
				System.out.println(">s" + i + " " + seqID + " " + (firstPos+sequenceLength));
				System.out.println(Utils.convertToFasta(firstSeq));
				continue;
			}
			
			// compare number of shared kmers out of total to another sequence
			// from
			// same position
			int offset = (int) ((this.requestedLength * 2) - this.overlap);
			int secondPos = (firstPos + offset) % sequence.length();
			String secondSeq = getSequence(sequenceLength, secondPos, sequence,
					(this.halfError ? 0 : errorRate), firstAdj, errors, (this.halfError ? 0 :insertionPercentage),
					(this.halfError ? 0 : deletionPercentage), (this.halfError ? 0 : subPercentage), true);
			if (this.verbose) {
				System.err.println("Given seq " + firstPos + " of len " + sequence.length() + " and offset " + secondPos + " due to offset " + offset);
				System.err.println(">" + seqID + "_" + firstPos + "\n" + firstSeq);
				System.err.println(">" + seqID + "_" + secondPos + "\n" + secondSeq);
			}
			if (firstSeq.length() != secondSeq.length() || firstSeq.length() != this.requestedLength) {
				System.err.println("Error wrong length first: " + firstSeq.length() + " second: " + secondSeq.length() + " requested " + this.requestedLength);
				System.exit(1);
			}
			this.sharedJaccard.add(compareKmers(firstSeq, secondSeq));
			this.sharedMinHash.add(compareMinHash(firstSeq, secondSeq));
			this.sharedMerCounts.add(this.sharedCount);

			// compare number of shared kmers out of total to another sequence
			// from a
			// non-overlapping position
			// get a non-overlapping position
			if (this.reference != null) {
				sequence = null;
				int secondSeqID = 0;
				while (sequence == null
						|| sequence.length() < 2 * sequenceLength) {
					secondSeqID = generator.nextInt(sequences.length);
					sequence = sequences[secondSeqID];
				}
				secondPos = generator.nextInt(sequence.length());
				while (seqID == secondSeqID && Utils
						.getRangeOverlap(firstPos, firstPos + sequenceLength,
								secondPos, secondPos + sequenceLength) > 0) {
					secondPos = generator.nextInt(sequence.length());
				}
				// generate error for second sequence
				secondSeq = getSequence(sequenceLength, secondPos, sequence,
						(this.halfError ? 0 : errorRate), firstAdj, errors, (this.halfError ? 0 : insertionPercentage),
						(this.halfError ? 0 : deletionPercentage), (this.halfError ? 0 : subPercentage), true);
			} else {
				secondPos = 0;
				secondSeq = buildRandomSequence(sequenceLength);
			}

			if (firstSeq.length() != secondSeq.length() || firstSeq.length() != this.requestedLength) {
				System.err.println("Error wrong length " + firstSeq.length());
				System.exit(1);
			}
			// System.err.println("First: "+firstSeq.length());
			// System.err.println("Second: "+secondSeq.length());

			this.randomJaccard.add(compareKmers(firstSeq, secondSeq));
			this.randomMinHash.add(compareMinHash(firstSeq, secondSeq));
			this.randomMerCounts.add(this.sharedCount);
		}

		if (this.randomJaccard.size() != this.randomMerCounts.size()
				|| this.sharedJaccard.size() != this.sharedMerCounts.size()
				|| this.sharedJaccard.size() != this.randomJaccard.size()) {
			System.err.println("Error trial number not consistent!");
		}

		if (this.sharedMerCounts.size() == 0) { 
			return;
		}
		
		for (int i = 0; i < this.totalTrials; i++) {
			System.out.println(this.sharedMerCounts.get(i) + "\t"
					+ this.sharedJaccard.get(i) + "\t"
					+ this.sharedMinHash.get(i) + "\t"
					+ this.randomMerCounts.get(i) + "\t"
					+ this.randomJaccard.get(i) + "\t"
					+ this.randomMinHash.get(i));
		}
		System.out.print("Shared mer counts stats: ");
		outputStats(this.sharedMerCounts, System.out);
		System.out.println();

		System.out.print("Shared jaccard stats: ");
		outputStats(this.sharedJaccard, System.out);
		System.out.println();

		System.out.print("Shared MinHash jaccard stats: ");
		outputStats(this.sharedMinHash, System.out);
		System.out.println();

		System.out.print("Random mer counts stats: ");
		outputStats(this.randomMerCounts, System.out);
		System.out.println();

		System.out.print("Random jaccard stats: ");
		outputStats(this.randomJaccard, System.out);
		System.out.println();

		System.out.print("Random MinHash jaccard stats: ");
		outputStats(this.randomMinHash, System.out);
		System.out.println();
	}
}
