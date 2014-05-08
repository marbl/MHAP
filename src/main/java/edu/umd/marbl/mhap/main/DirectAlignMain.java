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

import java.io.IOException;
import java.util.HashSet;
import java.util.Locale;

import edu.umd.marbl.mhap.direct.DirectHashSearch;
import edu.umd.marbl.mhap.direct.SequenceDirectHashStreamer;
import edu.umd.marbl.mhap.direct.SequenceDirectHashes;
import edu.umd.marbl.mhap.general.AbstractSequenceHashStreamer;
import edu.umd.marbl.mhap.general.AbstractSequenceSearchMain;
import edu.umd.marbl.mhap.utils.Utils;

public final class DirectAlignMain extends AbstractSequenceSearchMain<DirectHashSearch, SequenceDirectHashes>
{
	private final HashSet<Integer> filter;

	private final int kmerSize;

	private final int minStoreLength;

	private final int numHashes;

	private final int numMinMatches;

	private final double acceptScore;

	private final double maxShift;

	private static final double DEFAULT_FILTER_CUTOFF = 1.0e-5;

	private static final int DEFAULT_ORDERED_KMER_SIZE = 12;

	private static final int DEFAULT_KMER_SIZE = 16;

	private static final double DEFAULT_MAX_SHIFT_PERCENT = 0.2;

	private static final int DEFAULT_MIN_STORE_LENGTH = 0;

	private static final boolean DEFAULT_NO_SELF = false;

	private static final int DEFAULT_NUM_MIN_MATCHES = 10;

	private static final int DEFAULT_NUM_THREADS = Runtime.getRuntime().availableProcessors() * 2;

	private static final int DEFAULT_NUM_WORDS = 0;

	private static final double DEFAULT_ACCEPT_SCORE = 0.05;

	public static void main(String[] args) throws Exception
	{
		// set the locale
		Locale.setDefault(Locale.US);

		String inFile = null;
		String toFile = null;

		int kmerSize = DEFAULT_KMER_SIZE;
		int numHashes = DEFAULT_NUM_WORDS;
		int numMinMatches = DEFAULT_NUM_MIN_MATCHES;
		int numThreads = DEFAULT_NUM_THREADS;
		boolean noSelf = DEFAULT_NO_SELF;
		String filterFile = null;
		double filterThreshold = DEFAULT_FILTER_CUTOFF;
		double maxShift = DEFAULT_MAX_SHIFT_PERCENT;
		int minStoreLength = DEFAULT_MIN_STORE_LENGTH;
		String processFile = null;
		double acceptScore = DEFAULT_ACCEPT_SCORE;

		for (int i = 0; i < args.length; i++)
		{
			if (args[i].trim().equalsIgnoreCase("-k"))
			{
				kmerSize = Integer.parseInt(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("-s"))
			{
				inFile = args[++i];
			}
			else if (args[i].trim().equalsIgnoreCase("-q"))
			{
				toFile = args[++i];
			}
			else if (args[i].trim().equalsIgnoreCase("-p"))
			{
				processFile = args[++i];
			}
			else if (args[i].trim().equalsIgnoreCase("-f"))
			{
				filterFile = args[++i];
			}
			else if (args[i].trim().equalsIgnoreCase("--num-hashes"))
			{
				numHashes = Integer.parseInt(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--min-store-length"))
			{
				minStoreLength = Integer.parseInt(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--filter-threshold"))
			{
				filterThreshold = Double.parseDouble(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--num-min-matches"))
			{
				numMinMatches = Integer.parseInt(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--threshold"))
			{
				acceptScore = Double.parseDouble(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--max-shift"))
			{
				maxShift = Double.parseDouble(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--num-threads"))
			{
				numThreads = Integer.parseInt(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--no-self"))
			{
				noSelf = true;
			}
		}
		if (inFile == null && processFile == null)
		{
			printUsage("Error: no input or process files specified");
		}

		System.err.println("Running with input fasta: " + inFile);
		System.err.println("Running with process directory: " + processFile);
		System.err.println("Running with to directory or file: " + toFile);
		System.err.println("Running with kmer filter file: " + filterFile);
		System.err.println("kmer size:\t" + kmerSize);
		System.err.println("kmer filter percent cutoff:\t" + filterThreshold);
		System.err.println("num hashed words:\t" + numHashes);
		System.err.println("num min matches:\t" + numMinMatches);
		System.err.println("min hashed seq length:\t" + minStoreLength);
		System.err.println("max shift:\t" + maxShift);
		System.err.println("threshold:\t" + acceptScore);
		System.err.println("number of threads:\t" + numThreads);
		System.err.println("compute alignment to self of -s file:\t" + !noSelf);

		long startTime = System.nanoTime();

		// read in the kmer filter set
		HashSet<Integer> filter = null;
		if (filterFile != null)
		{
			System.err.println("Reading in filter file " + filterFile + ".");
			filter = Utils.createKmerFilter(filterFile, filterThreshold, kmerSize);
			System.err.println("Time (s) to read filter file: " + (System.nanoTime() - startTime) * 1.0e-9);
		}

		// start the main program
		DirectAlignMain main = new DirectAlignMain(processFile, inFile, toFile, noSelf, numHashes, kmerSize,
				numMinMatches, numThreads, filter, minStoreLength, maxShift, acceptScore);

		main.computeMain();
	}

	public static void printUsage(String error)
	{
		if (error != null)
		{
			System.err.println(error);
		}
		System.err
				.println("Usage 1 DirectAlignMain -s<fasta/dat from/self file> [-q<fasta/dat to file>] [-f<kmer filter list, must be sorted>]");
		System.err
				.println("Usage 2 DirectAlignMain -p<directory of fasta files> -q <output directory> [-f<kmer filter list, must be sorted>]");
		System.err.println("Options: ");
		System.err.println("\t -k [int merSize], default: " + DEFAULT_KMER_SIZE);
		System.err.println("\t  --memory [do not store kmers in memory]");
		System.err.println("\t  --num-hashes [int # random k-mers to select, 0 means take all (deterministic)], default: " + DEFAULT_NUM_WORDS);
		System.err.println("\t  --min-store-length [int # of minimum sequence length that is hashed], default: "
				+ DEFAULT_MIN_STORE_LENGTH);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_ACCEPT_SCORE);
		System.err
		.println("\t  --max-shift [double fraction of the overlap size where shift in k-mer match is still considered valid], default: "
				+ DEFAULT_MAX_SHIFT_PERCENT);
		System.err.println("\t  --num-min-matches [int # hashes that maches before performing local alignment], default: "
				+ DEFAULT_NUM_MIN_MATCHES);
		System.err.println("\t  --num-threads [int # threads to use for computation], default (2 x #cores): "
				+ DEFAULT_NUM_THREADS);
		System.err.println("\t  --no-self [do not compute results to self], default: " + DEFAULT_NO_SELF);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_ACCEPT_SCORE);
		System.err
				.println("\t  --max-shift [int # max sequence shift allowed for a valid kmer relative to median value], default: "
						+ DEFAULT_MAX_SHIFT_PERCENT);
		System.exit(1);
	}

	public DirectAlignMain(String processFile, String inFile, String toFile, boolean noSelf,
			int numHashes, int kmerSize, int numMinMatches, int numThreads, HashSet<Integer> filter, int minStoreLength,
			double maxShift, double acceptScore)
	{
		super(processFile, inFile, toFile, noSelf, numThreads);
		this.numHashes = numHashes;
		this.kmerSize = kmerSize;
		this.numMinMatches = numMinMatches;
		this.minStoreLength = minStoreLength;
		this.filter = filter;
		this.maxShift = maxShift;
		this.acceptScore = acceptScore;
	}

	@Override
	public DirectHashSearch getMatchSearch(AbstractSequenceHashStreamer<SequenceDirectHashes> hashStreamer) throws IOException
	{
		return new DirectHashSearch(hashStreamer, this.numHashes, this.numMinMatches, this.numThreads, false,
				this.minStoreLength, this.maxShift, this.acceptScore);
	}

	@Override
	public SequenceDirectHashStreamer getSequenceHashStreamer(String file, int offset) throws IOException
	{
		SequenceDirectHashStreamer seqStreamer;
		if (file.endsWith(".dat"))
			seqStreamer = new SequenceDirectHashStreamer(file, offset);
		else
			seqStreamer = new SequenceDirectHashStreamer(file, this.kmerSize, DEFAULT_ORDERED_KMER_SIZE, this.filter, offset);

		return seqStreamer;
	}

	@Override
	protected void outputFinalStat(DirectHashSearch matchSearch)
	{
		System.err.println("Total matches found: " + matchSearch.getMatchesProcessed());
		System.err.println("Hash table Shannon normalized entropy: " + matchSearch.hashTableNormalizedEnthropy());
		System.err.println("Average number of matches per lookup: " + (double) matchSearch.getMatchesProcessed()
				/ (double) matchSearch.getNumberSequencesSearched());
		System.err.println("Average number of table elements processed per lookup: " + (double) matchSearch.getNumberElementsProcessed()
				/ (double) (matchSearch.getNumberSequencesSearched()));
		System.err.println("Average number of table elements processed per match: " + (double) matchSearch.getNumberElementsProcessed()
				/ (double) (matchSearch.getMatchesProcessed()));
		System.err.println("Average % of hashed sequences hit per lookup: " + (double) matchSearch.getNumberSequencesHit()
				/ (double) (matchSearch.size() * matchSearch.getNumberSequencesSearched()) * 100.0);
		System.err.println("Average % of hashed sequences hit that are matches: "
				+ (double) matchSearch.getMatchesProcessed() / (double) matchSearch.getNumberSequencesHit() * 100.0);
		System.err.println("Average % of hashed sequences fully compared that are matches: " 
				+ (double)matchSearch.getMatchesProcessed()/(double)matchSearch.getNumberSequencesFullyCompared()*100.0);
		System.err.flush();
	}
}