package com.secret.fastalign.main;

import java.io.IOException;
import java.util.HashSet;
import java.util.Locale;

import com.secret.fastalign.general.AbstractSequenceHashStreamer;
import com.secret.fastalign.general.AbstractSequenceSearchMain;
import com.secret.fastalign.minhash.MinHashSearch;
import com.secret.fastalign.minhash.SequenceMinHashStreamer;
import com.secret.fastalign.minhash.SequenceMinHashes;
import com.secret.fastalign.utils.Utils;

public final class FastAlignMain extends AbstractSequenceSearchMain<MinHashSearch, SequenceMinHashes>
{
	private final HashSet<Integer> filter;

	private final int kmerSize;

	private final int minStoreLength;

	private final int numHashes;

	private final int numMinMatches;

	private final int subSequenceSize;

	private final double acceptScore;

	private final int maxShift;

	private static final double DEFAULT_FILTER_CUTOFF = 1.0e-5;

	private static final int DEFAULT_KMER_SIZE = 16;

	private static final int DEFAULT_ORDERED_KMER_SIZE = 12;

	private static final boolean DEFAULT_LARGE_MEMORY = true;

	private static final int DEFAULT_MAX_SHIFT_ALLOWED = 800;

	private static final int DEFAULT_MIN_STORE_LENGTH = 0;

	private static final boolean DEFAULT_NO_SELF = false;

	private static final int DEFAULT_NUM_MIN_MATCHES = 3;

	private static final int DEFAULT_NUM_THREADS = Runtime.getRuntime().availableProcessors() * 2;

	private static final int DEFAULT_NUM_WORDS = 1024;

	private static final int DEFAULT_SUB_SEQUENCE_SIZE = 5000;

	private static final double DEFAULT_ACCEPT_SCORE = 0.03;

	public static void main(String[] args) throws Exception
	{
		// set the locale
		Locale.setDefault(Locale.US);

		String inFile = null;
		String toFile = null;

		int kmerSize = DEFAULT_KMER_SIZE;
		int numHashes = DEFAULT_NUM_WORDS;
		int numMinMatches = DEFAULT_NUM_MIN_MATCHES;
		int subSequenceSize = DEFAULT_SUB_SEQUENCE_SIZE;
		boolean storeInMemory = DEFAULT_LARGE_MEMORY;
		int numThreads = DEFAULT_NUM_THREADS;
		boolean noSelf = DEFAULT_NO_SELF;
		String filterFile = null;
		double filterThreshold = DEFAULT_FILTER_CUTOFF;
		int maxShift = DEFAULT_MAX_SHIFT_ALLOWED;
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
			else if (args[i].trim().equalsIgnoreCase("--subsequence-size"))
			{
				subSequenceSize = Integer.parseInt(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--threshold"))
			{
				acceptScore = Double.parseDouble(args[++i]);
			}
			else if (args[i].trim().equalsIgnoreCase("--max-shift"))
			{
				maxShift = Integer.parseInt(args[++i]);
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
		System.err.println("subsequence size:\t" + subSequenceSize);
		System.err.println("max shift:\t" + maxShift);
		System.err.println("threshold:\t" + acceptScore);
		System.err.println("number of threads:\t" + numThreads);
		System.err.println("use large amount of memory:\t" + storeInMemory);
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
		FastAlignMain main = new FastAlignMain(processFile, inFile, toFile, noSelf, subSequenceSize, numHashes, kmerSize,
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
				.println("Usage 1 FastAlignMain -s<fasta/dat from/self file> [-q<fasta/dat to file] [-f<kmer filter list]");
		System.err
				.println("Usage 2 FastAlignMain -p<directory of fasta files> -q <output directory> [-f<kmer filter list]");
		System.err.println("Options: ");
		System.err.println("\t -k [int merSize], default: " + DEFAULT_KMER_SIZE);
		System.err.println("\t  --memory [do not store kmers in memory]");
		System.err.println("\t  --num-hashes [int # hashes], default: " + DEFAULT_NUM_WORDS);
		System.err.println("\t  --min-store-length [int # of minimum sequence length that is hashed], default: "
				+ DEFAULT_MIN_STORE_LENGTH);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_ACCEPT_SCORE);
		System.err
				.println("\t  --max-shift [int # max sequence shift allowed for a valid kmer relative to median value], default: "
						+ DEFAULT_MAX_SHIFT_ALLOWED);
		System.err.println("\t  --num-min-matches [int # hashes that maches before performing local alignment], default: "
				+ DEFAULT_NUM_MIN_MATCHES);
		System.err.println("\t  --num-threads [int # threads to use for computation], default (2 x #cores): "
				+ DEFAULT_NUM_THREADS);
		System.err.println("\t  --subsequence-size [int size of maximum minhashed sequence], default: "
				+ DEFAULT_SUB_SEQUENCE_SIZE);
		System.err.println("\t  --no-self [do not compute results to self], default: " + DEFAULT_NO_SELF);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_ACCEPT_SCORE);
		System.err
				.println("\t  --max-shift [int # max sequence shift allowed for a valid kmer relative to median value], default: "
						+ DEFAULT_MAX_SHIFT_ALLOWED);
		System.exit(1);
	}

	public FastAlignMain(String processFile, String inFile, String toFile, boolean noSelf, int subSequenceSize,
			int numHashes, int kmerSize, int numMinMatches, int numThreads, HashSet<Integer> filter, int minStoreLength,
			int maxShift, double acceptScore)
	{
		super(processFile, inFile, toFile, noSelf, numThreads);
		this.subSequenceSize = subSequenceSize;
		this.numHashes = numHashes;
		this.kmerSize = kmerSize;
		this.numMinMatches = numMinMatches;
		this.minStoreLength = minStoreLength;
		this.filter = filter;
		this.maxShift = maxShift;
		this.acceptScore = acceptScore;
	}

	@Override
	public MinHashSearch getMatchSearch(AbstractSequenceHashStreamer<SequenceMinHashes> hashStreamer) throws IOException
	{
		return new MinHashSearch(hashStreamer, this.numHashes, this.numMinMatches, this.numThreads, false,
				this.minStoreLength, this.maxShift, this.acceptScore);
	}

	@Override
	public SequenceMinHashStreamer getSequenceHashStreamer(String file, int offset) throws IOException
	{
		SequenceMinHashStreamer seqStreamer;
		if (file.endsWith(".dat"))
			seqStreamer = new SequenceMinHashStreamer(file, offset);
		else
			seqStreamer = new SequenceMinHashStreamer(file, this.kmerSize, this.numHashes, this.subSequenceSize,
					DEFAULT_ORDERED_KMER_SIZE, this.filter, offset);

		return seqStreamer;
	}

	@Override
	protected void outputFinalStat(MinHashSearch matchSearch)
	{
		System.err.println("Total matches found: " + matchSearch.getMatchesProcessed());
		System.err.println("Average number of matches per lookup: " + (double) matchSearch.getMatchesProcessed()
				/ (double) matchSearch.getNumberSequencesSearched() * 100.0);
		System.err.println("Average % of hashed sequences hit per lookup: " + (double) matchSearch.getNumberSequencesHit()
				/ (double) (matchSearch.size() * matchSearch.getNumberSequencesSearched()) * 100.0);
		System.err.println("Average % of hashed sequences hit that are returned: "
				+ (double) matchSearch.getMatchesProcessed() / (double) matchSearch.getNumberSequencesHit() * 100.0);
		System.err.println("Average % of hashed sequences fully compared that are matches: " 
				+ (double)matchSearch.getMatchesProcessed()/(double)matchSearch.getNumberSequencesFullyCompared()*100.0);
		System.err.flush();
	}
}