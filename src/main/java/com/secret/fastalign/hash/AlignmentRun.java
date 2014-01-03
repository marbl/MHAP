package com.secret.fastalign.hash;
import java.util.ArrayList;
import java.util.List;

public class AlignmentRun {

	private static final int DEFAULT_NUM_HASHES = 1000;

	private static final int DEFAULT_KMER_SIZE = 10;


	private static final int DEFAULT_SKIP = 500;
	private static final int DEFAULT_THRESHOLD = 1;
	private static final String[] fastaSuffix = {"fna", "contigs", "final", "fasta", "fa"};

	public static void main(String[] args) throws Exception {
		String inFile = null;
		int kmerSize = DEFAULT_KMER_SIZE;
		int threshold = DEFAULT_THRESHOLD;
		int numHashes = DEFAULT_NUM_HASHES; 
		int maxSkip = DEFAULT_SKIP;

		for (int i = 0; i < args.length; i++) {
			if (args[i].trim().equalsIgnoreCase("-k")) {
				kmerSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("-s")) {
				inFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("--num-hashes")) {
				numHashes = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--threshold")) {
				threshold = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--max-skip")) {
				maxSkip = Integer.parseInt(args[++i]);
			}
		}
		if (inFile == null) {
			printUsage("Error: no input fasta file specified");
		}

		System.err.println("Running with input fasta: " + inFile);
		System.err.println("kmer size:\t" + kmerSize);
		System.err.println("threshold:\t" + threshold);
		System.err.println("num hashes\t" + numHashes);
		System.err.println("max skip\t" + maxSkip);

		// read and index the kmers
		long startTime = System.nanoTime();

		FastaData data = new FastaData(inFile, fastaSuffix, kmerSize);

		System.err.println("Time (s) to read: " + (System.nanoTime() - startTime)*1.0e-9);

		// compute hashes
		startTime = System.nanoTime();

		MinHash minHash = new MinHash(numHashes);
		minHash.addData(data);

		System.err.println("Time (s) to hash: " + (System.nanoTime() - startTime)*1.0e-9);

		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		startTime = System.nanoTime();

		//find out the scores
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (Sequence seq : data.getSequences())
		{
			//get the matches
			List<MatchResult> matches = minHash.findMatches(seq, 0.0);
			
			//added to the list of solutions
			results.addAll(matches);
		}
		
		System.err.println("Time (s) to score: " + (System.nanoTime() - startTime)*1.0e-9);
		
		System.out.println("Found "+results.size()+" matches:");
		
		//output result
		for (MatchResult match : results)
		{
			System.out.format("Sequence match (%s - %s) with identity score %f.\n", match.getFromId(), match.getToId(), match.getScore());
		}
	}

	public static void printUsage(String error) {
		if (error != null) {
			System.err.println(error);
		}
		System.err.println("Usage buildMulti <-s fasta file>");
		System.err.println("Options: ");
		System.err.println("\t -k [int merSize], default: " + DEFAULT_KMER_SIZE);
		System.err.println("\t  --num-hashes [int # hashes], default: " + DEFAULT_NUM_HASHES);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_THRESHOLD);
		System.err.println("\t --max-skip [int bp maximum distance to nearest minimum value when guessing overlap positions], default: " + DEFAULT_SKIP);
		System.exit(1);
	}
}
