package com.secret.fastalign.main;
import java.util.ArrayList;
import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.MatchResult;
import com.secret.fastalign.minhash.MinHashSearch;

public class FastAlignMain 
{	
	private static final int DEFAULT_NUM_WORDS = 256;

	private static final int DEFAULT_KMER_SIZE = 14;

	private static final double DEFAULT_THRESHOLD = 0.07;
	
	private static final boolean DEFAULT_LARGE_MEMORY = true;

	public static void main(String[] args) throws Exception {
		String inFile = null;
		int kmerSize = DEFAULT_KMER_SIZE;
		double threshold = DEFAULT_THRESHOLD;
		int numWords = DEFAULT_NUM_WORDS; 
		boolean storeInMemory = DEFAULT_LARGE_MEMORY;

		for (int i = 0; i < args.length; i++) {
			if (args[i].trim().equalsIgnoreCase("-k")) {
				kmerSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("-s")) {
				inFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("--num-hashes")) {
				numWords = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--threshold")) {
				threshold = Double.parseDouble(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("-memory")) {
				storeInMemory = false;
			}
		}
		if (inFile == null) {
			printUsage("Error: no input fasta file specified");
		}

		System.err.println("Running with input fasta: " + inFile);
		System.err.println("kmer size:\t" + kmerSize);
		System.err.println("threshold:\t" + threshold);
		System.err.println("num hashed words:\t" + numWords);
		
		// read and index the kmers
		long startTime = System.nanoTime();

		FastaData data = new FastaData(inFile, kmerSize);
		
		System.err.println("Read in "+data.size()+" sequences.");

		System.err.println("Time (s) to read: " + (System.nanoTime() - startTime)*1.0e-9);
		
		//System.err.println("Press Enter");
		//System.in.read();
		
		MinHashSearch hashSearch = new MinHashSearch(kmerSize, numWords, data, storeInMemory, false);

		System.err.println("Time (s) to hash: " + (System.nanoTime() - startTime)*1.0e-9);

		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		startTime = System.nanoTime();

		ArrayList<MatchResult> results = hashSearch.findMatches(threshold);
		
		System.err.println("Time (s) to score and output: " + (System.nanoTime() - startTime)*1.0e-9);
		
		//sort to get the best scores on top
		//Collections.sort(results);		

		//System.err.println("Found "+results.size()+" matches:");
		
		//output result
		for (MatchResult match : results)
		{
			System.out.println(match);
		}		
	}

	public static void printUsage(String error) {
		if (error != null) {
			System.err.println(error);
		}
		System.err.println("Usage FastAlignMain <-s fasta file>");
		System.err.println("Options: ");
		System.err.println("\t -k [int merSize], default: " + DEFAULT_KMER_SIZE);
		System.err.println("\t -memory [store kmers in memory] default: " + DEFAULT_LARGE_MEMORY);
		System.err.println("\t  --num-hashes [int # hashes], default: " + DEFAULT_NUM_WORDS);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_THRESHOLD);
		System.exit(1);
	}
}