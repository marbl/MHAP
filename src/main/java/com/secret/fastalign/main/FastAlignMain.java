package com.secret.fastalign.main;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Collections;

import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.minhash.MinHashSearch;
import com.secret.fastalign.utils.FastAlignRuntimeException;

public final class FastAlignMain 
{	
	private static final int DEFAULT_NUM_WORDS = 1024;

	private static final int DEFAULT_KMER_SIZE = 16;

	private static final double DEFAULT_THRESHOLD = 0.04;
			
	private static final int DEFAULT_NUM_MIN_MATCHES = 2;

	private static final int DEFAULT_SUB_SEQUENCE_SIZE = 5000;
	
	private static final int DEFAULT_NUM_THREADS = Runtime.getRuntime().availableProcessors()*2;
	
	private static final boolean DEFAULT_LARGE_MEMORY = true;

	private static final boolean DEFAULT_NO_SELF = false;

	public static void main(String[] args) throws Exception 
	{
		String inFile = null;
		String toFile = null;
		
		int kmerSize = DEFAULT_KMER_SIZE;
		double threshold = DEFAULT_THRESHOLD;
		int numWords = DEFAULT_NUM_WORDS;
		int numMinMatches = DEFAULT_NUM_MIN_MATCHES;
		int subSequenceSize = DEFAULT_SUB_SEQUENCE_SIZE; 
		boolean storeInMemory = DEFAULT_LARGE_MEMORY;
		int numThreads = DEFAULT_NUM_THREADS;
		boolean noSelf = DEFAULT_NO_SELF;

		for (int i = 0; i < args.length; i++) {
			if (args[i].trim().equalsIgnoreCase("-k")) {
				kmerSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("-s")) {
				inFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("-q")) {
				toFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("--num-hashes")) {
				numWords = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--threshold")) {
				threshold = Double.parseDouble(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--num-min-matches")) {
				numMinMatches = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--subsequence-size")) {
				subSequenceSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--num-threads")) {
				numThreads = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--memory")) {
				storeInMemory = false;
			} else if (args[i].trim().equalsIgnoreCase("--no-self")) {
				noSelf = true;
			}
		}
		if (inFile == null) {
			printUsage("Error: no input fasta file specified");
		}

		System.err.println("Running with input fasta: " + inFile);
		System.err.println("kmer size:\t" + kmerSize);
		System.err.println("threshold:\t" + threshold);
		System.err.println("num hashed words:\t" + numWords);
		System.err.println("num min matches:\t" + numMinMatches);
		System.err.println("subsequence size:\t" + subSequenceSize);
		System.err.println("number of threads:\t" + numThreads);
		System.err.println("use large amount of memory:\t" + storeInMemory);
		System.err.println("compute alignment to self of -s file:\t" + !noSelf);
		
		long startTotalTime = System.nanoTime();

		// read and index the kmers
		FastaData data = new FastaData(inFile);
		//System.err.println("Read in "+data.currentCacheSize()+" sequences.");
		
		//System.err.println("Press Enter");
		//System.in.read();
		
		long startTime = System.nanoTime();
		MinHashSearch hashSearch = new MinHashSearch(data, kmerSize, numWords, numMinMatches, subSequenceSize, numThreads, storeInMemory, false);
		System.err.println("Processed "+data.getNumberProcessed()+" sequences.");
		System.err.println("Time (s) to read and hash from file: " + (System.nanoTime() - startTime)*1.0e-9);

		long startTotalScoringTime = System.nanoTime();

		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		if (toFile==null)
		{
			startTime = System.nanoTime();
			hashSearch.findMatches(threshold);
			System.err.println("Time (s) to score and output to self: " + (System.nanoTime() - startTime)*1.0e-9);
		}
		else
		{
			File file = new File(toFile);
			
			if (!file.exists())
				throw new FastAlignRuntimeException("To-file does not exist.");
			
			ArrayList<File> toFiles = new ArrayList<>();
			
			//if not dictory just add the file
			if (!file.isDirectory())
			{
				toFiles.add(file);
			}
			else
			{			
				//read the directory content
				File[] fileList = file.listFiles(new FilenameFilter()
				{				
					@Override
					public boolean accept(File dir, String name)
					{
						if (!name.startsWith("."))
							return true;
						
						return false;
					}
				});
				
				for (File cf : fileList)
					toFiles.add(cf);
			}

			//sort the files in alphabetical order
			Collections.sort(toFiles);

			//first perform to self
			startTime = System.nanoTime();
			if (!noSelf)
			{
				hashSearch.findMatches(threshold);
				System.out.flush();
				System.err.println("Time (s) to score and output to self: " + (System.nanoTime() - startTime)*1.0e-9);
			}

			//no do to all files
			for (File cf : toFiles)
			{			
				// read and index the kmers
				data = new FastaData(cf.getCanonicalPath());
				System.err.println("Opened fasta file "+cf.getCanonicalPath()+".");
	
				//match the file
				startTime = System.nanoTime();
				hashSearch.findMatches(data, threshold);
				
				System.out.flush();
				System.err.println("Processed "+data.getNumberProcessed()+" to sequences.");
				System.err.println("Time (s) to score, hash to-file, and output: " + (System.nanoTime() - startTime)*1.0e-9);
			}
		}

		System.err.println("Total scoring time (s): " + (System.nanoTime() - startTotalScoringTime)*1.0e-9);
		System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);
		System.err.println("Total matches found: "+hashSearch.getMatchesProcessed());
		System.err.println("Average number of matches per lookup: " + (double)hashSearch.getMatchesProcessed()/(double)hashSearch.getNumberSequencesSearched()*100.0);
		System.err.println("Average % of hashed sequences hit per lookup: " + (double)hashSearch.getNumberSequencesHit()/(double)(hashSearch.size()*hashSearch.getNumberSequencesSearched())*100.0);
		System.err.println("Average % of hashed sequences hit that are matches: " + (double)hashSearch.getMatchesProcessed()/(double)hashSearch.getNumberSequencesHit()*100.0);
		System.err.println("Average % of hashed sequences fully compared that are matches: " + (double)hashSearch.getMatchesProcessed()/(double)hashSearch.getNumberSequencesFullyCompared()*100.0);
		System.err.flush();
	}

	public static void printUsage(String error) {
		if (error != null) {
			System.err.println(error);
		}
		System.err.println("Usage FastAlignMain -s<fasta from/self file> [-q<fasta to file]");
		System.err.println("Options: ");
		System.err.println("\t -k [int merSize], default: " + DEFAULT_KMER_SIZE);
		System.err.println("\t  --memory [do not store kmers in memory]");
		System.err.println("\t  --num-hashes [int # hashes], default: " + DEFAULT_NUM_WORDS);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_THRESHOLD);
		System.err.println("\t  --num-min-matches [int # hashes that maches before performing local alignment], default: " + DEFAULT_NUM_MIN_MATCHES);
		System.err.println("\t  --num-threads [int # threads to use for computation], default (2 x #cores): " + DEFAULT_NUM_THREADS);
		System.err.println("\t  --subsequence-size [int size of maximum minhashed sequence], default: " + DEFAULT_SUB_SEQUENCE_SIZE);
		System.err.println("\t  --no-self [do not compute results to self], default: "+DEFAULT_NO_SELF);
		System.exit(1);
	}
}