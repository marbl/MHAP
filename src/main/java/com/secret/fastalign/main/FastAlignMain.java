package com.secret.fastalign.main;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Locale;

import com.secret.fastalign.general.SequenceMinHashStreamer;
import com.secret.fastalign.minhash.MinHashSearch;
import com.secret.fastalign.utils.FastAlignRuntimeException;
import com.secret.fastalign.utils.Utils;

public final class FastAlignMain 
{	
	private static final double DEFAULT_FILTER_CUTOFF = 1.0e-5;

	private static final int DEFAULT_KMER_SIZE = 14;

	private static final boolean DEFAULT_LARGE_MEMORY = true;
			
	private static final int DEFAULT_MAX_SHIFT_ALLOWED = 800;

	private static final int DEFAULT_MIN_STORE_LENGTH = 0;
	
	private static final boolean DEFAULT_NO_SELF = false;
	
	private static final int DEFAULT_NUM_MIN_MATCHES = 4;

	private static final int DEFAULT_NUM_THREADS = Runtime.getRuntime().availableProcessors()*2;
	
	private static final int DEFAULT_NUM_WORDS = 1024;

	private static final int DEFAULT_SUB_SEQUENCE_SIZE = 5000;

	private static final double DEFAULT_THRESHOLD = 0.03;

	public static void main(String[] args) throws Exception 
	{
		//set the locale
	  Locale.setDefault(Locale.US);
		
		String inFile = null;
		String toFile = null;
		
		int kmerSize = DEFAULT_KMER_SIZE;
		double threshold = DEFAULT_THRESHOLD;
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

		for (int i = 0; i < args.length; i++) {
			if (args[i].trim().equalsIgnoreCase("-k")) {
				kmerSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("-s")) {
				inFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("-q")) {
				toFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("-p")) {
				processFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("-f")) {
				filterFile = args[++i];
			} else if (args[i].trim().equalsIgnoreCase("--num-hashes")) {
				numHashes = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--min-store-length")) {
				minStoreLength = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--threshold")) {
				threshold = Double.parseDouble(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--filter-threshold")) {
				filterThreshold = Double.parseDouble(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--num-min-matches")) {
				numMinMatches = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--subsequence-size")) {
				subSequenceSize = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--max-shift")) {
				maxShift = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--num-threads")) {
				numThreads = Integer.parseInt(args[++i]);
			} else if (args[i].trim().equalsIgnoreCase("--memory")) {
				storeInMemory = false;
			} else if (args[i].trim().equalsIgnoreCase("--no-self")) {
				noSelf = true;
			}
		}
		if (inFile == null && processFile==null) {
			printUsage("Error: no input or process files specified");
		}

		System.err.println("Running with input fasta: " + inFile);
		System.err.println("Running with process directory: " + processFile);
		System.err.println("Running with to directory or file: " + toFile);
		System.err.println("Running with kmer filter file: " + filterFile);
		System.err.println("kmer size:\t" + kmerSize);
		System.err.println("threshold:\t" + threshold);
		System.err.println("kmer filter percent cutoff:\t" + filterThreshold);
		System.err.println("num hashed words:\t" + numHashes);
		System.err.println("num min matches:\t" + numMinMatches);
		System.err.println("min hashed seq length:\t" + minStoreLength);
		System.err.println("subsequence size:\t" + subSequenceSize);
		System.err.println("max shift:\t" + maxShift);
		System.err.println("number of threads:\t" + numThreads);
		System.err.println("use large amount of memory:\t" + storeInMemory);
		System.err.println("compute alignment to self of -s file:\t" + !noSelf);
		
		//System.err.println("Press Enter");
		//System.in.read();
		
		long startTotalTime = System.nanoTime();		
		long startTime = System.nanoTime();

		//read in the kmer filter set
		HashSet<Integer> filter = null;
		if (filterFile!=null)
		{
			System.err.println("Reading in filter file "+filterFile+".");
			filter = Utils.createKmerFilter(filterFile, filterThreshold, kmerSize);
			System.err.println("Time (s) to read filter file: " + (System.nanoTime() - startTime)*1.0e-9);
		}

		long processTime = System.nanoTime();
		
		//if processing a directory
		if (processFile!=null)
		{
			File file = new File(processFile);			
			if (!file.exists())
				throw new FastAlignRuntimeException("Process file does not exist.");

			if (toFile==null)
				throw new FastAlignRuntimeException("Target directory option -q must be defined.");
			
			File toDirectory = new File(toFile);
			if (!toDirectory.exists() || !toDirectory.isDirectory())
				throw new FastAlignRuntimeException("Target directory doesn't exit.");
			
			//allocate directory files
			ArrayList<File> processFiles = new ArrayList<>();
			
			//if not dictory just add the file
			if (!file.isDirectory())
			{
				processFiles.add(file);
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
					processFiles.add(cf);				
			}
			
			for (File pf : processFiles)
			{
				startTime = System.nanoTime();
				
				SequenceMinHashStreamer seqStreamer = new SequenceMinHashStreamer(pf.getAbsolutePath(), kmerSize, numHashes, subSequenceSize, MinHashSearch.DEFAULT_SUB_KMER_SIZE, filter);
				
				String outputString = pf.getName();
				int i = outputString.lastIndexOf('.');
				if (i>0)
					outputString = outputString.substring(0, i);
				
				//combine with the directory name
				outputString = toDirectory.getPath()+File.separator+outputString+".dat";
				
				//store the file to disk
				seqStreamer.writeToBinary(outputString, false, numThreads);
				
				System.err.println("Read, hashed, and stored file "+pf.getPath()+" to "+outputString+".");
				System.err.println("Time (s): " + (System.nanoTime() - startTime)*1.0e-9);
			}
			
			System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);

			return;
		}
		
		// read and index the kmers
		SequenceMinHashStreamer seqStreamer;
		if (inFile.endsWith(".dat"))
		  seqStreamer = new SequenceMinHashStreamer(inFile);		
		else
	  seqStreamer = new SequenceMinHashStreamer(inFile, kmerSize, numHashes, subSequenceSize, MinHashSearch.DEFAULT_SUB_KMER_SIZE, filter);
		
		//create search object
		MinHashSearch hashSearch = new MinHashSearch(seqStreamer, numHashes, numMinMatches, numThreads, false, maxShift, minStoreLength);

		System.err.println("Processed "+seqStreamer.getNumberProcessed()+" unique sequences (fwd and rev).");
		System.err.println("Time (s) to read and hash from file: " + (System.nanoTime() - processTime)*1.0e-9);

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
				if (cf.getName().endsWith(".dat"))
					seqStreamer = new SequenceMinHashStreamer(cf.getCanonicalPath());
				else
					seqStreamer = new SequenceMinHashStreamer(cf.getCanonicalPath(), kmerSize, numHashes, subSequenceSize, MinHashSearch.DEFAULT_SUB_KMER_SIZE, filter);
				System.err.println("Opened fasta file "+cf.getCanonicalPath()+".");
	
				//match the file
				startTime = System.nanoTime();
				hashSearch.findMatches(seqStreamer, threshold);
				
				System.out.flush();
				System.err.println("Processed "+seqStreamer.getNumberProcessed()+" to sequences.");
				System.err.println("Time (s) to score, hash to-file, and output: " + (System.nanoTime() - startTime)*1.0e-9);
			}
		}
		
		//flush output
		System.out.flush();

		System.err.println("Total scoring time (s): " + (System.nanoTime() - startTotalScoringTime)*1.0e-9);
		System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);
		System.err.println("Total sequences searched: " + hashSearch.getNumberSequencesSearched());
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
		System.err.println("Usage 1 FastAlignMain -s<fasta from/self file> [-q<fasta to file] [-f<kmer filter list]");
		System.err.println("Usage 2 FastAlignMain -p<directory of fasta files> -q <output directory> [-f<kmer filter list]");
		System.err.println("Options: ");
		System.err.println("\t -k [int merSize], default: " + DEFAULT_KMER_SIZE);
		System.err.println("\t  --memory [do not store kmers in memory]");
		System.err.println("\t  --num-hashes [int # hashes], default: " + DEFAULT_NUM_WORDS);
		System.err.println("\t  --min-store-length [int # of minimum sequence length that is hashed], default: " + DEFAULT_MIN_STORE_LENGTH);
		System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_THRESHOLD);
		System.err.println("\t  --max-shift [int # max sequence shift allowed for a valid kmer relative to median value], default: " + DEFAULT_MAX_SHIFT_ALLOWED);
		System.err.println("\t  --num-min-matches [int # hashes that maches before performing local alignment], default: " + DEFAULT_NUM_MIN_MATCHES);
		System.err.println("\t  --num-threads [int # threads to use for computation], default (2 x #cores): " + DEFAULT_NUM_THREADS);
		System.err.println("\t  --subsequence-size [int size of maximum minhashed sequence], default: " + DEFAULT_SUB_SEQUENCE_SIZE);
		System.err.println("\t  --no-self [do not compute results to self], default: "+DEFAULT_NO_SELF);
		System.exit(1);
	}
}