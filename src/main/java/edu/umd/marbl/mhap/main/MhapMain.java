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

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import edu.umd.marbl.mhap.general.FastaData;
import edu.umd.marbl.mhap.general.Sequence;
import edu.umd.marbl.mhap.general.SequenceId;
import edu.umd.marbl.mhap.sketch.CountMin;
import edu.umd.marbl.mhap.sketch.KmerCounts;
import edu.umd.marbl.mhap.sketch.MinHashSearch;
import edu.umd.marbl.mhap.sketch.SequenceSketchStreamer;
import edu.umd.marbl.mhap.utils.MhapRuntimeException;
import edu.umd.marbl.mhap.utils.PackageInfo;
import edu.umd.marbl.mhap.utils.ParseOptions;
import edu.umd.marbl.mhap.utils.Utils;

public final class MhapMain
{
	private final double acceptScore;
	private final HashSet<Long> filter;
	private final String inFile;
	private final int kmerSize;
	private final double maxShift;
	private final int minStoreLength;
	private final boolean noSelf;
	private final int numHashes;
	private final int numMinMatches;
	protected final int numThreads;
	private final String processFile;
	private final int subSequenceSize;	
	private final String toFile;
	private final boolean weighted;
	
	private final KmerCounts kmerCounter;

	private static final double DEFAULT_ACCEPT_SCORE = 0.04;

	private static final double DEFAULT_FILTER_CUTOFF = 1.0e-5;

	private static final int DEFAULT_KMER_SIZE = 16;

	private static final double DEFAULT_MAX_SHIFT_PERCENT = 0.2;

	private static final int DEFAULT_MIN_STORE_LENGTH = 0;

	private static final int DEFAULT_NUM_MIN_MATCHES = 3;
	
	private static final int DEFAULT_NUM_THREADS = Runtime.getRuntime().availableProcessors() * 2;
	private static final int DEFAULT_NUM_WORDS = 1024;
	private static final int DEFAULT_ORDERED_KMER_SIZE = 12;
	private static final int DEFAULT_SUB_SEQUENCE_SIZE = 100000;
	public static void main(String[] args) throws Exception
	{
		// set the locale
		Locale.setDefault(Locale.US);
		
		ParseOptions options = new ParseOptions();
		options.addStartTextLine("MHAP: MinHash Alignment Protocol. A tool for overlapping long-read sequences in bioinformatics.");
		options.addStartTextLine("\tVersion: "+PackageInfo.VERSION+", Build time: "+PackageInfo.BUILD_TIME);		
		options.addStartTextLine("\tUsage 1 (direct execution): java -server -Xmx<memory> -jar <MHAP jar> -s<fasta/dat from/self file> [-q<fasta/dat to file>] [-f<kmer filter list, must be sorted>]");
		options.addStartTextLine("\tUsage 2 (generate precomputed binaries): java -server -Xmx<memory> -jar <MHAP jar> -p<directory of fasta files> -q <output directory> [-f<kmer filter list, must be sorted>]");
		options.addOption("-s", "Usage 1 only. The FASTA or binary dat file (see Usage 2) of reads that will be stored in a box, and that all subsequent reads will be compared to.", "");
		options.addOption("-q", "Usage 1: The FASTA file of reads, or a directory of files, that will be compared to the set of reads in the box (see -s). Usage 2: The output directory for the binary formatted dat files.", "");
		options.addOption("-p", "Usage 2 only. The directory containing FASTA files that should be converted to binary format for storage.", "");
		options.addOption("-f", "k-mer filter file used for filtering out highly repetative k-mers. Must be sorted in descending order of frequency (second column).", "");
		options.addOption("-k", "[int], k-mer size used for MinHashing. The k-mer size for second stage filter is always set to "+DEFAULT_ORDERED_KMER_SIZE+", and currently cannot be modified.", DEFAULT_KMER_SIZE);
		options.addOption("--num-hashes", "[int], number of min-mers to be used in MinHashing.", DEFAULT_NUM_WORDS);
		options.addOption("--threshold", "[double], the threshold similarity score cutoff for the second stage sort-merge filter. This is based on the average number of k-mers matching in the overlapping region.", DEFAULT_ACCEPT_SCORE);
		options.addOption("--filter-threshold", "[double], the cutoff at which the k-mer in the k-mer filter file is considered repetitive. This value for a specific k-mer is specified in the second column in the filter file. If no filter file is provided, this option is ignored.", DEFAULT_FILTER_CUTOFF);
		options.addOption("--max-seq-size", "[int], not currently used.", DEFAULT_SUB_SEQUENCE_SIZE);
		options.addOption("--max-shift", "[double], region size to the left and right of the estimated overlap, as derived from the median shift and sequence length, where a k-mer matches are still considered valid. Second stage filter only.", DEFAULT_MAX_SHIFT_PERCENT);
		options.addOption("--num-min-matches", "[int], minimum # min-mer that must be shared before computing second stage filter. Any sequences below that value are considered non-overlapping.", DEFAULT_NUM_MIN_MATCHES);
		options.addOption("--num-threads", "[int], number of threads to use for computation. Typically set to 2 x #cores.", DEFAULT_NUM_THREADS);
		options.addOption("--weighted", "Perform weighted MinHashing.", true);
		options.addOption("--min-store-length", "[int], The minimum length of the read that is stored in the box. Used to filter out short reads from FASTA file.", DEFAULT_MIN_STORE_LENGTH);
		options.addOption("--no-self", "Do not compute the overlaps between sequences inside a box. Should be used when the to and from sequences are coming from different files.", false);
		options.addOption("--store-full-id", "Store full IDs as seen in FASTA file, rather than storing just the sequence position in the file. Some FASTA files have long IDS, slowing output of results. IDs not stored in compressed files.", false);
		options.addOption("--pacbio_fast", "Set all the parameters for the PacBio fast setting. This is the current best guidance, and could change at any time without warning.", false);
		options.addOption("--pacbio_sensitive", "Set all the parameters for the PacBio sensitive settings. This is the current best guidance, and could change at any time without warning.", false);
		options.addOption("--pacbio_experimental", "Set all the parameters for the PacBio experimental settings. This is the current best guidance, and could change at any time without warning.", false);
		
		if (!options.process(args))
			System.exit(0);
		
		//set the defaults for different type of data
		if (options.get("--pacbio_fast").getBoolean() || options.get("--pacbio_sensitive").getBoolean() || options.get("--pacbio_experimental").getBoolean())
		{
			if (!options.get("-k").isSet())
			{
				options.setOptions("-k", 16);
				if (options.get("--pacbio_experimental").getBoolean())
					options.setOptions("-k", 14);
			}
			
			if (!options.get("--num-min-matches").isSet())
			{
				options.setOptions("--num-min-matches", 3);
				if (options.get("--pacbio_experimental").getBoolean())
					options.setOptions("--num-min-matches", 1);
			}
			
			if (options.get("--pacbio_fast").getBoolean() && options.get("--pacbio_sensitive").getBoolean())
			{
				System.out.println("Two default sequence type parameters cannot be set at the same time.");
				System.out.println(options.helpMenuString());
				System.exit(1);
			}
			
			if (!options.get("--num-hashes").isSet())
			{
				if (options.get("--pacbio_fast").getBoolean())
					options.setOptions("--num-hashes", 512);
				else
				if (options.get("--pacbio_sensitive").getBoolean())
					options.setOptions("--num-hashes", 1256);
				else
				if (options.get("--pacbio_experimental").getBoolean())
					options.setOptions("--num-hashes", 1256);
			}
		}		
		
		if (options.get("-s").getString().isEmpty() && options.get("-p").getString().isEmpty())
		{
			System.out.println("Please set the -s or the -p options. See options below:");
			System.out.println(options.helpMenuString());
			System.exit(1);
		}
		
		if (!options.get("-p").getString().isEmpty() && options.get("-q").getString().isEmpty() )
		{
			System.out.println("Please set the -q option. See options below:");
			System.out.println(options.helpMenuString());
			System.exit(1);
		}
		
		//check for file existance
		if (!options.get("-p").getString().isEmpty() && !new File(options.get("-p").getString()).exists())
		{
			System.out.println("Could not find requested file/folder: "+options.get("-p").getString());
			System.exit(1);
		}

		//check for file existance
		if (!options.get("-s").getString().isEmpty() && !new File(options.get("-s").getString()).exists())
		{
			System.out.println("Could not find requested file/folder: "+options.get("-s").getString());
			System.exit(1);
		}
		
		//check for file existance
		if (!options.get("-q").getString().isEmpty() && !new File(options.get("-q").getString()).exists())
		{
			System.out.println("Could not find requested file/folder: "+options.get("-q").getString());
			System.exit(1);
		}
		
		//check for file existance
		if (!options.get("-f").getString().isEmpty() && !new File(options.get("-f").getString()).exists())
		{
			System.out.println("Could not find requested file/folder: "+options.get("-f").getString());
			System.exit(1);
		}
		
		//check range
		if (options.get("--num-threads").getInteger()<=0)
		{
			System.out.println("Number of threads must be positive.");
			System.exit(1);
		}

		//check range
		if (options.get("-k").getInteger()<=0)
		{
			System.out.println("k-mer size must be positive.");
			System.exit(1);
		}
		
		//check range
		if (options.get("--num-min-matches").getInteger()<=0)
		{
			System.out.println("Minimum number of matches must be positive.");
			System.exit(1);
		}

		//check range
		if (options.get("--min-store-length").getInteger()<0)
		{
			System.out.println("The minimum read length stored must be >=0.");
			System.exit(1);
		}

		//check range
		if (options.get("--max-shift").getDouble()<-1.0)
		{
			System.out.println("The minimum shift must be greater than -1.");
			System.exit(1);
		}
		
		//check range
		if (options.get("--threshold").getDouble()<0.0)
		{
			System.out.println("The second stage filter cutoff must be >=0.");
			System.exit(1);
		}

		//check other options
		//TODO move into the class
		if (options.get("--store-full-id").getBoolean())
			SequenceId.STORE_FULL_ID = true;
		else
			SequenceId.STORE_FULL_ID = false;

		
		//printing the options used
		System.err.println("Running with these settings:");
		System.err.println("Version = "+PackageInfo.VERSION);
		System.err.println("Build time = "+PackageInfo.BUILD_TIME);
		System.err.println(options);

		// start the main program
		MhapMain main = new MhapMain(options);

		//execute main computation code
		main.computeMain();
	}
	

	public MhapMain(ParseOptions options) throws IOException 
	{
		this.processFile = options.get("-p").getString();
		this.inFile = options.get("-s").getString();
		this.toFile = options.get("-q").getString();
		this.noSelf = options.get("--no-self").getBoolean();
		this.numThreads = options.get("--num-threads").getInteger();
		
		this.subSequenceSize = options.get("--max-seq-size").getInteger();
		this.numHashes = options.get("--num-hashes").getInteger();
		this.kmerSize = options.get("-k").getInteger();
		this.numMinMatches = options.get("--num-min-matches").getInteger();
		this.minStoreLength = options.get("--min-store-length").getInteger();
		this.maxShift = options.get("--max-shift").getDouble();
		this.acceptScore = options.get("--threshold").getDouble();
		this.weighted = options.get("--weighted").getBoolean();
	
		// read in the kmer filter set
		String filterFile = options.get("-f").getString();
		
		if (!filterFile.isEmpty())
		{
			long startTime = System.nanoTime();
			System.err.println("Reading in filter file " + filterFile + ".");
			try
			{
				this.filter = Utils.createKmerFilter(filterFile, options.get("--filter-threshold").getDouble(), this.kmerSize, 0);
			}
			catch (Exception e)
			{
				throw new MhapRuntimeException("Could not parse k-mer filter file.", e);
			}
			System.err.println("Time (s) to read filter file: " + (System.nanoTime() - startTime) * 1.0e-9);
			
			this.kmerCounter = null;
		}
		else
		{
			this.filter = null;
			this.kmerCounter = recordFastaKmerCounts(inFile, options.get("--filter-threshold").getDouble());
		}

	}
	
	public KmerCounts recordFastaKmerCounts(String file, double filterCutoff) throws IOException
	{
		System.err.println("Computing k-mer counts...");
		
		final FastaData data = new FastaData(this.inFile, 0);
		
		final CountMin<Long> countMin = new CountMin<>(1.0e-5, 1.0-1.0e-5, 0);
		//System.err.println(countMin.getDepth()+" "+countMin.getWidth());
		
		// figure out number of cores
		ExecutorService execSvc = Executors.newFixedThreadPool(this.numThreads);

		final AtomicInteger counter = new AtomicInteger();
		for (int iter = 0; iter < this.numThreads; iter++)
		{
			Runnable task = new Runnable()
			{
				@Override
				public void run()
				{
					try
					{
						Sequence seq = data.dequeue();
						while (seq != null)
						{
							//get the kmers integers
							long[] kmerHashes = Utils.computeSequenceHashesLong(seq.getString(), MhapMain.this.kmerSize, 0);
							
							//store the values
							for (long val : kmerHashes)
								countMin.add(val);

							//get the kmers integers for reverse compliment
							kmerHashes = Utils.computeSequenceHashesLong(seq.getReverseCompliment().getString(), MhapMain.this.kmerSize, 0);
							
							//store the values
							for (long val : kmerHashes)
								countMin.add(val);

							int currCount = counter.addAndGet(2);
							if (currCount % 5000 == 0)
								System.err.println("Kmers counted for " + currCount + " sequences (including reverse compliment)...");

							seq = data.dequeue();
						}
					}
					catch (IOException e)
					{
						throw new MhapRuntimeException(e);
					}
				}
			};

			// enqueue the task
			execSvc.execute(task);
		}

		// shutdown the service
		execSvc.shutdown();
		try
		{
			execSvc.awaitTermination(365L, TimeUnit.DAYS);
		}
		catch (InterruptedException e)
		{
			execSvc.shutdownNow();
			throw new MhapRuntimeException("Unable to finish all tasks.");
		}
		
		System.err.println("Computed k-mer counts for "+counter.get()+" sequences.");
		
		return new KmerCounts(countMin, counter.get(), filterCutoff);
	}

	public void computeMain() throws IOException
	{
		long startTotalTime = System.nanoTime();		
		long startTime = System.nanoTime();
		long processTime = System.nanoTime();
		
		//if processing a directory
		if (this.processFile!=null && !this.processFile.isEmpty())
		{
			System.err.println("Processing FASTA files for binary compression...");
			
			File file = new File(this.processFile);			
			if (!file.exists())
				throw new MhapRuntimeException("Process file does not exist.");

			if (this.toFile==null || this.toFile.isEmpty())
				throw new MhapRuntimeException("Target directory must be defined.");
			
			File toDirectory = new File(this.toFile);
			if (!toDirectory.exists() || !toDirectory.isDirectory())
				throw new MhapRuntimeException("Target directory doesn't exit.");
			
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
				
				SequenceSketchStreamer seqStreamer = getSequenceHashStreamer(pf.getAbsolutePath(), 0);
				
				String outputString = pf.getName();
				int i = outputString.lastIndexOf('.');
				if (i>0)
					outputString = outputString.substring(0, i);
				
				//combine with the directory name
				outputString = toDirectory.getPath()+File.separator+outputString+".dat";
				
				//store the file to disk
				seqStreamer.writeToBinary(outputString, false, this.numThreads);
				
				System.err.println("Processed "+seqStreamer.getNumberProcessed()+" sequences (fwd and rev).");
				System.err.println("Read, hashed, and stored file "+pf.getPath()+" to "+outputString+".");
				System.err.println("Time (s): " + (System.nanoTime() - startTime)*1.0e-9);
			}
			
			System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);

			return;
		}
		
		System.err.println("Processing files for storage in reverse index...");

		// read and index the kmers
		int seqNumberProcessed = 0;
				
		//create search object
		SequenceSketchStreamer seqStreamer = getSequenceHashStreamer(this.inFile, seqNumberProcessed);
		MinHashSearch hashSearch = getMatchSearch(seqStreamer);

		seqNumberProcessed += seqStreamer.getNumberProcessed()/2;
		System.err.println("Processed "+seqStreamer.getNumberProcessed()+" unique sequences (fwd and rev).");
		System.err.println("Time (s) to read and hash from file: " + (System.nanoTime() - processTime)*1.0e-9);

		long startTotalScoringTime = System.nanoTime();

		//System.err.println("Press Enter...");
		//System.in.read();
		
		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		if (this.toFile==null || this.toFile.isEmpty())
		{
			startTime = System.nanoTime();
			hashSearch.findMatches();
			System.err.println("Time (s) to score and output to self: " + (System.nanoTime() - startTime)*1.0e-9);
		}
		else
		{
			File file = new File(this.toFile);
			
			if (!file.exists())
				throw new MhapRuntimeException("To-file does not exist.");
			
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
			if (!this.noSelf)
			{
				hashSearch.findMatches();
				System.out.flush();
				System.err.println("Time (s) to score and output to self: " + (System.nanoTime() - startTime)*1.0e-9);
			}

			//no do to all files
			for (File cf : toFiles)
			{			
				// read and index the kmers
				seqStreamer = getSequenceHashStreamer(cf.getAbsolutePath(), seqNumberProcessed);
				System.err.println("Opened fasta file "+cf.getCanonicalPath()+".");
	
				//match the file
				startTime = System.nanoTime();
				hashSearch.findMatches(seqStreamer);
				
				//flush to get the output
				System.out.flush();
				
				seqNumberProcessed += seqStreamer.getNumberProcessed();
				System.err.println("Processed "+seqStreamer.getNumberProcessed()+" to sequences.");
				System.err.println("Time (s) to score, hash to-file, and output: " + (System.nanoTime() - startTime)*1.0e-9);
			}
		}
		
		//flush output
		System.out.flush();
		
		//output time
		System.err.println("Total scoring time (s): " + (System.nanoTime() - startTotalScoringTime)*1.0e-9);
		System.err.println("Total time (s): " + (System.nanoTime() - startTotalTime)*1.0e-9);
		
		//output final stats
		outputFinalStat(hashSearch);
	}

	public MinHashSearch getMatchSearch(SequenceSketchStreamer hashStreamer) throws IOException
	{
		return new MinHashSearch(hashStreamer, this.numHashes, this.numMinMatches, this.numThreads, false,
				this.minStoreLength, this.maxShift, this.acceptScore);
	}
	
	public SequenceSketchStreamer getSequenceHashStreamer(String file, int offset) throws IOException
	{
		SequenceSketchStreamer seqStreamer;
		if (file.endsWith(".dat"))
			seqStreamer = new SequenceSketchStreamer(file, offset);
		else
			seqStreamer = new SequenceSketchStreamer(file, this.kmerSize, this.numHashes, this.subSequenceSize,
					DEFAULT_ORDERED_KMER_SIZE, this.filter, this.kmerCounter, this.weighted, offset);

		return seqStreamer;
	}

	protected void outputFinalStat(MinHashSearch matchSearch)
	{
		System.err.println("MinHash search time (s): " + matchSearch.getMinHashSearchTime());
		//System.err.println("Sort-merge search time (s): " + matchSearch.getSortMergeTime());
		System.err.println("Total matches found: " + matchSearch.getMatchesProcessed());
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
