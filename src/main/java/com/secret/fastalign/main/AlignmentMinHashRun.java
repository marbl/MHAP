package com.secret.fastalign.main;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.LogManager;

import com.secret.fastalign.data.FastaData;
import com.secret.fastalign.data.Sequence;
import com.secret.fastalign.minhash.MinHash;

public class AlignmentMinHashRun {

	private static final int DEFAULT_NUM_HASHES = 5000;

	private static final int DEFAULT_KMER_SIZE = 10;

	private static final double DEFAULT_DATA_ERROR = 0.25;

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

		LogManager.getLogManager().reset();
		
		double kmerError = MinHash.probabilityKmerMatches(DEFAULT_DATA_ERROR, kmerSize);
		
		System.out.println("Probability of shared kmer in equal string: "+kmerError);
		
		// read and index the kmers
		long startTime = System.nanoTime();

		FastaData data = new FastaData(inFile, fastaSuffix, kmerSize);
		
		System.out.println("Read in "+data.size()+" sequences.");

		System.err.println("Time (s) to read: " + (System.nanoTime() - startTime)*1.0e-9);

		//System.out.println("Press Enter:");
		//System.in.read();

		// compute hashes
		startTime = System.nanoTime();

		MinHash minHash = new MinHash(numHashes, kmerSize);
		
		int skip = 0;
		for (Sequence seq : data.getSequences())
		{
			if (skip%5==0)
				minHash.addSequence(seq);
			
			skip++;
		}
		
		/*
		int[] s = new int[10000+kmerSize];
		Random generator = new Random(1);
		for (int iter=0; iter<s.length; iter++)
			s[iter] = generator.nextInt(3);

		ArrayList<Sequence> mySeq = new ArrayList<>();
		for (int iter=0; iter<10; iter++)
		{
			Sequence seq = new Sequence(Utils.errorString(s, DEFAULT_DATA_ERROR), new SequenceId(iter));
			mySeq.add(seq);
			minHash.addSequence(seq);
		}
		*/

		System.err.println("Time (s) to hash: " + (System.nanoTime() - startTime)*1.0e-9);

		// now that we have the hash constructed, go through all sequences to recompute their min and score their matches
		startTime = System.nanoTime();

		//find out the scores
		ArrayList<MatchResult> results = new ArrayList<MatchResult>();
		for (Sequence seq : data.getSequences())
		//for (Sequence seq : mySeq)
		{
			//get the matches
			List<MatchResult> matches = minHash.findMatches(seq, kmerError);
			
			//added to the list of solutions
			results.addAll(matches);
		}
		
		//sort to get the best scores on top
		//Collections.sort(results);		
		Collections.shuffle(results);
		
		System.err.println("Time (s) to score: " + (System.nanoTime() - startTime)*1.0e-9);
		
		System.out.println("Found "+results.size()+" matches:");
		
		Matrix matrix = MatrixLoader.load("/Users/kberlin/Dropbox/Projects/fast-align/src/test/resources/com/secret/fastalign/matrix/score_matrix.txt");
		
		//output result
		int count = 0;
		double mean = 0;
		for (MatchResult match : results)
		{
			Sequence s1 = data.getSequence(match.getFromId());
			Sequence s2 = data.getSequence(match.getToId());
			
			//compute the actual match
			double score = jaligner.SmithWatermanGotoh.align(new jaligner.Sequence(s1.getString()), new jaligner.Sequence(s2.getString()), matrix, 5, 3).getScore();
			score = score/Math.min((double)s1.length(), (double)s2.length());
						
			//System.out.format("Sequence match (%s - %s) with identity score %f (SW=%f).\n", match.getFromId(), match.getToId(), match.getScore(), score);
			System.out.format("%f %f %s %s %d\n", match.getScore(), score, match.getFromId(), match.getToId(), match.getFromShift());
			
			mean += match.getScore();
			
			count++;
			
			if (count>100)
				break;
		}
		
		mean = mean/count;
		
		System.out.println("Mean: "+mean);
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
