package com.secret.fastalign.main;

import java.util.ArrayList;
import java.util.Random;
import java.util.GregorianCalendar;
import java.util.Calendar;
import java.util.HashSet;
import java.io.PrintStream;
import com.secret.fastalign.general.FastaData;
import com.secret.fastalign.general.SequenceId;
import com.secret.fastalign.utils.Utils;

public class KmerStatSimulator
{
	private int kmer = 12;
	private int overlap = 100;

	private ArrayList<Double> randomJaccard = new ArrayList<Double>();
	private ArrayList<Double> randomMerCounts = new ArrayList<Double>();
	private String reference = null;
	private double requestedLength = 5000;

	private ArrayList<String> sequences = new ArrayList<String>();
	private double sharedCount = 0;
	private ArrayList<Double> sharedJaccard = new ArrayList<Double>();
	private ArrayList<Double> sharedMerCounts = new ArrayList<Double>();

	private int totalTrials = 100;

	private static Random generator = null;
	public static int seed = 0;

	public static void main(String[] args) throws Exception
	{
		if (args.length < 3)
		{
			printUsage();
			System.exit(1);
		}

		KmerStatSimulator f = new KmerStatSimulator();
		f.totalTrials = Integer.parseInt(args[0]);
		f.requestedLength = Double.parseDouble(args[2]);
		f.kmer = Integer.parseInt(args[1]);
		f.overlap = Integer.parseInt(args[3]);
		if (args.length > 7)
		{
			f.reference = args[7];
		}
		if (f.overlap > f.requestedLength)
		{
			System.err.println("Cannot have overlap > sequence length");
			System.exit(1);
		}
		
		f.simulate(Double.parseDouble(args[4]), Double.parseDouble(args[5]), Double.parseDouble(args[6]));
	}

	public static void printUsage()
	{
		System.err
				.println("Example usage: simulateSharedKmers <#trials> <kmer size> <seq length> <overlap length> <insertion> <subst> <del> [refernce genome]");
	}

	@SuppressWarnings("unused")
	public KmerStatSimulator()
	{
		if (false)
		{
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

	private String buildRandomSequence(int length)
	{
		StringBuilder st = new StringBuilder();

		for (int i = 0; i < length; i++)
		{
			st.append(getRandomBase(null));
		}
		return st.toString();
	}

	public double compareKmers(String first, String second)
	{
		HashSet<String> firstSeqs = new HashSet<String>();
		HashSet<String> totalSeqs = new HashSet<String>();
		HashSet<String> shared = new HashSet<String>();

		for (int i = 0; i <= first.length() - this.kmer; i++)
		{
			String fmer = first.substring(i, i + this.kmer);
			firstSeqs.add(fmer);
			totalSeqs.add(fmer);
		}

		for (int i = 0; i <= second.length() - this.kmer; i++)
		{
			String smer = second.substring(i, i + this.kmer);

			if (firstSeqs.contains(smer))
			{
				shared.add(smer);
			}
			else
			{
				totalSeqs.add(smer);
			}
		}
		this.sharedCount = shared.size();
		return shared.size() / (double) totalSeqs.size();
	}

	private char getRandomBase(Character toExclude)
	{
		Character result = null;

		while (result == null)
		{
			double base = generator.nextDouble();
			if (base < 0.25)
			{
				result = 'A';
			}
			else if (base < 0.5)
			{
				result = 'C';
			}
			else if (base < 0.75)
			{
				result = 'G';
			}
			else
			{
				result = 'T';
			}

			if (toExclude != null && toExclude == result)
			{
				result = null;
			}
		}

		return result;
	}

	@SuppressWarnings("unused")
	private String getSequence(int firstLen, int firstPos, String sequence, double errorRate, StringBuilder profile,
			StringBuilder realErrorStr)
	{
		return getSequence(firstLen, firstPos, sequence, errorRate, profile, realErrorStr, 0.792, 0.122, 0.086);
	}

	private String getSequence(int firstLen, int firstPos, String sequence, double errorRate, StringBuilder profile,
			StringBuilder realErrorStr, double insertionRate, double deletionRate, double substitutionRate)
	{
		StringBuilder firstSeq = new StringBuilder();
		firstSeq.append(sequence.substring(firstPos, Math.min(sequence.length(), (firstPos + firstLen))));
		if (firstSeq.length() < firstLen)
		{
			firstSeq.append(sequence.substring(0, Math.min(sequence.length(), firstLen - firstSeq.length())));
		}

		// now mutate
		int realError = 0;
		for (int i = 0; i < firstSeq.length();)
		{
			if (generator.nextDouble() < errorRate)
			{
				double errorType = generator.nextDouble();
				if (errorType < substitutionRate)
				{ // mismatch
					// switch base
					firstSeq.setCharAt(i, getRandomBase(firstSeq.charAt(i)));
					realError++;
					i++;
				}
				else if (errorType < insertionRate + substitutionRate)
				{ // insert
					firstSeq.insert(i + 1, Character.toLowerCase(getRandomBase(null)));
					// profile.insert(i+1,"X");
					realError++;
					i += 2;
				}
				else
				{ // delete
					firstSeq.deleteCharAt(i);
					//firstSeq.setCharAt(i, 'D');
					// profile.setCharAt(i, '-');
					realError++;
				}
			}
			else
			{
				i++;
			}
		}
		realErrorStr.append((double) realError / firstLen);
		return firstSeq.toString();
	}

	private void outputStats(ArrayList<Double> values, PrintStream out)
	{
		double mean = 0.0;
		double variance = 0.0;
		int N = 0;
		for (double d : values)
		{
			N++;
			double delta = (d - mean);
			mean += (delta / N);
			variance += delta * (d - mean);
		}
		variance /= (double)N;
		double stdev = Math.sqrt(variance);
		out.print(mean + "\t" + stdev);
	}

	public void simulate(double insertionRate, double delRate, double subRate) throws Exception
	{
		double errorRate = insertionRate + delRate + subRate;
		double insertionPercentage = insertionRate / errorRate;
		double deletionPercentage = delRate / errorRate;
		double subPercentage = subRate / errorRate;

		if (errorRate < 0 || errorRate > 1)
		{
			System.err.println("Error rate must be between 0 and 1");
			System.exit(1);
		}
		System.err.println("Started...");

		FastaData data = null;
		if (this.reference != null)
		{
			data = new FastaData(this.reference);
		}
		System.err.println("Loaded reference");

		for (int i = 0; i < this.totalTrials; i++)
		{
			if (i % 1000 == 0)
			{
				System.err.println("Done " + i + "/" + this.totalTrials);
			}
			int firstLen = (int) this.requestedLength;
			int firstPos = 0;
			String sequence = buildRandomSequence(firstLen * 2);

			if (this.reference != null)
			{
				sequence = null;
				while (sequence == null || sequence.length() < 2 * firstLen)
				{
					// pick a sequence from our reference
					int seqID = generator.nextInt(this.sequences.size())+1;
					data.getSequence(new SequenceId(seqID, true));
				}

				// now pick a position
				firstPos = generator.nextInt(sequence.length());
			}
			// simulate sequence with error
			StringBuilder firstAdj = new StringBuilder();
			StringBuilder errors = new StringBuilder();
			String firstSeq = getSequence(firstLen, firstPos, sequence, errorRate, firstAdj, errors, insertionPercentage,
					deletionPercentage, subPercentage);

			// compare number of shared kmers out of total to another sequence from
			// same position
			int offset = (int) (firstLen - this.overlap);
			int secondPos = (firstPos + offset) % sequence.length();
			String secondSeq = getSequence(firstLen, secondPos, sequence, errorRate, firstAdj, errors, insertionPercentage,
					deletionPercentage, subPercentage);
			this.sharedJaccard.add(compareKmers(firstSeq, secondSeq));
			this.sharedMerCounts.add(this.sharedCount);

			// compare number of shared kmers out of total to another sequence from a
			// non-overlapping position
			// get a non-overlapping position
			if (this.reference != null)
			{
				secondPos = generator.nextInt(sequence.length());
				while (Utils.getRangeOverlap(firstPos, firstPos + firstLen, secondPos, secondPos + firstLen) > 0)
				{
					secondPos = generator.nextInt(sequence.length());
				}
			}
			else
			{
				secondPos = 0;
				sequence = buildRandomSequence(firstLen);
			}
			
			//generate error for second sequence
			secondSeq = getSequence(firstLen, secondPos, sequence, errorRate, firstAdj, errors, insertionPercentage,
					deletionPercentage, subPercentage);
			
			System.err.println("Hi "+firstSeq.length());
			System.err.println("Hi2 "+secondSeq.length());
				
			this.randomJaccard.add(compareKmers(firstSeq, secondSeq));
			this.randomMerCounts.add(this.sharedCount);
		}

		if (this.randomJaccard.size() != this.randomMerCounts.size() || this.sharedJaccard.size() != this.sharedMerCounts.size()
				|| this.sharedJaccard.size() != this.randomJaccard.size())
		{
			System.err.println("Error trial number not consistent!");
		}

		for (int i = 0; i < this.totalTrials; i++)
		{
			System.out.println(this.sharedMerCounts.get(i) + "\t" + this.sharedJaccard.get(i) + "\t" + this.randomMerCounts.get(i) + "\t"
					+ this.randomJaccard.get(i));
		}
		System.out.print("Shared mer counts stats: ");
		outputStats(this.sharedMerCounts, System.out);
		System.out.println();

		System.out.print("Shared jaccard stats: ");
		outputStats(this.sharedJaccard, System.out);
		System.out.println();

		System.out.print("Random mer counts stats: ");
		outputStats(this.randomMerCounts, System.out);
		System.out.println();

		System.out.print("Random jaccard stats: ");
		outputStats(this.randomJaccard, System.out);
		System.out.println();
	}
}
