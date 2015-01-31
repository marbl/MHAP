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
package edu.umd.marbl.mhap.utils;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.io.File;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.nio.ByteBuffer;

import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;

public final class Utils
{

	public enum ToProtein
	{
		AAA("K"), AAC("N"), AAG("K"), AAT("N"), ACA("T"), ACC("T"), ACG("T"), ACT("T"), AGA("R"), AGC("S"), AGG("R"), AGT(
				"S"), ATA("I"), ATC("I"), ATG("M"), ATT("I"), CAA("Q"), CAC("H"), CAG("Q"), CAT("H"), CCA("P"), CCC("P"), CCG(
				"P"), CCT("P"), CGA("R"), CGC("R"), CGG("R"), CGT("R"), CTA("L"), CTC("L"), CTG("L"), CTT("L"), GAA("E"), GAC(
				"D"), GAG("E"), GAT("D"), GCA("A"), GCC("A"), GCG("A"), GCT("A"), GGA("G"), GGC("G"), GGG("G"), GGT("G"), GTA(
				"V"), GTC("V"), GTG("V"), GTT("V"), TAA("X"), TAC("Y"), TAG("X"), TAT("Y"), TCA("S"), TCC("S"), TCG("S"), TCT(
				"S"), TGA("X"), TGC("C"), TGG("W"), TGT("C"), TTA("L"), TTC("F"), TTG("L"), TTT("F");

		/*
		 * Ala/A GCU, GCC, GCA, GCG Leu/L UUA, UUG, CUU, CUC, CUA, CUG Arg/R
		 * CGU, CGC, CGA, CGG, AGA, AGG Lys/K AAA, AAG Asn/N AAU, AAC Met/M AUG
		 * Asp/D GAU, GAC Phe/F UUU, UUC Cys/C UGU, UGC Pro/P CCU, CCC, CCA, CCG
		 * Gln/Q CAA, CAG Ser/S UCU, UCC, UCA, UCG, AGU, AGC Glu/E GAA, GAG
		 * Thr/T ACU, ACC, ACA, ACG Gly/G GGU, GGC, GGA, GGG Trp/W UGG His/H
		 * CAU, CAC Tyr/Y UAU, UAC Ile/I AUU, AUC, AUA Val/V GUU, GUC, GUA, GUG
		 * START AUG STOP UAG, UGA, UAA
		 */
		private String other;

		ToProtein(String other)
		{
			this.other = other;
		}

		public String getProtein()
		{
			return this.other;
		}
	}

	public enum Translate
	{
		A("T"), B("V"), C("G"), D("H"), G("C"), H("D"), K("M"), M("K"), N("N"), R("Y"), S("S"), T("A"), V("B"), W("W"), Y(
				"R");

		private String other;

		Translate(String other)
		{
			this.other = other;
		}

		public String getCompliment()
		{
			return this.other;
		}
	}

	public static final int BUFFER_BYTE_SIZE = 8388608; // 8MB
	public static final NumberFormat DECIMAL_FORMAT = new DecimalFormat("############.########");
	public static final int FASTA_LINE_LENGTH = 60;

	public static final int MBYTES = 1048576;

	public static int checkForEnd(String line, int brackets)
	{
		if (line.startsWith("{"))
		{
			brackets++;
		}
		if (line.startsWith("}"))
		{
			brackets--;
		}
		if (brackets == 0)
		{
			return -1;
		}

		return brackets;
	}

	public static long[] computeHashes(String item, int numWords, int seed)
	{
		long[] hashes = new long[numWords];

		for (int word = 0; word < numWords; word += 2)
		{
			HashFunction hashFunc = Hashing.murmur3_128(seed + word);
			Hasher hasher = hashFunc.newHasher();
			hasher.putUnencodedChars(item);

			// get the two longs out
			HashCode hc = hasher.hash();
			ByteBuffer bb = ByteBuffer.wrap(hc.asBytes());
			hashes[word] = bb.getLong(0);
			if (word + 1 < numWords)
				hashes[word + 1] = bb.getLong(8);
		}

		return hashes;
	}

	public final static int[] computeHashesInt(Object obj, int numWords, int seed)
	{
		if (obj instanceof Integer)
			return computeHashesIntInt((Integer) obj, numWords, seed);
		if (obj instanceof Long)
			return computeHashesIntLong((Long) obj, numWords, seed);
		if (obj instanceof Double)
			return computeHashesIntDouble((Double) obj, numWords, seed);
		if (obj instanceof Float)
			return computeHashesIntFloat((Float) obj, numWords, seed);
		if (obj instanceof String)
			return computeHashesIntString((String) obj, numWords, seed);

		throw new MhapRuntimeException("Cannot hash class type " + obj.getClass().getCanonicalName());
	}

	public final static int[] computeHashesIntDouble(double obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];

		HashFunction hf = Hashing.murmur3_32(seed);

		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putDouble(obj).putInt(iter).hash();

			hashes[iter] = hc.asInt();
		}

		return hashes;
	}

	public final static int[] computeHashesIntFloat(float obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];

		HashFunction hf = Hashing.murmur3_32(seed);

		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putFloat(obj).putInt(iter).hash();

			hashes[iter] = hc.asInt();
		}

		return hashes;
	}

	public final static int[] computeHashesIntInt(int obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];

		HashFunction hf = Hashing.murmur3_32(seed);

		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putInt(obj).putInt(iter).hash();

			hashes[iter] = hc.asInt();
		}

		return hashes;
	}

	public final static int[] computeHashesIntLong(long obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];

		HashFunction hf = Hashing.murmur3_32(seed);

		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putLong(obj).putInt(iter).hash();

			hashes[iter] = hc.asInt();
		}

		return hashes;
	}

	public final static int[] computeHashesIntString(String obj, int numWords, int seed)
	{
		int[] hashes = new int[numWords];

		HashFunction hf = Hashing.murmur3_32(seed);

		for (int iter = 0; iter < numWords; iter++)
		{
			HashCode hc = hf.newHasher().putUnencodedChars(obj).putInt(iter).hash();

			hashes[iter] = hc.asInt();
		}

		return hashes;
	}

	public final static long[][] computeKmerHashes(final String seq, final int kmerSize, final int numWords, final int seed)
	{
		final int numberKmers = seq.length()-kmerSize+1;

		if (numberKmers < 1)
			throw new MhapRuntimeException("Kmer size bigger than string length.");

		// get the rabin hashes
		final long[] rabinHashes = computeSequenceHashesLong(seq, kmerSize, seed);

		final long[][] hashes = new long[rabinHashes.length][numWords];

		// Random rand = new Random(0);
		for (int iter = 0; iter < rabinHashes.length; iter++)
		{
			// rand.setSeed(rabinHashes[iter]);
			long x = rabinHashes[iter];

			for (int word = 0; word < numWords; word++)
			{
				// hashes[iter][word] = rand.nextLong();

				// XORShift Random Number Generators
				x ^= (x << 21);
				x ^= (x >>> 35);
				x ^= (x << 4);
				hashes[iter][word] = x;
			}
		}
		
		return hashes;
	}
	
	public final static long[][] computeKmerHashesExact(final String seq, final int kmerSize, final int numWords, final int seed)
	{
		HashFunction hf = Hashing.murmur3_128(seed);

		long[][] hashes = new long[seq.length() - kmerSize + 1][numWords];
		for (int iter = 0; iter < hashes.length; iter++)
		{
			String subStr = seq.substring(iter, iter + kmerSize);
			
			for (int word=0; word<numWords; word++)
			{
				HashCode hc = hf.newHasher().putUnencodedChars(subStr).putInt(word).hash();
				hashes[iter][word] = hc.asLong();
			}
		}
		
		return hashes;
	}
	
	public final static int[] computeSequenceHashes(final String seq, final int kmerSize)
	{
		// RollingSequenceHash rabinHash = new RollingSequenceHash(kmerSize);
		// final int[] rabinHashes = rabinHash.hashInt(seq);

		HashFunction hf = Hashing.murmur3_32(0);

		int[] hashes = new int[seq.length() - kmerSize + 1];
		for (int iter = 0; iter < hashes.length; iter++)
		{
			HashCode hc = hf.newHasher().putUnencodedChars(seq.substring(iter, iter + kmerSize)).hash();
			hashes[iter] = hc.asInt();
		}

		return hashes;
	}
	
	public final static long[] computeSequenceHashesLong(final String seq, final int kmerSize, final int seed)
	{
		HashFunction hf = Hashing.murmur3_128(seed);

		long[] hashes = new long[seq.length() - kmerSize + 1];
		for (int iter = 0; iter < hashes.length; iter++)
		{
			HashCode hc = hf.newHasher().putUnencodedChars(seq.substring(iter, iter + kmerSize)).hash();
			hashes[iter] = hc.asLong();
		}

		return hashes;
	}

	// add new line breaks every FASTA_LINE_LENGTH characters
	public final static String convertToFasta(String supplied)
	{
		StringBuffer converted = new StringBuffer();
		int i = 0;
		String[] split = supplied.trim().split("\\s+");
		if (split.length > 1)
		{ // process as a qual
			int size = 0;
			for (i = 0; i < split.length; i++)
			{
				converted.append(split[i]);
				size += split[i].length();
				if (i != (split.length - 1))
				{
					if (size >= FASTA_LINE_LENGTH)
					{
						size = 0;
						converted.append("\n");
					}
					else
					{
						converted.append(" ");
					}
				}
			}
		}
		else
		{
			for (i = 0; (i + FASTA_LINE_LENGTH) < supplied.length(); i += FASTA_LINE_LENGTH)
			{
				converted.append(supplied.substring(i, i + FASTA_LINE_LENGTH));
				converted.append("\n");
			}
			converted.append(supplied.substring(i, supplied.length()));
		}
		return converted.toString();
	}
	
	public final static int countLetterInRead(String fasta, String letter)
	{
		return countLetterInRead(fasta, letter, false);
	}

	public final static int countLetterInRead(String fasta, String letter, Boolean caseSensitive)
	{
		String ungapped = Utils.getUngappedRead(fasta);
		int len = ungapped.length();
		if (len == 0)
		{
			return -1;
		}

		int increment = letter.length();
		int count = 0;

		for (int i = 0; i <= ungapped.length() - increment; i += increment)
		{
			if (letter.equals(ungapped.substring(i, i + increment)) && caseSensitive)
			{
				count++;
			}
			if (letter.equalsIgnoreCase(ungapped.substring(i, i + increment)) && !caseSensitive)
			{
				count++;
			}
		}
		return count;
	}

	public final static HashSet<Long> createKmerFilter(String fileName, double maxPercent, int kmerSize, int seed)
			throws IOException
	{
		File file = new File(fileName);

		// make sure don't leak resources
		try (BufferedReader bf = new BufferedReader(new FileReader(file), BUFFER_BYTE_SIZE);)
		{
			// generate hashset
			ArrayList<Long> filterArray = new ArrayList<Long>();

			String line = bf.readLine();
			while (line != null)
			{
				String[] str = line.split("\\s+", 3);

				if (str.length < 2)
					throw new MhapRuntimeException(
							"Kmer filter file must have at least two column [kmer kmer_percent].");

				double percent = Double.parseDouble(str[1]);

				// if greater, add to hashset
				if (percent > maxPercent)
				{
					long[] minHash = Utils.computeSequenceHashesLong(str[0], kmerSize, seed);

					if (minHash.length > 1)
						System.err.println("Warning filter file kmer size larger than setting!");

					for (long val : minHash)
						filterArray.add(val);
				}
				else
					break;

				// read the next line
				line = bf.readLine();
			}
			return new HashSet<Long>(filterArray);
		}
	}

	public final static int[] errorString(int[] s, double readError)
	{
		int[] snew = s.clone();

		Random generator = new Random();
		for (int iter = 0; iter < s.length; iter++)
		{
			if (generator.nextDouble() < readError)
				while (snew[iter] == s[iter])
					snew[iter] = generator.nextInt(3);
		}

		return snew;
	}

	public final static int[] generateInts(long seed, int numWords)
	{
		int[] hashes = new int[numWords];

		long x = seed;
		for (int word = 0; word < numWords; word++)
		{
			// hashes[iter][word] = rand.nextInt();

			// XORShift Random Number Generators
			x ^= (x << 21);
			x ^= (x >>> 35);
			x ^= (x << 4);
			hashes[word] = (int) x;
		}

		return hashes;
	}

	public final static long[] generateLongs(long seed, int numWords)
	{
		long[] hashes = new long[numWords];

		long x = seed;
		for (int word = 0; word < numWords; word++)
		{
			// hashes[iter][word] = rand.nextInt();

			// XORShift Random Number Generators
			x ^= (x << 21);
			x ^= (x >>> 35);
			x ^= (x << 4);
			hashes[word] = x;
		}

		return hashes;
	}

	public final static BufferedReader getFile(String fileName, String postfix) throws Exception
	{
		String[] array = new String[1];
		array[0] = postfix;

		return getFile(fileName, array);
	}

	public final static BufferedReader getFile(String fileName, String[] postfix) throws IOException
	{
		BufferedReader bf = null;

		if (fileName.endsWith("bz2"))
		{
			// open file as a pipe
			System.err.println("Running command " + "bzip2 -dc " + new File(fileName).getAbsolutePath() + " |");
			Process p = Runtime.getRuntime().exec("bzip2 -dc " + new File(fileName).getAbsolutePath() + " |");
			bf = new BufferedReader(new InputStreamReader(p.getInputStream()), BUFFER_BYTE_SIZE);
			System.err.println(bf.ready());
		}
		else if (fileName.endsWith("gz"))
		{
			// open file as a pipe
			System.err.println("Runnning comand " + "gzip -dc " + new File(fileName).getAbsolutePath() + " |");
			Process p = Runtime.getRuntime().exec("gzip -dc " + new File(fileName).getAbsolutePath() + " |");
			bf = new BufferedReader(new InputStreamReader(p.getInputStream()), BUFFER_BYTE_SIZE);
			System.err.println(bf.ready());
		}
		else
		{
			int i = 0;
			for (i = 0; i < postfix.length; i++)
			{
				if (fileName.endsWith(postfix[i]))
				{
					bf = new BufferedReader(new FileReader(fileName), BUFFER_BYTE_SIZE);
					break;
				}
			}
			if (i == postfix.length)
			{
				System.err.println("Unknown file format " + fileName + " Skipping!");
			}
		}

		return bf;
	}

	public final static String getID(String line)
	{
		String ids[] = line.split(":");
		int commaPos = ids[1].indexOf(",");
		if (commaPos != -1)
		{
			return ids[1].substring(1, commaPos).trim();
		}
		else
		{
			return ids[1];
		}
	}

	public final static double getLetterPercentInRead(String fasta, String letter)
	{
		int ungappedLen = getUngappedRead(fasta).length();
		int count = countLetterInRead(fasta, letter);

		return count / (double) ungappedLen;
	}

	public final static int getOvlSize(int readA, int readB, int ahang, int bhang)
	{
		if ((ahang <= 0 && bhang >= 0) || (ahang >= 0 && bhang <= 0))
		{
			return -1;
		}

		if (ahang < 0)
		{
			return readA - Math.abs(bhang);
		}
		else
		{
			return readA - ahang;
		}
	}

	public final static int getRangeOverlap(int startA, int endA, int startB, int endB)
	{
		int minA = Math.min(startA, endA);
		int minB = Math.min(startB, endB);
		int maxA = Math.max(startA, endA);
		int maxB = Math.max(startB, endB);

		int start = Math.max(minA, minB);
		int end = Math.min(maxA, maxB);

		return (end - start + 1);
	}

	public final static String getUngappedRead(String fasta)
	{
		fasta = fasta.replaceAll("N", "");
		fasta = fasta.replaceAll("-", "");

		assert (fasta.length() >= 0);

		return fasta;
	}

	public final static String getValue(String line, String key)
	{
		if (line.startsWith(key))
		{
			return line.split(":")[1];
		}

		return null;
	}

	public final static <H> double hashEfficiency(HashMap<Integer, ArrayList<H>> c)
	{
		double e = hashEnthropy(c);
		double log2inv = 1.0 / Math.log(2);
		double scaling = Math.log(c.size()) * log2inv;

		return e / scaling;
	}

	public final static <H> double hashEnthropy(HashMap<Integer, ArrayList<H>> c)
	{
		double sum = 0.0;
		double log2inv = 1.0 / Math.log(2);

		double[] p = new double[c.size()];
		int size = 0;
		int count = 0;
		for (ArrayList<H> elem : c.values())
		{
			size += elem.size();
			p[count++] = elem.size();
		}

		for (int iter = 0; iter < p.length; iter++)
		{
			double val = p[iter] / (double) size;
			sum -= val * Math.log(val) * log2inv;
		}

		return sum;
	}

	public final static boolean isAContainedInB(int startA, int endA, int startB, int endB)
	{
		int minA = Math.min(startA, endA);
		int minB = Math.min(startB, endB);
		int maxA = Math.max(startA, endA);
		int maxB = Math.max(startB, endB);

		return (minB < minA && maxB > maxA);
	}

	public final static Pair<Double, Double> linearRegression(int[] a, int[] b, int size)
	{
		// take one pass and compute means
		int xy = 0;
		int x = 0;
		int y = 0;
		int x2 = 0;

		for (int iter = 0; iter < size; iter++)
		{
			xy += a[iter] * b[iter];
			x += a[iter];
			y += b[iter];
			x2 += a[iter] * a[iter];
		}

		double Ninv = 1.0 / (double) size;

		double beta = ((double) xy - Ninv * (double) (x * y)) / ((double) x2 - Ninv * (double) (x * x));
		double alpha = Ninv * ((double) y - beta * (double) x);

		return new Pair<Double, Double>(alpha, beta);
	}

	public final static double mean(double[] a, int size)
	{
		double x = 0.0;
		for (int iter = 0; iter < size; iter++)
			x += a[iter];

		return x / (double) size;
	}

	public final static double mean(int[] a, int size)
	{
		int x = 0;
		for (int iter = 0; iter < size; iter++)
			x += a[iter];

		return x / (double) size;
	}

	public final static double pearsonCorr(int[] a, int[] b, int size)
	{
		if (size < 2)
			return 0.0;

		double meana = mean(a, size);
		double meanb = mean(b, size);
		double stda = std(a, size, meana);
		double stdb = std(b, size, meanb);

		double r = 0.0;
		for (int iter = 0; iter < size; iter++)
		{
			r += ((double) a[iter] - meana) * ((double) b[iter] - meanb) / (stda * stdb);
		}

		return r / (double) (size - 1);
	}

	// adapted form
	// http://blog.teamleadnet.com/2012/07/quick-select-algorithm-find-kth-element.html
	public final static int quickSelect(int[] array, int k, int length)
	{
		if (array == null || length <= k)
			return Integer.MAX_VALUE;

		int from = 0;
		int to = length - 1;

		// if from == to we reached the kth element
		while (from < to)
		{
			int r = from;
			int w = to;
			int mid = array[(r + w) / 2];

			// stop if the reader and writer meets
			while (r < w)
			{
				if (array[r] >= mid)
				{
					// put the large values at the end
					int tmp = array[w];
					array[w] = array[r];
					array[r] = tmp;
					w--;
				}
				else
				{
					// the value is smaller than the pivot, skip
					r++;
				}
			}

			// if we stepped up (r++) we need to step one down
			if (array[r] > mid)
				r--;

			// the r pointer is on the end of the first k elements
			if (k <= r)
			{
				to = r;
			}
			else
			{
				from = r + 1;
			}
		}

		return array[k];
	}

	public final static String rc(String supplied)
	{
		StringBuilder st = new StringBuilder();
		for (int i = supplied.length() - 1; i >= 0; i--)
		{
			char theChar = supplied.charAt(i);

			if (theChar != '-')
			{
				Translate t = Translate.valueOf(Character.toString(theChar).toUpperCase());
				st.append(t.getCompliment());
			}
			else
			{
				st.append("-");
			}
		}
		return st.toString();
	}

	public final static double std(double[] a, int size, double mean)
	{
		double x = 0.0;
		for (int iter = 0; iter < size; iter++)
		{
			double val = a[iter] - mean;
			x += val * val;
		}

		return Math.sqrt(x / (double) (size - 1));
	}

	public final static double std(int[] a, int size, double mean)
	{
		double x = 0.0;
		for (int iter = 0; iter < size; iter++)
		{
			double val = (double) a[iter] - mean;
			x += val * val;
		}

		return Math.sqrt(x / (double) (size - 1));
	}

	public final static String toProtein(String genome, boolean isReversed, int frame)
	{
		StringBuilder result = new StringBuilder();

		if (isReversed)
		{
			genome = rc(genome);
		}
		genome = genome.replaceAll("-", "");

		for (int i = frame; i < (genome.length() - 3); i += 3)
		{
			String codon = genome.substring(i, i + 3);
			String protein = ToProtein.valueOf(codon).getProtein();
			result.append(protein);
		}

		return result.toString();
	}

	public static String toString(double[][] A)
	{
		StringBuilder s = new StringBuilder();

		s.append("[");

		for (double[] a : A)
		{
			if (a != null)
			{
				for (int iter = 0; iter < a.length - 1; iter++)
					s.append("" + a[iter] + ",");

				if (a.length > 0)
					s.append("" + a[a.length - 1]);
			}
			s.append("\n");
		}
		s.append("]");

		return new String(s);
	}

	public static String toString(long[][] A)
	{
		StringBuilder s = new StringBuilder();

		s.append("[");

		for (long[] a : A)
		{
			if (a != null)
			{
				for (int iter = 0; iter < a.length - 1; iter++)
					s.append("" + a[iter] + ",");

				if (a.length > 0)
					s.append("" + a[a.length - 1]);
			}
			s.append("\n");
		}
		s.append("]");

		return new String(s);
	}
}
