package edu.umd.marbl.mhap.main;

import edu.umd.marbl.mhap.align.AlignElementDoubleSketch;
import edu.umd.marbl.mhap.align.AlignElementString;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.align.Alignment;
import edu.umd.marbl.mhap.impl.MinHashBitSequenceSubSketches;
import edu.umd.marbl.mhap.impl.OverlapInfo;
import edu.umd.marbl.mhap.sketch.MinHashBitSketch;
import edu.umd.marbl.mhap.sketch.OrderedNGramHashes;
import edu.umd.marbl.mhap.utils.RandomSequenceGenerator;

public class AlignmentTry
{

	public static void main(String[] args)
	{
		String a = "bcdefghij1234567890";
		String b = "abcdefghij1234567890";
		
		RandomSequenceGenerator generator = new RandomSequenceGenerator();
		a = generator.generateRandomSequence(2000);
		b = a.substring(800, 1800);
		a = generator.addPacBioError(a);
		b = generator.addPacBioError(b);
		//b = generator.generateRandomSequence(1400);
		//b = a;
		
		Aligner<AlignElementString> aligner = new Aligner<AlignElementString>(true, -2.0, -1*Float.MAX_VALUE, 0.0);
		
		Alignment<AlignElementString> alignment = aligner.localAlignSmithWaterGotoh(new AlignElementString(a), new AlignElementString(b));
		
		System.err.println(alignment.getOverlapScore(5));
		
		System.out.println(alignment.outputAlignment());
		
		System.err.println("A1="+alignment.getA1());
		System.err.println("B1="+alignment.getB1());
		System.err.println("A2="+alignment.getA2());
		System.err.println("B2="+alignment.getB2());
		
		MinHashBitSequenceSubSketches m1 = new MinHashBitSequenceSubSketches(a, 7, 200, 20);
		MinHashBitSequenceSubSketches m2 = new MinHashBitSequenceSubSketches(b, 7, 200, 20);
		
		OverlapInfo info = m1.getOverlapInfo(new Aligner<AlignElementDoubleSketch<MinHashBitSketch>>(true, 0.00, 0.0, -0.52), m2);
				
		System.err.println("Compressed=");
		System.err.println(info.rawScore);
		System.err.println(info.a1);
		System.err.println(info.b1);
		System.err.println(info.a2);
		System.err.println(info.b2);		
		
		OverlapInfo info2 = m2.getOverlapInfo(new Aligner<AlignElementDoubleSketch<MinHashBitSketch>>(true, 0.00, 0.0, -0.52), m1);
		System.err.println("Swap=");
		System.err.println(info2.rawScore);
		System.err.println(info2.a1);
		System.err.println(info2.b1);
		System.err.println(info2.a2);
		System.err.println(info2.b2);		

		
		System.exit(1);
		
		OrderedNGramHashes hashes1 = new OrderedNGramHashes(a, 10);
		OrderedNGramHashes hashes2 = new OrderedNGramHashes(b, 10);
		
		System.err.println("Ordered=");
		System.err.println(hashes1.getOverlapInfo(hashes2, .2).a1);
		System.err.println(hashes1.getOverlapInfo(hashes2, .2).b1);
		System.err.println(hashes1.getOverlapInfo(hashes2, .2).a2);
		System.err.println(hashes1.getOverlapInfo(hashes2, .2).b2);
		
		/*
		SimHash s1 = new SimHash(a, kmerSize, 100);
		SimHash s2 = new SimHash(b, kmerSize, 100);
		
		MinHashSketch h1 = new MinHashSketch(a, kmerSize, 8000);
		MinHashSketch h2 = new MinHashSketch(b, kmerSize, 8000);
		
		MinHashBitSketch hb1 = new MinHashBitSketch(a, kmerSize, 100);
		MinHashBitSketch hb2 = new MinHashBitSketch(b, kmerSize, 100);

		System.err.println(s1.jaccard(s2));
		System.err.println(h1.jaccard(h2));
		System.err.println(hb1.jaccard(hb2));
		*/
	}

}
