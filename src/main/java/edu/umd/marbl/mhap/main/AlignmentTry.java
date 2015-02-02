package edu.umd.marbl.mhap.main;

import edu.umd.marbl.mhap.align.AlignElementSketch;
import edu.umd.marbl.mhap.align.AlignElementString;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.align.Alignment;
import edu.umd.marbl.mhap.impl.OverlapInfo;
import edu.umd.marbl.mhap.sketch.MinHashSketch;
import edu.umd.marbl.mhap.sketch.MinHashBitSketch;
import edu.umd.marbl.mhap.sketch.MinHashSketchSequence;
import edu.umd.marbl.mhap.sketch.OrderedNGramHashes;
import edu.umd.marbl.mhap.sketch.SimHash;
import edu.umd.marbl.mhap.utils.RandomSequenceGenerator;

public class AlignmentTry
{

	public static void main(String[] args)
	{
		String a = "bcdefghij1234567890";
		String b = "abcdefghij1234567890";
		
		int kmerSize = 16;
		
		RandomSequenceGenerator generator = new RandomSequenceGenerator(0);
		a = "Z"+generator.generateRandomSequence(2000)+"Z";
		b = a.substring(600, 1400);
		a = generator.addPacBioError(a);
		b = generator.addPacBioError(b);
		
		//String c= a;
		//a = b;
		//b = c;
		
		Aligner<AlignElementString> aligner = new Aligner<AlignElementString>(true, -2.0, -1.00);
		
		Alignment<AlignElementString> alignment = aligner.customAlignSmithWaterGotoh(new AlignElementString(a), new AlignElementString(b));
		
		System.out.println(alignment.outputAlignment());
		
		System.err.println("A1="+alignment.getA1());
		System.err.println("B1="+alignment.getB1());
		System.err.println("A2="+alignment.getA2());
		System.err.println("B2="+alignment.getB2());
		
		MinHashSketchSequence m1 = new MinHashSketchSequence(a, 8, 200, 20);
		MinHashSketchSequence m2 = new MinHashSketchSequence(b, 8, 200, 20);
		
		System.err.println("Size1="+m1.length());
		System.err.println("Size2="+m2.length());
		
		OverlapInfo info = m1.getOverlapInfo(new Aligner<AlignElementSketch<MinHashBitSketch>>(true, 0.0001, -10000.0), m2);
		System.err.println("Compressed=");
		System.err.println(info.a1);
		System.err.println(info.b1);
		System.err.println(info.a2);
		System.err.println(info.b2);		
		
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
