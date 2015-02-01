package edu.umd.marbl.mhap.main;

import edu.umd.marbl.mhap.align.AlignElementString;
import edu.umd.marbl.mhap.align.Aligner;
import edu.umd.marbl.mhap.align.Alignment;
import edu.umd.marbl.mhap.sketch.MinHash;
import edu.umd.marbl.mhap.sketch.MinHashBitSketch;
import edu.umd.marbl.mhap.sketch.SimHash;

public class AlignmentTry
{

	public static void main(String[] args)
	{
		String a = "1234567890abcdefghij";
		String b = "1234567890";
		
		Aligner<AlignElementString> aligner = new Aligner<AlignElementString>(true, -0.5, 0.0);
		
		Alignment<AlignElementString> alignment = aligner.localAlignSmithWater(new AlignElementString(b), new AlignElementString(a));
		
		System.out.println(alignment.outputAlignment());
		
		SimHash s1 = new SimHash(a, 1, 100);
		SimHash s2 = new SimHash(b, 1, 100);
		
		MinHash h1 = new MinHash(a, 1, 8000);
		MinHash h2 = new MinHash(b, 1, 8000);
		
		MinHashBitSketch hb1 = new MinHashBitSketch(a, 1, 100);
		MinHashBitSketch hb2 = new MinHashBitSketch(b, 1, 100);

		System.err.println(s1.jaccard(s2));
		System.err.println(h1.jaccard(h2));
		System.err.println(hb1.jaccard(hb2));
	}

}
