package edu.umd.marbl.mhap.utils;

import java.util.Random;

public class CharacterHash
{
	public int hashvalues[] = new int[1 << 16];

	static CharacterHash charhash = new CharacterHash();

	public static CharacterHash getInstance()
	{
		return charhash;
	}

	public CharacterHash()
	{
		Random r = new Random(1);
		for (int k = 0; k < this.hashvalues.length; ++k)
			this.hashvalues[k] = r.nextInt();
	}

}
