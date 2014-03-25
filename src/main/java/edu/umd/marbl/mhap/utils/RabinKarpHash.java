package edu.umd.marbl.mhap.utils;

public final class RabinKarpHash
{

	private int BtoN;

	public int hashvalue;

	int n;

	final static int B = 31;
	//final static int B = 1664525;

	static CharacterHash hasher = CharacterHash.getInstance();

	// this is purely for testing purposes
	public static int nonRollingHash(CharSequence s)
	{
		int value = 0;
		for (int i = 0; i < s.length(); ++i)
		{
			char c = s.charAt(i);
			int z = hasher.hashvalues[c];
			value = B * value + z;
		}
		return value;
	}

	// myn is the length in characters of the blocks you want to hash
	public RabinKarpHash(int myn)
	{
		this.hashvalue = 0;
		this.n = myn;
		this.BtoN = 1;
		for (int i = 0; i < this.n; ++i)
		{
			this.BtoN *= B;
		}
	}

	// add new character (useful to initiate the hasher)
	// return 32 bits (not even universal)
	public int eat(char c)
	{
		this.hashvalue = B * this.hashvalue + hasher.hashvalues[c];
		return this.hashvalue;
	}

	// remove old character and add new one
	// return 32 bits (not even universal)
	public int update(char outchar, char inchar)
	{
		this.hashvalue = B * this.hashvalue + hasher.hashvalues[inchar] - this.BtoN * hasher.hashvalues[outchar];
		return this.hashvalue;
	}

}