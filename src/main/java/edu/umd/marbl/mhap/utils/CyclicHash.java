/**
 * Daniel Lemire, Owen Kaser: Recursive n-gram hashing is pairwise independent, at best, Computer Speech & Language, Volume 24, Issue 4, October 2010, Pages 698-710 http://arxiv.org/abs/0705.4676
 */
package edu.umd.marbl.mhap.utils;

public final class CyclicHash
{
	public int hashvalue;

	int myr;

	int n;

	private final static CharacterHash hasher = CharacterHash.getInstance();

	public final static int wordsize = 32;

	private final static int fastleftshift1(int x)
	{
		return (x << 1) | (x >>> (wordsize - 1));
	}

	// this is purely for testing purposes
	public final static int nonRollingHash(CharSequence s)
	{
		int value = 0;
		for (int i = 0; i < s.length(); ++i)
		{
			char c = s.charAt(i);
			int z = hasher.hashvalues[c];
			value = fastleftshift1(value) ^ z;
		}
		return value;
	}

	// myn is the length in characters of the blocks you want to hash
	public CyclicHash(int myn)
	{
		this.n = myn;
		if (this.n > wordsize)
		{
			throw new IllegalArgumentException();
		}

	}

	// add new character (useful to initiate the hasher)
	// to get a strongly universal hash value, you have to ignore the last or
	// first (n-1) bits.
	public final int eat(char c)
	{
		this.hashvalue = fastleftshift1(this.hashvalue);
		this.hashvalue ^= hasher.hashvalues[c];
		return this.hashvalue;
	}

	private final int fastleftshiftn(int x)
	{
		return (x << this.n) | (x >>> (wordsize - this.n));
	}

	// remove old character and add new one
	// to get a strongly universal hash value, you have to ignore the last or
	// first (n-1) bits.
	public final int update(char outchar, char inchar)
	{
		int z = fastleftshiftn(hasher.hashvalues[outchar]);
		this.hashvalue = fastleftshift1(this.hashvalue) ^ z ^ hasher.hashvalues[inchar];
		return this.hashvalue;
	}

}