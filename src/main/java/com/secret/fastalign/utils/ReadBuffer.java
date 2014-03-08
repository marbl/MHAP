package com.secret.fastalign.utils;

public final class ReadBuffer
{
	private byte[] buff = new byte[2];
	
	public final byte[] getBuffer(int size)
	{
		if (this.buff.length<size)
			this.buff = new byte[size];

		return this.buff;
	}
	
	
}