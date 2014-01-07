package com.secret.fastalign.utils;

public class FastAlignRuntimeException extends RuntimeException
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 56387323839744808L;

	public FastAlignRuntimeException()
	{
		super();
	}

	public FastAlignRuntimeException(String message, Throwable cause, boolean enableSuppression,
			boolean writableStackTrace)
	{
		super(message, cause, enableSuppression, writableStackTrace);
	}

	public FastAlignRuntimeException(String message, Throwable cause)
	{
		super(message, cause);
	}

	public FastAlignRuntimeException(String message)
	{
		super(message);
	}

	public FastAlignRuntimeException(Throwable cause)
	{
		super(cause);
	}

	
}
