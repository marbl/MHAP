package edu.umd.marbl.mhap.sketch;

public class ZeroNGramsFoundException extends Exception
{

	private final String seqString;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3655558540692106680L;

	public ZeroNGramsFoundException(String message, String seqString)
	{
		super(message);
		this.seqString = seqString;
	}

	public ZeroNGramsFoundException(String message, Throwable cause, boolean enableSuppression,
			boolean writableStackTrace, String seqString)
	{
		super(message, cause, enableSuppression, writableStackTrace);
		this.seqString = seqString;
	}

	public ZeroNGramsFoundException(String message, Throwable cause, String seqString)
	{
		super(message, cause);
		this.seqString = seqString;
	}

	public ZeroNGramsFoundException(Throwable cause, String seqString)
	{
		super(cause);
		this.seqString = seqString;
	}
	
	public String getSequenceString()
	{
		return this.seqString;
	}

}
