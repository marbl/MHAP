package edu.umd.marbl.mhap.sketch;

public class SketchRuntimeException extends RuntimeException
{

	/**
	 * 
	 */
	private static final long serialVersionUID = 8422390842382501317L;

	public SketchRuntimeException()
	{
	}

	public SketchRuntimeException(String message)
	{
		super(message);
	}

	public SketchRuntimeException(Throwable cause)
	{
		super(cause);
	}

	public SketchRuntimeException(String message, Throwable cause)
	{
		super(message, cause);
	}

	public SketchRuntimeException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace)
	{
		super(message, cause, enableSuppression, writableStackTrace); 
	}

}
