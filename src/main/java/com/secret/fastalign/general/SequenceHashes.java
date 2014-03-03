package com.secret.fastalign.general;

import java.io.Serializable;

public interface SequenceHashes extends Serializable
{
	public byte[] getAsByteArray();
	public SequenceId getSequenceId();
	public int getSequenceLength();
}
