package com.secret.fastalign.simhash;

import com.secret.fastalign.data.SequenceId;

public interface VectorHash<T extends VectorHash<T>>
{
	double jaccord(T val);

	SequenceId getSequenceId();
}
