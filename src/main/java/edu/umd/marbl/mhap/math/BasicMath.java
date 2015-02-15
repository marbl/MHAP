/* 
 * ARMOR package
 * 
 * This  software is distributed "as is", without any warranty, including 
 * any implied warranty of merchantability or fitness for a particular
 * use. The authors assume no responsibility for, and shall not be liable
 * for, any special, indirect, or consequential damages, or any damages
 * whatsoever, arising out of or in connection with the use of this
 * software.
 * 
 * Copyright (c) 2012 by Konstantin Berlin 
 * University Of Maryland
 * 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */
package edu.umd.marbl.mhap.math;

/**
 * The Class BasicMath.
 */
public final class BasicMath
{

	/** The Constant PI. */
	public static final double PI = Math.PI;

	/** The Constant TWOPI. */
	public static final double TWOPI = 2.0 * PI;

	public static double abs(double a)
	{
		return Math.abs(a);
	}

	public static double[] abs(double[] a)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = BasicMath.abs(a[iter]);

		return val;
	}
	
	/**
	 * Acos.
	 * 
	 * @param x
	 *            the x
	 * @return the double
	 */
	public static double acos(double x)
	{
		return FastMath.acos(x);
	}

	/**
	 * Adds the.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] add(final double[] a, final double b)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] + b;

		return val;
	}

	/**
	 * Adds the.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] add(final double[] a, final double[] b)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] + b[iter];

		return val;
	}

	/**
	 * Angle.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double
	 */
	public final static double angle(final double[] a, final double[] b)
	{
		double angle = acos(normalizedDotProduct(a, b));

		return angle;
	}

	public static double angleAbsolute(double[] a, double[] b)
	{
		return Math.min(Math.abs(angle(a, b)), Math.abs(angle(a, BasicMath.mult(b, -1.0))));
	}

	/**
	 * Asin.
	 * 
	 * @param x
	 *            the x
	 * @return the double
	 */
	public final static double asin(final double x)
	{
		return FastMath.asin(x);
	}

	public final static double[][] catColumns(final double[][] A, final double[][] B)
	{
		if (A.length != B.length)
			throw new MathRuntimeException("Number of rows must be equal in A and B.");

		double[][] C = new double[A.length][A[0].length + B[0].length];

		for (int row = 0; row < C.length; row++)
		{
			for (int column = 0; column < A[row].length; column++)
				C[row][column] = A[row][column];

			for (int column = 0; column < B[row].length; column++)
				C[row][A[row].length + column] = B[row][column];
		}

		return C;
	}

	/**
	 * Closest power of two.
	 * 
	 * @param a
	 *            the a
	 * @return the int
	 */
	public final static int closestPowerOfTwo(final int a)
	{
		int power = a == 0 ? 0 : 32 - Integer.numberOfLeadingZeros(a - 1);

		return 1 << power;
	}

	/**
	 * Cos.
	 * 
	 * @param angle
	 *            the angle
	 * @return the double
	 */
	public final static double cos(final double angle)
	{
		return FastMath.cos(angle);
	}

	/**
	 * Creates the identity matrix.
	 * 
	 * @param m
	 *            the m
	 * @param n
	 *            the n
	 * @return the double[][]
	 */
	public final static double[][] createIdentityMatrix(final int m, final int n)
	{
		double[][] A = new double[m][n];

		for (int iterRow = 0; iterRow < A.length; iterRow++)
		{
			for (int iterColumn = 0; iterColumn < A[iterRow].length; iterColumn++)
			{
				A[iterRow][iterColumn] = 0.0;
				if (iterRow == iterColumn)
					A[iterRow][iterColumn] = 1.0;
			}
		}

		return A;
	}

	/**
	 * Cube.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public final static double cube(final double a)
	{
		return a * a * a;
	}

	public final static double det(final double[][] A)
	{
		if (A == null || A.length != 3 || A[0].length != 3)
			throw new MathRuntimeException("Currently can only compute determinant of 3x3 matrix.");

		double det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) - A[1][0]
				* (A[0][1] * A[2][2] - A[2][1] * A[0][2]) + A[2][0] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);

		return det;
	}

	/**
	 * Divide.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] divide(final double[] a, final double b)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] / b;

		return val;
	}

	/**
	 * Divide.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] divide(final double[] a, final double[] b)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] / b[iter];

		return val;
	}

	/**
	 * Dot product.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double
	 */
	public final static double dotProduct(final double[] a, final double[] b)
	{
		if (a.length != b.length)
			throw new MathRuntimeException("Vector lengths must be equal.");

		double val = 0.0;

		for (int iter = 0; iter < a.length; iter++)
			val += a[iter] * b[iter];

		return val;
	}

	/**
	 * Euclidean distance.
	 * 
	 * @param x1
	 *            the x1
	 * @param y1
	 *            the y1
	 * @param z1
	 *            the z1
	 * @param x2
	 *            the x2
	 * @param y2
	 *            the y2
	 * @param z2
	 *            the z2
	 * @return the double
	 */
	public final static double euclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return sqrt(euclideanDistanceSquared(x1, y1, z1, x2, y2, z2));
	}

	/**
	 * Euclidean distance squared.
	 * 
	 * @param x1
	 *            the x1
	 * @param x2
	 *            the x2
	 * @return the double
	 */
	public final static double euclideanDistanceSquared(final double x1, final double x2)
	{
		double xdif = x2 - x1;

		return xdif * xdif;
	}

	/**
	 * Euclidean distance squared.
	 * 
	 * @param x1
	 *            the x1
	 * @param y1
	 *            the y1
	 * @param x2
	 *            the x2
	 * @param y2
	 *            the y2
	 * @return the double
	 */
	public final static double euclideanDistanceSquared(double x1, double y1, double x2, double y2)
	{
		double xdif = x2 - x1;
		double ydif = y2 - y1;

		return (xdif * xdif + ydif * ydif);
	}

	/**
	 * Euclidean distance squared.
	 * 
	 * @param x1
	 *            the x1
	 * @param y1
	 *            the y1
	 * @param z1
	 *            the z1
	 * @param x2
	 *            the x2
	 * @param y2
	 *            the y2
	 * @param z2
	 *            the z2
	 * @return the double
	 */
	public final static double euclideanDistanceSquared(double x1, double y1, double z1, double x2, double y2, double z2)
	{
		double xdif = x2 - x1;
		double ydif = y2 - y1;
		double zdif = z2 - z1;

		return (xdif * xdif + ydif * ydif + zdif * zdif);
	}

	public static boolean hasNaN(double[] x)
	{
		for (double val : x)
			if (Double.isNaN(val))
				return true;

		return false;
	}

	/**
	 * Checks if is identity matrix.
	 * 
	 * @param A
	 *            the a
	 * @return true, if is identity matrix
	 */
	public static boolean isIdentityMatrix(double[][] A)
	{
		if (A == null)
			return false;

		if (A.length != A[0].length)
			return false;

		for (int iterRow = 0; iterRow < A.length; iterRow++)
			for (int iterColumn = 0; iterColumn < A[iterRow].length; iterColumn++)
			{
				if (iterRow == iterColumn)
				{
					if (A[iterRow][iterColumn] != 1.0)
						return false;
				}
				else if (A[iterRow][iterColumn] != 0.0)
					return false;

			}

		return true;
	}

	public static boolean isNonNegative(double[] x)
	{
		for (double val : x)
			if (val < 0)
				return false;

		return true;
	}

	/*
	 * public static long nearestPow2(long x) { double logX =
	 * Math.log10(x)/Math.log10(2); if (Math.round(logX)<=logX) return x; else
	 * return (int)Math.pow(2, Math.floor(logX+1)); }
	 */

	public final static double laplanceProbabilty(double x, double b)
	{
		return 1.0/(2.0*b)*Math.exp(-Math.abs(x)/b);
	}

	/**
	 * Matrix to array.
	 * 
	 * @param A
	 *            the a
	 * @return the double[]
	 */
	public static double[] matrixToArray(double A[][])
	{
		double[] val = new double[A.length * A[0].length];

		for (int iterRow = 0; iterRow < A.length; iterRow++)
		{
			for (int iterColumn = 0; iterColumn < A[iterRow].length; iterColumn++)
				val[iterRow * A[0].length] = A[iterRow][iterColumn];
		}

		return val;
	}

	/**
	 * Max.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public final static double max(final double[] a)
	{
		double val = a[0];
		for (double elem : a)
			val = Math.max(val, elem);

		return val;
	}

	/**
	 * Min.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public final static double min(final double[] a)
	{
		double val = a[0];
		for (double elem : a)
			val = Math.min(val, elem);

		return val;
	}

	/**
	 * Mult.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] mult(final double[] a, final double b)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] * b;

		return val;
	}

	/**
	 * Mult.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] mult(final double[] a, final double[] b)
	{
		if (a == null || b == null)
			throw new MathRuntimeException("Arrays cannot be null.");

		if (a.length != b.length)
			throw new MathRuntimeException("Arrays must be of equal length.");

		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] * b[iter];

		return val;
	}

	/**
	 * Mult.
	 * 
	 * @param A
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[][]
	 */
	public final static double[][] mult(final double[][] A, final double b)
	{
		double[][] X = new double[A.length][A[0].length];

		for (int iterRow = 0; iterRow < A.length; iterRow++)
		{
			for (int iterColumn = 0; iterColumn < A[iterRow].length; iterColumn++)
			{
				X[iterRow][iterColumn] = A[iterRow][iterColumn] * b;
			}
		}

		return X;
	}

	public final static double[] mult(final double[][] A, final double[] b)
	{
		if (A == null || b == null)
			throw new java.lang.NullPointerException("Values cannot be null.");

		if (A[0].length != b.length)
			throw new MathRuntimeException("Matrix dimension [" + A.length + ", " + A[0].length
					+ "] does not match vector length " + b.length + ".");

		double[] x = new double[A.length];

		for (int iterRow = 0; iterRow < A.length; iterRow++)
		{
			x[iterRow] = 0.0;
			for (int iterColumn = 0; iterColumn < A[iterRow].length; iterColumn++)
			{
				x[iterRow] += A[iterRow][iterColumn] * b[iterColumn];
			}
		}

		return x;
	}

	public final static double[][] mult(final double[][] A, final double[][] B)
	{
		if (A == null || B == null)
			throw new java.lang.NullPointerException("Matrices cannot be null.");

		if (A[0].length != B.length)
			throw new MathRuntimeException("Matrices' dimensions do not match.");

		double[][] C = new double[A.length][B[0].length];

		for (int row = 0; row < A.length; row++)
		{
			for (int col = 0; col < B[0].length; col++)
			{
				C[row][col] = 0.0;
				for (int iter = 0; iter < B.length; iter++)
					C[row][col] += A[row][iter] * B[iter][col];
			}
		}

		return C;
	}

	public final static double[] multTranspose(final double[][] A, final double[] x)
	{
		double[] value = new double[A[0].length];

		for (int iterRow = 0; iterRow < A[0].length; iterRow++)
		{
			value[iterRow] = 0.0;
			for (int iterColumn = 0; iterColumn < A.length; iterColumn++)
			{
				value[iterRow] += A[iterColumn][iterRow] * x[iterColumn];
			}
		}

		return value;
	}

	public final static double[][] multTranspose(double[][] A, double[][] B)
	{
		if (A == null || B == null)
			throw new java.lang.NullPointerException("Matrices cannot be null.");

		if (A.length != B.length)
			throw new MathRuntimeException("Matrices' dimensions do not match.");

		double[][] C = new double[A[0].length][B[0].length];

		for (int colA = 0; colA < A[0].length; colA++)
		{
			for (int colB = 0; colB < B[0].length; colB++)
			{
				C[colA][colB] = 0.0;
				for (int iter = 0; iter < A.length; iter++)
					C[colA][colB] += A[iter][colA] * B[iter][colB];
			}
		}

		return C;
	}

	/**
	 * Nearest multiple.
	 * 
	 * @param n
	 *            the n
	 * @param base
	 *            the base
	 * @return the int
	 */
	public final static int nearestMultiple(int n, int base)
	{
		int x = n / base;

		if (x * base == n)
			return n;
		else
			return x * base + base;
	}

	public final static int[] nonZeroIndicies(final double[] x, final double absTolerance)
	{
		// count the number of elements
		int size = 0;
		for (int iter = 0; iter < x.length; iter++)
			if (Math.abs(x[iter]) > absTolerance)
				size++;

		// record the elements
		int[] list = new int[size];
		size = 0;
		for (int iter = 0; iter < x.length; iter++)
			if (Math.abs(x[iter]) > absTolerance)
			{
				list[size] = iter;
				size++;
			}

		return list;
	}

	public final static double[] nonZeroValues(final double[] x, final double absTolerance)
	{
		int[] list = nonZeroIndicies(x, absTolerance);
		double[] xnew = new double[list.length];

		int count = 0;
		for (int index : list)
		{
			xnew[count] = x[index];
			count++;
		}

		return xnew;
	}

	/*
	 * public static double[] legendrePolynomial(int n, double x) throws
	 * ArithmeticException { double P[] = new double[n + 1];
	 * 
	 * P[0] = 1.0; P[1] = x;
	 * 
	 * for (int m = 1; m < n - 1; m++) { P[m + 1] = ((2.0 * m + 1.0) * x * P[m]
	 * - m * P[m - 1]) / (m + 1.0); }
	 * 
	 * return P; }
	 */

	/*
	 * static public double[] normalizedlegendrePolynomial(int n, double x)
	 * throws ArithmeticException { double[] P = legendrePolynomial(n, x);
	 * 
	 * double norm = BasicMath.sqrt(n + .5); double sign = -1;
	 * 
	 * for (int m = 0; m < P.length; m++) { P[m] = sign * P[m] / norm; sign =
	 * -sign; }
	 * 
	 * return P; }
	 */

	/**
	 * Norm.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public static double norm(double[] a)
	{
		return BasicMath.sqrt(normSquared(a));
	}

	public final static double normalizedDotProduct(final double[] a, final double[] b)
	{
		return dotProduct(a, b) / (norm(a) * norm(b));
	}

	/**
	 * Norm squared.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public final static double normSquared(final double[] a)
	{
		double r = 0.0;

		for (double elem : a)
			r += elem * elem;

		return r;
	}

	/**
	 * Round to nearest.
	 * 
	 * @param x
	 *            the x
	 * @param n
	 *            the n
	 * @return the double
	 */
	public final static double roundToNearest(final double x, final int n)
	{
		double shift = Math.pow(10, n);
		return Math.round(x * shift) / shift;
	}

	/**
	 * Sin.
	 * 
	 * @param angle
	 *            the angle
	 * @return the double
	 */
	public final static double sin(final double angle)
	{
		return FastMath.sin(angle);
	}

	/**
	 * Sinc.
	 * 
	 * @param x
	 *            the x
	 * @return the double
	 */
	public final static double sinc(final double x)
	{
		return (x == 0 || (x < 1.0e-8 && x > -1.0e-8)) ? 1.0 : sin(x) / x;
	}

	/**
	 * Sqrt.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public final static double sqrt(final double a)
	{
		return FastMath.sqrt(a);
	}

	/**
	 * Square.
	 * 
	 * @param a
	 *            the a
	 * @return the double
	 */
	public final static double square(final double a)
	{
		return a * a;
	}

	/**
	 * Square.
	 * 
	 * @param a
	 *            the a
	 * @return the double[]
	 */
	public final static double[] square(final double[] a)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] * a[iter];

		return val;
	}

	/**
	 * Subtract.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] subtract(final double[] a, final double b)
	{
		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] - b;

		return val;
	}

	/**
	 * Subtract.
	 * 
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the double[]
	 */
	public final static double[] subtract(final double[] a, final double[] b)
	{
		if (a.length != b.length)
			throw new MathRuntimeException("Vectors must be of same length.");

		double[] val = new double[a.length];

		for (int iter = 0; iter < a.length; iter++)
			val[iter] = a[iter] - b[iter];

		return val;
	}

	public final static double sum(final double[] a)
	{
		if (a == null)
			return 0.0;

		double sum = 0.0;
		for (double val : a)
			sum += val;

		return sum;
	}

	public final static double[][] transpose(double[][] A)
	{
		if (A == null)
			return null;

		double[][] At = new double[A[0].length][A.length];

		for (int row = 0; row < A.length; row++)
			for (int col = 0; col < A[row].length; col++)
				At[col][row] = A[row][col];

		return At;
	}

}