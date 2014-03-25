package edu.umd.marbl.mhap.utils;

import java.lang.reflect.Array;

/**
 * The Class HashCodeUtil.
 */
public final class HashCodeUtil {

  /// PRIVATE ///
  /** The oD d_ prim e_ number. */
  private static final int fODD_PRIME_NUMBER = 37;

  /** The Constant SEED. */
  public static final int SEED = 23;

  /**
	 * First term.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @return the int
	 */
  private static int firstTerm( int aSeed ){
    return fODD_PRIME_NUMBER * aSeed;
  }

  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aBoolean
	 *          the a boolean
	 * @return the int
	 */
  public static int hash( int aSeed, boolean aBoolean ) {
    return firstTerm( aSeed ) + ( aBoolean ? 1 : 0 );
  }

  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aChar
	 *          the a char
	 * @return the int
	 */
  public static int hash( int aSeed, char aChar ) {
    //System.out.println("char...");
    return firstTerm( aSeed ) + aChar;
  }
  
  public final static int hash(int aSeed, char[] charArray, int start, int size)
  {
  	int hash = 0;
  	for (int iter=0; iter<size; iter++)
  	{
  		hash += firstTerm(aSeed) + charArray[start+iter];
  	}
  	
  	return hash;
  }

  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aDouble
	 *          the a double
	 * @return the int
	 */
  public static int hash( int aSeed , double aDouble ) {
    return hash( aSeed, Double.doubleToLongBits(aDouble) );
  }

  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aFloat
	 *          the a float
	 * @return the int
	 */
  public static int hash( int aSeed , float aFloat ) {
    return hash( aSeed, Float.floatToIntBits(aFloat) );
  }

  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aInt
	 *          the a int
	 * @return the int
	 */
  public static int hash( int aSeed , int aInt ) {
    /*
    * Implementation Note
    * Note that byte and short are handled by this method, through
    * implicit conversion.
    */
    //System.out.println("int...");
    return firstTerm( aSeed ) + aInt;
  }


  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aLong
	 *          the a long
	 * @return the int
	 */
  public static int hash( int aSeed , long aLong ) {
    //System.out.println("long...");
    return firstTerm(aSeed)  + (int)( aLong ^ (aLong >>> 32) );
  }

  /**
	 * Hash.
	 * 
	 * @param aSeed
	 *          the a seed
	 * @param aObject
	 *          the a object
	 * @return the int
	 */
  public static int hash( int aSeed , Object aObject ) {
    int result = aSeed;
    if ( aObject == null) {
      result = hash(result, 0);
    }
    else if ( ! isArray(aObject) ) {
      result = hash(result, aObject.hashCode());
    }
    else {
      int length = Array.getLength(aObject);
      for ( int idx = 0; idx < length; ++idx ) {
        Object item = Array.get(aObject, idx);
        //recursive call!
        result = hash(result, item);
      }
    }
    return result;
  }

  /**
	 * Checks if is array.
	 * 
	 * @param aObject
	 *          the a object
	 * @return true, if is array
	 */
  private static boolean isArray(Object aObject){
    return aObject.getClass().isArray();
  }
} 
