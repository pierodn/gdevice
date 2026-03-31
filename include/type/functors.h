#pragma once

//#include "utils/math/algebra/interpolators.h"


// TODO f(x,y)
// TODO f(x,y), f(image)


// 
// Weight smoothers
//
// for cubic-like interpolation, still using 2 points instead of 4.
//
// Examples:
//
// lerp( scurve3(a), y1, y2 )
// lerp( cosine(a), y1, y2 )
// 
// refer to: http://libnoise.sourceforge.net/noisegen/index.html
//
//



//template<typename T>
inline float zero( float )
{
	return float(0);
}

//template<typename T>
inline float unit( float )
{
	return float(1);
}

inline float unit( float, float )
{
	return float(1);
}

template<typename T, T k>
inline T constant( float )
{
	return k;
}
/*
template<typename T>
inline T identity( T x )
{
	return T(x);
}

template<typename T>
inline T cosine( T a )
{
	return T((1-cos(a*PI))/2);
}

template<typename T>
inline T scurve3( T a )
{
	return T( a*a*(3-2*a) );
}

template<typename T>
inline T scurve5( T a )
{
	float a3 = a*a*a;
	float a4 = a3*a;
	float a5 = a4*a;
	return T(6*a5 - 15*a4 + 10*a3);
}
*/
/*
template< typename T, T A(T), T B(T) >
inline T add( T x, T y )
{
	return A(x)+B(y);
}
*/

template< float A(float), float B(float) >
inline float add( float x, float y )
{
	return A(x)+B(y);
}
