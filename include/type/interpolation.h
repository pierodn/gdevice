#pragma once

// TODO use generic typename for domain (add typename S)

//
// 1D
//

template <typename S, typename T> inline T first( S x, T y1, T y2)
{
	return y1;
}

template <typename S, typename T> inline T nearest( S x, T y1, T y2)
{
	return x<0.5 ? y1 : y2;
}

template <typename S, typename T> inline T lerp( S x, const T& y1, const T& y2) 
{
	return (1-x)*y1 + x*y2;
}

template <typename S, typename T> inline T coserp( S x, const T& y1, const T& y2) 
{
	return lerp( (1-cos(PI*x))/2, y1, y2 );
}

template <typename S, typename T> inline T hermite3( S x, const T& y1, const T& y2) 
{
	return lerp( x*x*(3-2*x), y1, y2 );
}
/*
// PROP
float scurve( float x, float average, float steepness ) 
{
	return 1/(1+exp(-steepness*(x-average)));
}
float gaussian( float x, float average, float steepness ) 
{
	return exp(-steepness*(x-average)*(x-average));
}

float fast_gaussian( float x, float average, float width ) //iq
{
	x = fabsf(x - average);
    if( x>width ) return 0.0f;
    x /= width;
    return 1.0f - x*x*(3.0f-2.0f*x);
}


*/

// TODO quintic ??

/*
template <typename T> inline T hermite5( float x, T y1, T y2) 
{
	return lerp( x*x*x*(10-15*x+6*x*x), y1, y2 );
}

template <typename T> inline T berp( float x, T y1, T y2) 
{
	return y1 + (y2-y1) * (sin(x*PI*(0.2+2.5*x*x*x))*pow(1-x,2.2) + x) * (1+(1.2*(1-x)));
}

template <typename T> inline T bounce( float x, T y1, T y2) 
{
	//return lerp( abs( sin(2*PI*(x+1)*(x+1))*(1-x)), y1, y2 );
	return y1;
}

template <typename T> inline T bounce( float x, T y1, T y2)
{
	return y1;
}*/
/*
template <typename T> inline T smoothstep( float x, T min, T max) 
{
	x = clamp( x, min, max );
    float y1 = (x-min)/(max-min);
    float y2 = (x-min)/(max-min);
    return -2*y1 * y1 *y1 + 3*y2 * y2;
}
*/
template <typename T> inline T clerp( float x, T y1 , T y2 )
{
    float min = 0.0f;
    float max = 360.0f;
    float half = abs( (max-min)/2.0f );

	return	y2-y1 < -half ? y1 + (max-y1+y2)*x :
			y2-y1 > +half ? y1 - (max-y2+y1)*x :
							y1 + (y2-y1)*x;
}

// interpolates x in [0,1] over a cubic curve passing for 2 points 
// (-1,y0), (0,y1), (1,y2), (2,y3)
template <typename T> inline T cubic( float x, T y0, T y1, T y2, T y3 )
{
	T p = (y3 - y2) - (y0 - y1);
	T q = (y0 - y1) - p;
	T r = y2 - y0;
	T s = y1;
	return p*x*x*x + q*x*x + r*x + s;
}



//
// 2D
//

template <typename T> inline T first( float x, float y, T z11, T z12, T z21, T z22 )
{
	return z11;
}

template <typename T> inline T nearest( float x, float y, T z11, T z12, T z21, T z22 )
{
	return x<0.5 ? ( y<0.5 ? z11 : z12 ) : ( y<0.5 ? z21 : z22 );
}

template <typename T> inline T lerp( float x, float y, const T& z11, const T& z12, const T& z21, const T& z22 )
{
	return z11*(1-x)*(1-y) + z12*(1-x)*y + z21*x*(1-y) + z22*x*y;
}

template <typename T> inline T coserp( float x, float y, const T& z11, const T& z12, const T& z21, const T& z22 )
{
	return lerp( (1-cos(x*PI))/2, (1-cos(y*PI))/2, z11, z12, z21, z22 );
}



// TODO bilinear
// http://en.wikipedia.org/wiki/Bilinear_interpolation
// f(x,y) = f(0,0)*(1-x)*(1-y) + f(1,0)*x*(1-y) + f(0,1)*(1-x)*y + f(1,1)*x*y

// TODO biquadratic


// TODO bicubic 
// Perlin cubic interpolator using Catmull-Rom:  http://mrl.nyu.edu/~perlin/java/Bicubic.html 
// Perlin cubic and bicubic with more matrices: http://mrl.nyu.edu/~perlin/cubic/Cubic_java.html           




// http://www.paulinternet.nl/?page=bicubic

/**
 * Interpolates the value of a point in a two dimensional surface using bicubic interpolation.
 * The value is calculated using the position of the point and the values of the 16 surrounding points.
 * The coordinate system works this way: The point p[x][y] has position (x-1, y-1).
 * For example, if you call this function with x=0 and y=1, p[1][2] will be returned.
 * You can pass any value for x and y, but you usually want them to be between 0 and 1.
 * Note that the returned value can be more or less than any of the values of the surrounding points. 
 * 
 * @param p A 4x4 array containing the values of the 16 surrounding points
 * @param x The horizontal position, between 0 and 1
 * @param y The vertical position, between 0 and 1
 * @return the interpolated value
 */
/*
double bicubicInterpolate (double* p, double x, double y) {
	double a00 = p[1][1];
	double a01 = p[1][2] - p[1][1]/2 - p[1][0]/3 - p[1][3]/6;
	double a02 = p[1][0]/2 - p[1][1] + p[1][2]/2;
	double a03 = p[1][1]/2 - p[1][0]/6 - p[1][2]/2 + p[1][3]/6;
	double a10 = p[2][1] - p[1][1]/2 - p[0][1]/3 - p[3][1]/6;
	double a11 = p[0][0]/9 + p[0][1]/6 - p[0][2]/3 + p[0][3]/18 + p[1][0]/6 + p[1][1]/4 - p[1][2]/2 + p[1][3]/12 - p[2][0]/3 - p[2][1]/2 + p[2][2] - p[2][3]/6 + p[3][0]/18 + p[3][1]/12 - p[3][2]/6 + p[3][3]/36;
	double a12 = p[0][1]/3 - p[0][0]/6 - p[0][2]/6 - p[1][0]/4 + p[1][1]/2 - p[1][2]/4 + p[2][0]/2 - p[2][1] + p[2][2]/2 - p[3][0]/12 + p[3][1]/6 - p[3][2]/12;
	double a13 = p[0][0]/18 - p[0][1]/6 + p[0][2]/6 - p[0][3]/18 + p[1][0]/12 - p[1][1]/4 + p[1][2]/4 - p[1][3]/12 - p[2][0]/6 + p[2][1]/2 - p[2][2]/2 + p[2][3]/6 + p[3][0]/36 - p[3][1]/12 + p[3][2]/12 - p[3][3]/36;
	double a20 = p[0][1]/2 - p[1][1] + p[2][1]/2;
	double a21 = p[0][2]/2 - p[0][1]/4 - p[0][0]/6 - p[0][3]/12 + p[1][0]/3 + p[1][1]/2 - p[1][2] + p[1][3]/6 - p[2][0]/6 - p[2][1]/4 + p[2][2]/2 - p[2][3]/12;
	double a22 = p[0][0]/4 - p[0][1]/2 + p[0][2]/4 - p[1][0]/2 + p[1][1] - p[1][2]/2 + p[2][0]/4 - p[2][1]/2 + p[2][2]/4;
	double a23 = p[0][1]/4 - p[0][0]/12 - p[0][2]/4 + p[0][3]/12 + p[1][0]/6 - p[1][1]/2 + p[1][2]/2 - p[1][3]/6 - p[2][0]/12 + p[2][1]/4 - p[2][2]/4 + p[2][3]/12;
	double a30 = p[1][1]/2 - p[0][1]/6 - p[2][1]/2 + p[3][1]/6;
	double a31 = p[0][0]/18 + p[0][1]/12 - p[0][2]/6 + p[0][3]/36 - p[1][0]/6 - p[1][1]/4 + p[1][2]/2 - p[1][3]/12 + p[2][0]/6 + p[2][1]/4 - p[2][2]/2 + p[2][3]/12 - p[3][0]/18 - p[3][1]/12 + p[3][2]/6 - p[3][3]/36;
	double a32 = p[0][1]/6 - p[0][0]/12 - p[0][2]/12 + p[1][0]/4 - p[1][1]/2 + p[1][2]/4 - p[2][0]/4 + p[2][1]/2 - p[2][2]/4 + p[3][0]/12 - p[3][1]/6 + p[3][2]/12;
	double a33 = p[0][0]/36 - p[0][1]/12 + p[0][2]/12 - p[0][3]/36 - p[1][0]/12 + p[1][1]/4 - p[1][2]/4 + p[1][3]/12 + p[2][0]/12 - p[2][1]/4 + p[2][2]/4 - p[2][3]/12 - p[3][0]/36 + p[3][1]/12 - p[3][2]/12 + p[3][3]/36;

	double x2 = x * x;
	double x3 = x2 * x;
	double y2 = y * y;
	double y3 = y2 * y;

	return a00 + a01 * y + a02 * y2 + a03 * y3 +
	       a10 * x + a11 * x * y + a12 * x * y2 + a13 * x * y3 +
	       a20 * x2 + a21 * x2 * y + a22 * x2 * y2 + a23 * x2 * y3 +
	       a30 * x3 + a31 * x3 * y + a32 * x3 * y2 + a33 * x3 * y3;
}
*/



//
// 3D
//

template <typename T> inline T first( float x, float y, float z,
			T z111, T z112, T z121, T z122, T z211, T z212, T z221, T z222 )
{
	return z111;
}

template <typename T> inline T nearest( float x, float y, float z,
			T z111, T z112, T z121, T z122, T z211, T z212, T z221, T z222 )
{
	return x<0.5 ? (y<0.5 ? (z<0.5 ? z111 : z112 ) : (z<0.5 ? z121 : z122 )) :
		(y<0.5 ? (z<0.5 ? z211 : z212 ) : (z<0.5 ? z221 : z222 ));
}

// Trilinear 
// http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

template <typename T> inline T lerp( float x, float y, float z,
			T z111, T z112, T z121, T z122, T z211, T z212, T z221, T z222 )
{
	return 	z111*(1-x)*(1-y)*(1-z) +
			z112*(1-x)*(1-y)*z     +
			z121*(1-x)*y    *(1-z) +
			z122*(1-x)*y    *z     +
			z211*x    *(1-y)*(1-z) +
			z212*x    *(1-y)*z     +
			z221*x    *y    *(1-z) +
			z222*x    *y    *z;
}


// TODO triquadratic 
// TODO tricubic





template<typename T, float* abscissa, T* ordinate, float interpolator(float, T, T)> 
inline T interpolate(float x)
{
	float* x1 = abscissa;
	T*     y1 = ordinate;

	for(; *x1<=x && *x1<1.0; x1++, y1++ ); 

	float* x0 = x1-1;
	T*     y0 = y1-1;

	return interpolator( (x-*x0)/(*x1-*x0), *y0, *y1 );
}

template<float* points>
inline float cubicinterpolate(float x)
{
	// ---x0---x1-x-x2---x3---

	float* x2;
	for( x2=points; *x2<=x && *x2<1.0; x2+=2 );

	float* x1 = x2==points ? x2 : x2-2; 
	float* x0 = x1==points ? x1 : x1-2; //
	float* x3 = *x2==1.0   ? x2 : x2+2; //

	float* y0 = x0+1;
	float* y1 = x1+1;
	float* y2 = x2+1;
	float* y3 = x3+1;

	return //coserp( (x-*x1)/(*x2-*x1), *y1, *y2 );
		cubic( (x-*x1)/(*x2-*x1), *y0, *y1, *y2, *y3 );
}