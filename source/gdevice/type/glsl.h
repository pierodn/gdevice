#pragma once

//#include <cmath>

//#define __SSE__

/*
	Only some operations (specialized for T=float and N=4) are 
	explicitly optimized with SSE intrinsics. 
	All others rely on specific compiler optimizations such as 
	unrolling, RVO and SSE code generation. 
	So wisely arrange your compiler configuration.
*/


//
// Abstract types
//
template<class T, int N> union vec;
template<class T, int N, int M=N> union mat;

//
// CPU types
//
typedef char				int8;
typedef short				int16;
typedef int					int32;
typedef long long			int64;
typedef unsigned char		uint8;
typedef unsigned short		uint16;
typedef unsigned int		uint32;
typedef unsigned long long	uint64;
typedef unsigned char		byte;
typedef vec<byte,4> 		rgba;
typedef vec<float,8>		vec8;		// TEMP terrain

#define PRECISION 16
typedef long long			fixed;
typedef vec<fixed,2>		fpvec2;
typedef vec<fixed,3>		fpvec3;



//
// GPU types
//
typedef unsigned int	 uint;
typedef vec<float, 2>    vec2;
typedef vec<float, 3>    vec3;
typedef vec<float, 4>    vec4;
typedef vec<double,2>   dvec2;
typedef vec<double,3>   dvec3;
typedef vec<double,4>   dvec4;
typedef vec<bool,  2>   bvec2;
typedef vec<bool,  3>   bvec3;
typedef vec<bool,  4>   bvec4;
typedef vec<int,   2>   ivec2;
typedef vec<int,   3>   ivec3;
typedef vec<int,   4>   ivec4;
typedef vec<uint,  2>   uvec2; 
typedef vec<uint,  3>   uvec3;
typedef vec<uint,  4>   uvec4;
typedef mat<float, 2>    mat2;
typedef mat<float, 3>    mat3;
typedef mat<float, 4>    mat4;
typedef mat<float, 2,2>  mat2x2;
typedef mat<float, 2,3>  mat2x3;
typedef mat<float, 2,4>  mat2x4;
typedef mat<float, 3,2>  mat3x2;
typedef mat<float, 3,3>  mat3x3;
typedef mat<float, 3,4>  mat3x4;
typedef mat<float, 4,2>  mat4x2;
typedef mat<float, 4,3>  mat4x3;
typedef mat<float, 4,4>  mat4x4;
typedef mat<double,2>   dmat2;
typedef mat<double,3>   dmat3;
typedef mat<double,4>   dmat4;
typedef mat<double,2,2> dmat2x2;
typedef mat<double,2,3> dmat2x3;
typedef mat<double,2,4> dmat2x4;
typedef mat<double,3,2> dmat3x2;
typedef mat<double,3,3> dmat3x3;
typedef mat<double,3,4> dmat3x4;
typedef mat<double,4,2> dmat4x2;
typedef mat<double,4,3> dmat4x3;
typedef mat<double,4,4> dmat4x4;


//
// Constants
//
//#if !defined(PI)
const double PI			= 3.14159265358979323846;
//#endif 
const double LN2		= 0.693147180559945309417;
const double EPSILON	= 0.00001f;


//
// Abstract functions
//
template<class T> inline T radians(T x)			{ return (PI/180)*x; }
template<class T> inline T degrees(T x)			{ return (180/PI)*x; }
/*template<class T> inline T sin (T x)			{ return static_cast<T>(::sin(x)); }
template<class T> inline T cos (T x)			{ return static_cast<T>(::cos(x)); }
template<class T> inline T tan (T x)			{ return static_cast<T>(::tan(x)); }
template<class T> inline T asin(T x)			{ return static_cast<T>(::asin(x)); }
template<class T> inline T acos(T x)			{ return static_cast<T>(::acos(x)); }
template<class T> inline T atan(T x)			{ return static_cast<T>(::atan(x)); }
template<class T> inline T atan(T x, T y)		{ return static_cast<T>(::atan(x,y)); }
// TODO sinh, cosh, tanh, asinh, acosh, atanh
template<class T> inline T pow (T x, T y)		{ return static_cast<T>(::pow(x,y)); }
template<class T> inline T exp (T x)			{ return static_cast<T>(::exp(x)); }*/
template<class T> inline T log (T x)			{ return static_cast<T>(::log(x)); }
template<class T> inline T exp2(T x)			{ return static_cast<T>(::exp(x*LN2)); }
template<class T> inline T log2(T x)			{ return static_cast<T>(::log(x)/LN2); }
//template<class T> inline T sqrt(T x)			{ return static_cast<T>(::sqrt(x)); }
template<class T> inline T inversesqrt(T x)		{ return static_cast<T>(1/::sqrt(x)); }
//template<class T> inline T abs(T x)				{ return x<0 ? -x : x; }
template<class T> inline T sign(T x)			{ return x>0 ? 1.0 : x<0 ? -1.0 : 0.0; }
//template<class T> inline T floor(T x)			{ return static_cast<T>(::floor(x)); }
//template<class T> inline T ceil(T x)			{ return static_cast<T>(::ceil(x)); }
template<class T> inline T fract(const T x)			{ return static_cast<T>(x-::floor(x)); }
template<class T> inline T mod(T x, T y)		{ return static_cast<T>(::fmod(x,y)); }
template<class T> inline T clamp(T x,T x0=T(0),T x1=T(1)) { return x<x0 ? x0 : x>x1 ? x1 : x; }
template<class T> inline T mix(T x0, T x1, float a)	{ return x0*(1-a) + x1*a; }
template<class T> inline T step(T x0, T x)		{ return x<=x0 ? 0.0f : 1.0f; }
template<class T> inline T smoothstep(T x0, T x1, T x) { T t = clamp((x-x0)/(x1-x0), static_cast<T>(0), static_cast<T>(1)); return t*t*(3 - 2*t); }
template<class T> inline T length(T x)          { return ::sqrt(x*x); }
template<class T> inline T distance(T p0, T p1) { return length(p0-p1); }
template<class T> inline T dot(T x, T y)        { return x*y; }
template<class T> inline T normalize(T)         { return 1; }
template<class T> inline T faceforward(T n, T i, T nref) { return nref*i < 0 ? n : -n; }
template<class T> inline T reflect(T i, T n)    { return i - 2*dot(n,i)*n; }
template<class T> inline T refract(T i, T n, float eta) { T k = 1 - eta*eta*(1-dot(n,i)*dot(n,i)); return k<0 ? 0 : eta*i - (eta*dot(n,i) + sqrt(k))*n; }
#undef min
#undef max
template<class T> inline T min(T x, T y)		{ return x<y ? x : y; }
template<class T> inline T max(T x, T y)		{ return x>y ? x : y; }
template<class T> inline void swap(T& a, T& b)  { T t; t=a; a=b; b=t; }

//inline bool equal(float a, float b, float epsilon=EPSILON) { return fabs(a-b)<epsilon; }
template<class T> inline void equal(T a, T b, T epsilon=T(EPSILON)) { return abs(a-b)<epsilon; }
template<class T> inline T amod(T x, T y)		{ return (x/y - floor(x/y))*y; }

// TODO? lerp, coserp, hermite2, cubic..


//
// Definitions for generic T and N
//
template<class T, int N> union vec
{
	T array[N];

	inline T& operator[](int i) { return array[i]; }
    inline const T& operator[](int i) const { return array[i]; }

	inline vec() {}
	inline vec( const T& s ) { for( int i=0; i<N; i++ ) array[i] = s; }
	inline vec( const vec<T,N+1>& v ) { for( int i=0; i<N; i++ ) array[i] = v[i]; }
	inline vec( const vec<T,N-1>& v, const T s ) { for( int i=0; i<N-1; i++ ) array[i] = v[i]; array[N-1] = s; }
};

/* PROP1
// TODO use reinterpret_cast ? 

template<class T, int N, T f(const T&)> inline vec<T,N> functor1( const vec<T,N>& a )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = f(a[i]);
	return temp;
}
template<class T, int N, T f(const T&,const T&) > inline vec<T,N> functor2( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = f(a[i],b[i]);
	return temp;
}
template<class T, int N, T f(T&,T&,T&) > inline vec<T,N> 
	functor2( const vec<T,N>& a, const vec<T,N>& b, const vec<T,N>& c )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = f(a[i],b[i],c[i]);
	return temp;
}
*/

/* PROP2
template<class T, int N> inline vec<T,N> atan22( const vec<T,N>& a, const vec<T,N>& b ) { return functor<T,N,atan2>(a,b); }
*/

/* PROP3
VEC operator+	( VECREF v, SCALAR s)				RETURN_VEC( v[i]+s );
VEC abs			( VECREF v )						RETURN_VEC( fabs(v[i]) );
VEC clamp		( VECREF v, SCALAR v0, SCALAR v1 )	RETURN_VEC( clamp(v[i], v0[i], v1[i]) );
*/

template<class T, int N> inline vec<T,N> operator-( const vec<T,N>& a )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = -a[i];
	return temp;
}

template<class T, int N> inline vec<T,N> operator+( const vec<T,N>& a, const vec<T,N>& b )
{ 
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]+b[i];
	return temp;
}
template<class T, int N> inline vec<T,N> operator-( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]-b[i];
	return temp;
}
template<class T, int N> inline vec<T,N> operator*( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]*b[i];
	return temp;
}
template<class T, int N> inline vec<T,N> operator/( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]/b[i];
	return temp;
}

template<class T, int N> inline vec<T,N> operator+( const vec<T,N>& a, const T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]+s;
	return temp;
}
template<class T, int N> inline vec<T,N> operator-( const vec<T,N>& a, const T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]-s;
	return temp;
}
template<class T, int N> inline vec<T,N> operator*( const vec<T,N>& a, const T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]*s;
	return temp;
}
template<class T, int N> inline vec<T,N> operator/( const vec<T,N>& a, const T& s)
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]/s;
	return temp;
}

template<class T, int N> inline vec<T,N> operator+( const T& s, const vec<T,N>& a )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = s+a[i];
	return temp;
}
template<class T, int N> inline vec<T,N> operator-( const T& s, const vec<T,N>& a )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = s-a[i];
	return temp;
}
template<class T, int N> inline vec<T,N> operator*( const T& s, const vec<T,N>& a )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = s*a[i];
	return temp;
}
template<class T, int N> inline vec<T,N> operator/( const T& s, const vec<T,N>& a )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = s/a[i];
	return temp;
}

template<class T, int N> inline vec<T,N>& operator+=( vec<T,N>& a, const vec<T,N>& b)
{
	a = a + b;
	return a;
}
template<class T, int N> inline vec<T,N>& operator-=( vec<T,N>& a, const vec<T,N>& b)
{
	a = a - b;
	return a;
}
template<class T, int N> inline vec<T,N>& operator*=( vec<T,N>& a, const vec<T,N>& b)
{
	a = a * b;
	return a;
}
template<class T, int N> inline vec<T,N>& operator/=( vec<T,N>& a, const vec<T,N>& b)
{
	a = a / b;
	return a;
}

template<class T, int N> inline vec<T,N>& operator+=( vec<T,N>& a, const T& s )
{
	a = a + s;
	return a;
}
template<class T, int N> inline vec<T,N>& operator-=( vec<T,N>& a, const T& s )
{
	a = a - s;
	return a;
}
template<class T, int N> inline vec<T,N>& operator*=( vec<T,N>& a, const T& s )
{
	a = a * s;
	return a;
}
template<class T, int N> inline vec<T,N>& operator/=( vec<T,N>& a, const T& s )
{
	a = a / s;
	return a;
}

// vec functions
template<class T, int N> inline vec<T,N> radians( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = radians(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> degrees( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = degrees(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> sin( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = sin(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> cos( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = cos(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> tan( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = tan(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> asin( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = asin(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> acos( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = acos(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> atan( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = atan(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> atan( const vec<T,N>& y, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = atan2(y[i],x[i]);
	return temp;
}

template<class T, int N> inline vec<T,N> pow( const vec<T,N>& a, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = pow(a[i],x[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> exp2( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = exp2(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> log2( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = log2(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> exp( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = ::exp(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> log( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = ::log(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> sqrt( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = sqrt(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> inversesqrt( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = inversesqrt(v[i]);
	return temp;
}

template<class T, int N> inline vec<T,N> abs( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = fabs(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> sign( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = sign(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> floor( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = floor(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> ceil( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = ceil(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> fract( const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = fract(v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> mod( const vec<T,N>& v, T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = mod(v[i],s);
	return temp;
}
template<class T, int N> inline vec<T,N> mod( T& s, const vec<T,N>& v )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = mod(s,v[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> mod( const vec<T,N>& y, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = mod(y[i],x[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> min( const vec<T,N>& v, const T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = min(v[i],s);
	return temp;
}
template<class T, int N> inline vec<T,N> min( const vec<T,N>& y, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = min(y[i],x[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> max( const vec<T,N>& v, const T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = max(v[i],s);
	return temp;
}
template<class T, int N> inline vec<T,N> max( const vec<T,N>& y, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = max(y[i],x[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> clamp( const vec<T,N>& x, const T& x0, const T& x1 )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = clamp(x[i],x0,x1);
	return temp;
}
template<class T, int N> inline vec<T,N> clamp( const vec<T,N>& x, const vec<T,N>& x0, const vec<T,N>& x1 )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = clamp(x[i],x0[i],x1[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> mix( const vec<T,N>& a, const vec<T,N>& b, const float& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = mix(a[i],b[i],s);
	return temp;
}
template<class T, int N> inline vec<T,N> mix( const vec<T,N>& a, const vec<T,N>& b, const vec<float,N>& c )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = mix(a[i],b[i],c[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> step( const vec<T,N>& edge, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = step(edge[i],x[i]);
	return temp;
}
template<class T, int N> inline vec<T,N> smoothstep( const vec<T,N>& edge0, const vec<T,N>& edge1, const vec<T,N>& x )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = smoothstep(edge0[i],edge1[i],x[i]);
	return temp;
}

template<class T, int N> inline T length( const vec<T,N>& v )
{
	return static_cast<T>(sqrt(dot(v,v)));
}
template<class T, int N> inline T distance( const vec<T,N>& a, const vec<T,N>& b )
{
	return length(a-b);
}
template<class T, int N> inline T dot( const vec<T,N>& a, const vec<T,N>& b )
{
	T t = 0;
	for( int i=0; i<N; i++ ) t += a[i]*b[i];
	return t;
}

template<class T, int N> inline vec<T,N> normalize( const vec<T,N>& v )
{
	T len = length(v);
	if (len  == 0) return v;
	return v/len;
}
template<class T, int N> inline vec<T,N> faceforward( const vec<T,N>& n, const vec<T,N>& i, const vec<T,N>& nref )
{
	return dot(nref,i) < 0 ? n : -n;
}
template<class T, int N> inline vec<T,N> reflect( const vec<T,N>& i, const vec<T,N>& n )
{
	return i - 2*dot(n,i)*n;
}
template<class T, int N> inline vec<T,N> refract( const vec<T,N>& i, const vec<T,N>& n, const float& eta )
{
	T d = dot(n,i);
	T k = 1 - eta*eta*( 1 - d*d );
	return k<0 ? T(0) : eta*i - ( eta*d + sqrt(k) )*n;
}

template<class T, int N> inline bool operator==( const vec<T,N>& a, const vec<T,N>& b )
{
	for( int i=0; i<N; i++ ) if( a[i]!=b[i] ) return false;
	return true;
}
template<class T, int N> inline bool operator!=( const vec<T,N>& a, const vec<T,N>& b )
{
	for( int i=0; i<N; i++ ) if( a[i]==b[i] ) return false;
	return true;
}
template<class T, int N> inline bool equal( const vec<T,N>& a, const vec<T,N>& b, const T& epsilon=T(EPSILON) )
{
	for( int i=0; i<N; i++ ) if( abs(a[i]-b[i])>epsilon ) return false;
	return true;
}

template<class T, int N> inline vec<T,N> amod( const vec<T,N>& v, T& s )
{
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = amod(v[i],s);
	return temp;
}

// TODO noise1, noise2..


template<class T, int N, int M> union mat
{
	T array[N*M];

	inline T& operator[](int i)	{ return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline mat() {}
	inline mat( const T& s )		{ diag( vec<T,N>(s) ); }
	inline mat( const vec<T,N>& v ) { diag(v); }

	inline mat( const mat<T,N-1>& m ) 
	{ 
		for( int j=0; j<N-1; j++ )
		for( int i=0; i<N-1; i++ ) array[j*N+i] = m[j*N+i];

		for( int j=0; j<N-1; j++ ) array[j*N + N-1] = 0;
		for( int i=0; i<N-1; i++ ) array[(N-1)*N + i] = 0; 
		array[N*N-1] = T(1);
	}
	/*
// TEMP
	inline void diag( const vec<T,N>& v )
	{
		for( int i=0; i<N*N; i++) array[i] = T(0);
		for( int i=0; i<N; i++) array[i*N+i] = v[i];
	}*/
};

template<class T, int N> inline mat<T,N> operator+( const mat<T,N>& a, const mat<T,N>& b ) 
{
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = a[i]+b[i];
	return temp;
}
template<class T, int N> inline mat<T,N> operator-( const mat<T,N>& a, const mat<T,N>& b ) 
{
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = a[i]-b[i];
	return temp;
}
template<class T, int N> inline mat<T,N> operator+( const mat<T,N>& m, const T& s ) 
{
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = s+m[i];
	return temp;
}
template<class T, int N> inline mat<T,N> operator+( const T& s, const mat<T,N>& m )
{
    return m + s;
}
template<class T, int N> inline mat<T,N> operator-( const T& s, const mat<T,N>& m ) 
{
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = s-m[i];
	return temp;
}
template<class T, int N> inline mat<T,N> operator-( const mat<T,N>& m, const T& s ) 
{
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = m[i]-s;
	return temp;
}
template<class T, int N> inline mat<T,N> operator*( const T& s, const mat<T,N>& m ) 
{
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = s*m[i];
	return temp;
}
template<class T, int N> inline mat<T,N> operator*( const mat<T,N>& m, const T& s ) 
{ 
	return s*m; 
}

template<class T, int N> inline mat<T,N> operator/( const mat<T,N>& m, const T& s ) 
{ 
	mat<T,N> temp;
	for (int i = 0; i < N*N; i++) temp[i] = m[i]/s;
	return temp;
}

template<class T, int N> inline mat<T,N> operator/=( const mat<T,N>& m, const T& s ) 
{ 
	return m/s;
}

template<class T, int N> inline vec<T,N> operator*( const mat<T,N>& m, const vec<T,N>& v )
{
	vec<T,N> temp;
	for(int j=0; j<N; j++)
	{
		temp[j] = 0;
		for(int i=0; i<N; i++)
			temp[j] += m[j+N*i] * v[i];
	}
	return temp;
}

template<class T, int N> inline vec<T,N> operator*( const vec<T,N>& v, const mat<T,N>& m )
{
	vec<T,N> temp;
	for(int j=0; j<N; j++) 
	{
		temp[j] = 0;
		for(int i=0; i<N; i++)
			temp[j] += m[i+N*j] * v[i];
	}
	return temp;
}
template<class T, int N> inline vec<T,N>& operator*=( vec<T,N>& a, const mat<T,N>& m )
{
	a = a * m;
	return a;
}

// matrix-matrix multiplication
template<class T, int N> inline mat<T,N> operator*( const mat<T,N>& a, const mat<T,N>& b )
{
	mat<T,N> temp;
	for(int i=0; i<N; i++)
	for(int j=0; j<N; j++)
	{
		temp[j+i*N] = 0;
		for(int k=0; k<N; k++)
			temp[j+i*N] += a[j+k*N] * b[k+i*N];
	}
	return temp;
}


template<class T, int N> inline mat<T,N>& operator*=( mat<T,N>& a, const mat<T,N>& b )
{
	a = a * b;
	return a;
}

template<class T, int N, int M> inline mat<T,N,M> matrixcompmult( const mat<T,N,M>& a, const mat<T,N,M>& b )
{
	mat<T,N,M> temp;
	for(int i=0; i<N*M; i++) temp[i] = a[i]*b[i];
	return temp;
}

template<class T, int N, int M> inline mat<T,M,N> outerProduct( const vec<T,M>& a, const vec<T,N>& b )
{
	mat<T,N,M> temp;
	for(int j=0; j<M; j++)
	for(int i=0; i<N; i++)
		temp[j*N+i] = a[j]*b[i];
	return temp;
}

template<class T, int N, int M> inline mat<T,M,N> transpose( const mat<T,N,M>& a )
{
	mat<T,M,N> temp;
	for(int j=0; j<N; j++ )
	for(int i=0; i<M; i++ )
		temp[j*M+i] = a[i*N+j];
	return temp;
}
/*
template<class T, int N, int M> inline mat<T,M,N> determinant( const mat<T,N,M>& a )
{
	T temp;
	// TODO
	return temp;
}
template<class T, int N, int M> inline mat<T,M,N> inverse( const mat<T,N,M>& a )
{
	mat<T,M,N> temp;
	// TODO
	return temp;
}
*/
template<class T, int N, int M> inline bool equal( const mat<T,N,M>& a, const mat<T,N,M>& b, const T& epsilon=T(EPSILON) )
{
	for( int i=0; i<N*M; i++ ) if( !equal(a[i],b[i],epsilon) ) return false;
	return true;
}


//
// Partial specializations for N=2
//
template<class T> union vec<T,2>
{
    T array[2];
    struct { T x,y; };
	struct { T u,v; };
    struct { T s,t; };
	struct { T i,j; };
    struct { T width, height; };
	struct { T start, count; };

	inline T& operator[](int i) { return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline vec() {}
	inline vec( const T& s1 ) { x = y = s1; }
	inline vec( const T& s1, const T& s2 ) { x = s1; y = s2; }
	inline vec( const vec<T,3>& v ) { x = v[0]; y = v[1]; }
};

template<class T> union mat<T,2>
{
	T array[4];
	struct { T a,b,c,d; };

	inline T& operator[](int i)	{ return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline mat() {}
	inline mat(const T& s) { set( s, 0, 0, s ); }
	inline mat(const T& m0, const T& m1, const T& m2, const T& m3) { set( m0, m1, m2,m3 ); }
	inline mat(const vec<T,2>& v) { set( v[0], 0, 0, v[1] ); }

	inline void set( const T& m0, const T& m1, const T& m2, const T& m3 )
	{
		array[0] = m0;  array[1] = m1;
		array[2] = m2;  array[3] = m3;
	}
};

//
// Partial specializations for N=3
//
template<class T> union vec<T,3>
{
    T array[3];
    struct { T x,y,z; };
    struct { T r,g,b; };
    struct { T s,t,p; };
    struct { vec<T,2> xy; T z; };
    struct { vec<T,2> st; T p; };
    struct { vec<T,2> rg; T b; };
    struct { T x; vec<T,2> yz; };
    struct { T r; vec<T,2> gb; };
    struct { T s; vec<T,2> tp; };
    struct { T pitch, yaw, roll; };
    struct { T theta, phi, radius; };
    struct { T width, height, depth; };
	struct { T longitude, latitude, altitude; };

	inline T& operator[](int i) { return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline vec() {}
	inline vec( const T& s1  ) { x = y = z = s1; }
	inline vec( const T& s1, const T& s2, const T& s3=T(0) ) { x = s1; y = s2; z = s3; }
	inline vec( const vec<T,4>& b ) { x = b.x; y = b.y; z = b.z; }
	inline vec( const vec<T,2>& b, const T s1 ) { x = b.x; y = b.y; z = s1; }
};

template<class T> inline vec<T,3> cross( const vec<T,3>& a, const vec<T,3>& b )
{
	vec<T,3> temp;
	temp[0] = a[1]*b[2] - a[2]*b[1];
	temp[1] = a[2]*b[0] - a[0]*b[2];
	temp[2] = a[0]*b[1] - a[1]*b[0];
	return temp;
}

template<class T> union mat<T,3>
{
	T array[9];
	struct { vec<T,3> xAxis,yAxis,zAxis; };
	struct { vec<T,3> position,rotation,scale; };

	inline T& operator[](int i) { return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline mat() {}
	inline mat(const T& s)  { set( s, 0, 0, 0, s, 0, 0, 0, s ); }
	inline mat(const T& m0, const T& m1, const T& m2, 
			   const T& m3, const T& m4, const T& m5, 
			   const T& m6, const T& m7, const T& m8) {  set( m0, m1, m2, m3, m4, m5, m6, m7, m8 ); }
	
	inline mat(const vec<T,3>& v) { set( v[0], 0, 0, 0, v[1], 0,    0,    0, v[2] ); }
	inline mat(const vec<T,2>& v) { set(    0, 0, 0, 0,    0, 0, v[0], v[1], 1 ); }
	inline mat(const mat<T,2>& m) { set( m[0], m[1], 0, m[2], m[3], 0, 0, 0, 1 ); }
	inline mat(const mat<T,4>& m) { set( m[0], m[1], m[2], m[4], m[5], m[6], m[8], m[9], m[10]); }
	inline mat(const mat<T,2>& m, const vec<T,2>& v) { set( m[0], m[1], 0, m[2], m[3], 0, v[0], v[1], 1 ); }

	inline void set( const T& m0, const T& m1, const T& m2, 
					 const T& m3, const T& m4, const T& m5, 
					 const T& m6, const T& m7, const T& m8 )
	{
		array[0] = m0;  array[1] = m1;  array[2] = m2;
		array[3] = m3;  array[4] = m4;  array[5] = m5;
		array[6] = m6;  array[7] = m7;  array[8] = m8;
	}
};

//
// Partial specializations for N=4
//
template<class T> union vec<T,4>
{
    T array[4];
    struct { T x,y,z,w; };
    struct { T r,g,b,a; }; 
    struct { T s,t,p,q; };
	struct { vec<T,3> xyz; T w; };
    struct { vec<T,3> rgb; T a; };
    struct { vec<T,3> stp; T q; };
    struct { T x; vec<T,3> yzw; };
    struct { T r; vec<T,3> gba; };
    struct { T s; vec<T,3> tpq; };
    struct { T x; vec<T,2> yz; T w; };
    struct { T r; vec<T,2> gb; T a; };
    struct { T s; vec<T,2> tp; T q; }; 
    struct { vec<T,2> xy; vec<T,2> zw; };
    struct { vec<T,2> rg; vec<T,2> ba; };
    struct { vec<T,2> st; vec<T,2> pq; };

	inline T& operator[](int i) { return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline vec() {}
	inline vec( const T& s1 ) { x = y = z = w = s1; }
	inline vec( const T& s1, const T& s2, const T& s3=T(0), const T& s4=T(0) ) { x = s1; y = s2; z = s3; w = s4; }
//	inline vec( const vec<T,5>& v ) { x = v[0]; y = v[1]; z = v[2]; w = v[3]; }
	inline vec( const vec<T,3>& v, const T s1 ) { x = v[0]; y = v[1]; z = v[2]; w = s1; }
	inline vec( const vec<T,2>& v, const T s1, const T s2 ) { x = v[0]; y = v[1]; z = s1; w = s2; }
	inline vec( const T s, const vec<T,3>& v  ) { x = s; y = v[0]; z = v[1]; w = v[2]; }
};

template<class T> inline vec<T,4> cross( const vec<T,4>& a, const vec<T,4>& b )
{
	vec<T,4> temp;
	temp[0] = a[1]*b[2] - a[2]*b[1];
	temp[1] = a[2]*b[0] - a[0]*b[2];
	temp[2] = a[0]*b[1] - a[1]*b[0];
	temp[3] = 1;
	return temp;
}

template<class T> union mat<T,4>
{
	T array[16]; 
	struct { vec<T,4> xAxis,yAxis,zAxis,translation; };

	inline T& operator[](int i) { return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline mat() {}

	inline mat( const T& s )
	{
		set( s,0,0,0, 0,s,0,0, 0,0,s,0, 0,0,0,s );
	}

	inline mat( const T& m0,  const T& m1,  const T& m2,  const T& m3, 
				const T& m4,  const T& m5,  const T& m6,  const T& m7, 
				const T& m8,  const T& m9,  const T& m10, const T& m11, 
				const T& m12, const T& m13, const T& m14, const T& m15 )
	{
		set( m0,m1,m2,m3, m4,m5,m6,m7, m8,m9,m10,m11, m12,m13,m14,m15 );
	}

	inline mat(const vec<T,4>& v) { set( v[0],0,0,0, 0,v[1],0,0, 0,0,v[2],0, 0,0,0,v[3] ); }

	// for homogeneous coordinate
	inline mat(const vec<T,3>& v) { set( v[0],0,0,0, 0,v[1],0,0, 0,0,v[2],0, 0,0,0,1 ); }
	inline mat(const mat<T,3>& m) { set( m[0],m[1],m[2],0, m[3],m[4],m[5],0, m[6],m[7],m[8],0, 0,0,0,1 ); }

	inline void set( const T& m0,  const T& m1,  const T& m2,  const T& m3, 
					 const T& m4,  const T& m5,  const T& m6,  const T& m7, 
					 const T& m8,  const T& m9,  const T& m10, const T& m11, 
					 const T& m12, const T& m13, const T& m14, const T& m15 )
	{
		array[0]  = m0;  array[1]  = m1;  array[2]  = m2;  array[3]  = m3;
		array[4]  = m4;  array[5]  = m5;  array[6]  = m6;  array[7]  = m7;
		array[8]  = m8;  array[9]  = m9;  array[10] = m10; array[11] = m11;
		array[12] = m12; array[13] = m13; array[14] = m14; array[15] = m15;
	}
};

// TEMP
template<class T> union vec<T,8>
{
    T array[8];
	struct { vec<T,4> color; vec<T,4> detail; };
   
	inline T& operator[](int i) { return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline vec() {}
	inline vec( const T& s1 ) { left = vec<T,4>(s1); right = vec<T,4>(s1); }
	inline vec( const vec<T,4>& c, const vec<T,4>& d ) { color = c; detail = d; }
};




//
// Partial specializations for T==bool
//
template<class T, int N> inline vec<bool,N> lessThan( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i] < b[i];
	return temp;
}
template<class T, int N> inline vec<bool,N> lessThanEqual( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i] <= b[i];
	return temp;
}
template<class T, int N> inline vec<bool,N> greaterThan( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i] > b[i];
	return temp;
}
template<class T, int N> inline vec<bool,N> greaterThanEqual( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i] >= b[i];
	return temp;
}
template<class T, int N> inline vec<bool,N> equal( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = equal(a[i],b[i]); //a[i] == b[i];
	return temp;
}
template<class T, int N> inline vec<bool,N> notEqual( const vec<T,N>& a, const vec<T,N>& b )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = !equal(a[i],b[i]); //a[i] != b[i];
	return temp;
}
template<int N> inline bool any( const vec<bool,N>& v )
{
	for( int i=0; i<N; i++ ) if (v[i]) return true;
	return false;
}
template<int N> inline bool all( const vec<bool,N>& v )
{
	for( int i=0; i<N; i++ ) if (!v[i]) return false;
	return true;
}
#ifndef __GNUC__
template<int N> inline vec<bool,N> not( const vec<bool,N>& v )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = !v[i];
	return temp;
}
#endif
// this operator!() can be used instead of the above not()
template<int N> inline vec<bool,N> operator!( const vec<bool,N>& v )
{
	vec<bool,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = !v[i];
	return temp;
}



//
// Partial specializations for T=float
//

/* TODO verify 
template<> inline float sqrt(float x)
{
	float temp;
	_asm
	{
		movupd xmm0, x;
		sqrtss xmm0, xmm0;
		movupd temp, xmm0;
	}
	return temp;
}
*/



//
// Partial specializations for T=float and N=2
//

#if 0
template<> inline vec<float,2> operator*( const vec<float,2>& a, const vec<float,2>& b )
{ 
	vec<float,2> temp; 
/*	__asm {
		mov   eax,a
		movss xmm0,[eax]a.x
		movss xmm1,[eax]a.y
		mov   eax,b
		movss xmm2,[eax]b.x
		movss xmm4,[eax]b.y
		mulss xmm0,xmm2
		mulss xmm1,xmm4
		movss [eax]temp.x,xmm0
		movss [eax]temp.y,xmm1
	} */
	temp.x = a.x * b.x;
	temp.y = a.y * b.y;
	return temp;
}
#endif


// TODO transformation factory: translation matrix, rotation matrix, scale matrix

inline mat2 rotate( float angle )
{
	mat2 temp;

	float A  = (float) cos( angle * PI / (-180.0) ); // meno !?
	float B  = (float) sin( angle * PI / (-180.0) );

	temp[0] =  A;
	temp[1] = -B;
	temp[2] =  B;
	temp[3] =  A;

	return temp;
}

inline dmat2 rotate( double angle )
{
	dmat2 temp;

	double A  = cos( angle * PI / (-180.0f) ); // meno !?
	double B  = sin( angle * PI / (-180.0f) );

	temp[0] =  A;
	temp[1] = -B;
	temp[2] =  B;
	temp[3] =  A;

	return temp;
}


inline float det( const mat2& m )
{
	return m.a*m.d - m.b*m.c;
}



//
// Partial specializations for T=float and N=3
//

/* TODO debug
inline vec3 normalize( const vec3& v )
{
	vec3 temp;
	__asm {
		mov   eax,v
		movss xmm0,[eax]v.x
		movss xmm1,[eax]v.y
		movss xmm2,[eax]v.z
		movss xmm4,xmm0
		movss xmm5,xmm1
		movss xmm6,xmm2

		mulss xmm0,xmm0
		mulss xmm1,xmm1
		mulss xmm2,xmm2

		addss xmm0,xmm1
		addss xmm2,xmm0

		rsqrtss xmm2,xmm2

		mulss xmm4,xmm2
		mulss xmm5,xmm2
		mulss xmm6,xmm2

		movss [eax]temp.x,xmm4
		movss [eax]temp.y,xmm5
		movss [eax]temp.z,xmm6
	} 
	return temp;
}
*/



inline mat3 rotate( float pitch, float yaw, float roll )
{
	mat3 temp;
	float Cx, Sx, Cy, Sy, Cz, Sz, CxSy, SxSy;

	Cx  = (float) cos( pitch * PI / (-180.0) );
	Sx  = (float) sin( pitch * PI / (-180.0) );
	Cy  = (float) cos(   yaw * PI / (-180.0) );
	Sy  = (float) sin(   yaw * PI / (-180.0) );
	Cz  = (float) cos(  roll * PI / (-180.0) );
	Sz  = (float) sin(  roll * PI / (-180.0) );
	CxSy = Cx * Sy;
	SxSy = Sx * Sy;

	temp[0] =    Cy * Cz;
	temp[1] =  SxSy * Cz + Cx * Sz;
	temp[2] = -CxSy * Cz + Sx * Sz;
	temp[3] =   -Cy * Sz;
	temp[4] = -SxSy * Sz + Cx * Cz;
	temp[5] =  CxSy * Sz + Sx * Cz;
	temp[6] =    Sy;
	temp[7] =   -Sx * Cy;
	temp[8] =    Cx * Cy;

	return temp;
}

inline mat3 rotate( const vec3& rotation )
{
	return rotate( rotation.pitch, rotation.yaw, rotation.roll );
}

inline mat3 transform( float rotation, float scale, float translation )
{
	return mat3( mat2(scale)*rotate(rotation), vec2(translation) );
}


// TODO templatize
inline dmat3 rotate( double pitch, double yaw, double roll )
{
	dmat3 temp;
	double Cx, Sx, Cy, Sy, Cz, Sz, CxSy, SxSy;

	Cx  = (double) cos( pitch * PI / (-180.0) );
	Sx  = (double) sin( pitch * PI / (-180.0) );
	Cy  = (double) cos(   yaw * PI / (-180.0) );
	Sy  = (double) sin(   yaw * PI / (-180.0) );
	Cz  = (double) cos(  roll * PI / (-180.0) );
	Sz  = (double) sin(  roll * PI / (-180.0) );
	CxSy = Cx * Sy;
	SxSy = Sx * Sy;

	temp[0] =    Cy * Cz;
	temp[1] =  SxSy * Cz + Cx * Sz;
	temp[2] = -CxSy * Cz + Sx * Sz;
	temp[3] =   -Cy * Sz;
	temp[4] = -SxSy * Sz + Cx * Cz;
	temp[5] =  CxSy * Sz + Sx * Cz;
	temp[6] =    Sy;
	temp[7] =   -Sx * Cy;
	temp[8] =    Cx * Cy;

	return temp;
}

inline dmat3 rotate( const dvec3& rotation )
{
	return rotate( rotation.pitch, rotation.yaw, rotation.roll );
}

inline dmat3 transform( double rotation, double scale, double translation )
{
	return dmat3( dmat2(scale)*rotate(rotation), dvec2(translation) );
}




inline float determinant( const mat3& m )
{
	return 
		  m[0] * (m[4] * m[8] - m[7] * m[5])
		- m[3] * (m[1] * m[8] - m[7] * m[2])
		+ m[6] * (m[1] * m[5] - m[4] * m[2]);
}


inline mat3 inverse( const mat3& m )
{
	mat3 temp;
	temp[0] = + (m[4] * m[8] - m[5] * m[7]);
	temp[1] = - (m[1] * m[8] - m[2] * m[7]);
	temp[2] = + (m[1] * m[5] - m[2] * m[4]);
	temp[3] = - (m[3] * m[8] - m[5] * m[6]);
	temp[4] = + (m[0] * m[8] - m[2] * m[6]);
	temp[5] = - (m[0] * m[5] - m[2] * m[3]);
	temp[6] = + (m[3] * m[7] - m[4] * m[6]);
	temp[7] = - (m[0] * m[7] - m[1] * m[6]);
	temp[8] = + (m[0] * m[4] - m[1] * m[3]);
	temp /= determinant(m);

	return temp;
}

inline mat3 inverseTranspose( const mat3& m )
{
	mat3 temp;
	temp[0] = + (m[4] * m[8] - m[5] * m[7]);
	temp[1] = - (m[3] * m[8] - m[5] * m[6]);
	temp[2] = + (m[3] * m[7] - m[4] * m[6]);
	temp[3] = - (m[1] * m[8] - m[2] * m[7]);
	temp[4] = + (m[0] * m[8] - m[2] * m[6]);
	temp[5] = - (m[0] * m[7] - m[1] * m[6]);
	temp[6] = + (m[1] * m[5] - m[2] * m[4]);
	temp[7] = - (m[0] * m[5] - m[2] * m[3]);
	temp[8] = + (m[0] * m[4] - m[1] * m[3]);
	temp /= determinant(m);

	return temp;
}

// TODO inline mat4 inverse( const mat4& m )



//
// Partial specializations for T==float and N=4
// 

#if defined(__SSE__)
	#include <xmmintrin.h>

#if defined(_MSC_VER)
	#define ALIGN(x)            __declspec(align(x))
	#define ALIGNED_MALLOC(s,a) _aligned_malloc(s,a)
	#define ALIGNED_FREE(s)     _aligned_free(s)
#else
	#define ALIGN(x)            __attribute__ ((aligned(x)))
	#define ALIGNED_MALLOC(s,a) __mingw_aligned_malloc(s,a)
	#define ALIGNED_FREE(s)     __mingw_aligned_free(s)
#endif

#define ASSERT_ALIGNED(p,a)   assert( ((unsigned long)(p) & (a-1))==0 );


template<> union vec<float,4>
{
	ALIGN(16) __m128 m;
    float array[4];
    struct { float x,y,z,w; };
    struct { float r,g,b,a; }; 
    struct { float s,t,p,q; };
	struct { vec3 xyz; float w; };
    struct { vec3 rgb; float a; };
    struct { vec3 stp; float q; };
    struct { float x; vec3 yzw; };
    struct { float r; vec3 gba; };
    struct { float s; vec3 tpq; };
    struct { float x; vec2 yz; float w; };
    struct { float r; vec2 gb; float a; };
    struct { float s; vec2 tp; float q; }; 
    struct { vec2 xy; vec2 zw; };
    struct { vec2 rg; vec2 ba; };
    struct { vec2 st; vec2 pq; };

	inline float& operator[](int i) { return array[i]; }
	inline const float& operator[](int i) const { return array[i]; }

	inline vec() {}
	inline vec( const float& s1 ) { x = y = z = w = s1; }
	inline vec( const float& s1, const float& s2, const float& s3=0, const float& s4=0 ) { x = s1; y = s2; z = s3; w = s4; }
	inline vec( const vec3& v, const float s1 ) { x = v[0]; y = v[1]; z = v[2]; w = s1; }
	inline vec( const vec2& v, const float s1, const float s2 ) { x = v[0]; y = v[1]; z = s1; w = s2; }

	static void* operator new( size_t size )		{ return ALIGNED_MALLOC(size,16); }
    static void* operator new[]( size_t size )		{ return ALIGNED_MALLOC(size,16); }
	static void  operator delete( void* pointer )	{ ALIGNED_FREE(pointer); }
    static void  operator delete[]( void* pointer )	{ ALIGNED_FREE(pointer); }
};

// http://cboard.cprogramming.com/cplusplus-programming/126493-matrix-operations-using-objects.html
/* DISABLED FOR DEBUGGING
inline vec4 operator-( const vec4& a )					{ vec4 t,u; u.m = _mm_set_ps(0,0,0,0); t.m = _mm_sub_ps(u.m,a.m); return t; }

inline vec4 operator+( const vec4& a, const vec4& b )	{ vec4 t; t.m = _mm_add_ps(a.m,b.m); return t; }
inline vec4 operator-( const vec4& a, const vec4& b )	{ vec4 t; t.m = _mm_sub_ps(a.m,b.m); return t; }
inline vec4 operator*( const vec4& a, const vec4& b )	{ vec4 t; t.m = _mm_mul_ps(a.m,b.m); return t; }
inline vec4 operator/( const vec4& a, const vec4& b )	{ vec4 t; t.m = _mm_div_ps(a.m,b.m); return t; }

inline vec4 operator+( const vec4& a, const float& s )	{ vec4 t,u; u.m = _mm_set_ps(s,s,s,s); t.m = _mm_add_ps(a.m,u.m); return t; }
inline vec4 operator-( const vec4& a, const float& s )	{ vec4 t,u; u.m = _mm_set_ps(s,s,s,s); t.m = _mm_sub_ps(a.m,u.m); return t; }
inline vec4 operator*( const vec4& a, const float& s )	{ vec4 t,u; u.m = _mm_set_ps(s,s,s,s); t.m = _mm_mul_ps(a.m,u.m); return t; }
inline vec4 operator/( const vec4& a, const float& s )	{ vec4 t,u; u.m = _mm_set_ps(s,s,s,s); t.m = _mm_div_ps(a.m,u.m); return t; }
 
inline vec4 operator+( const float& s, const vec4& a )	{ return a+s; }
inline vec4 operator-( const float& s, const vec4& a )	{ vec4 t,u; u.m = _mm_set_ps(s,s,s,s); t.m = _mm_sub_ps(u.m,a.m); return t; }
inline vec4 operator*( const float& s, const vec4& a )	{ return a*s; }
inline vec4 operator/( const float& s, const vec4& a )	{ vec4 t,u; u.m = _mm_set_ps(s,s,s,s); t.m = _mm_div_ps(u.m,a.m); return t; }
*/
/*
inline float dot( const vec4& a, const vec4& b )
{
// http://www.devmaster.net/forums/showthread.php?t=14569
	T t = 0;
	for( int i=0; i<N; i++ ) t += a[i]*b[i];
	return t;
}
*/
// http://www.infinity-universe.com/opengl/sselib.h
// MMX normalization: http://www.gamedev.net/community/forums/topic.asp?topic_id=389251

// TODO normalize, sqrt, dot, etc.
// TODO lerp

template<> union mat<float,4>
{
	ALIGN(16) __m128 m[4];
	float array[16]; 
	struct { vec4 r0,r1,r2,r3; }; // in realtà sono colonne !?

	inline float& operator[](int i) { return array[i]; }
	inline const float& operator[](int i) const { return array[i]; }

	inline mat() {}

	inline mat( const float& s )
	{
		set( s,0,0,0, 0,s,0,0, 0,0,s,0, 0,0,0,s );
	}

	inline mat( const float& m0,  const float& m1,  const float& m2,  const float& m3, 
				const float& m4,  const float& m5,  const float& m6,  const float& m7, 
				const float& m8,  const float& m9,  const float& m10, const float& m11, 
				const float& m12, const float& m13, const float& m14, const float& m15 )
	{
		set( m0,m1,m2,m3, m4,m5,m6,m7, m8,m9,m10,m11, m12,m13,m14,m15 );
	}

	inline void set( const float& m0,  const float& m1,  const float& m2,  const float& m3, 
					 const float& m4,  const float& m5,  const float& m6,  const float& m7, 
					 const float& m8,  const float& m9,  const float& m10, const float& m11, 
					 const float& m12, const float& m13, const float& m14, const float& m15 )
	{
		m[0] = _mm_set_ps(m0,m1,m2,m3);
		m[1] = _mm_set_ps(m4,m5,m6,m7);
		m[2] = _mm_set_ps(m8,m9,m10,m11);
		m[3] = _mm_set_ps(m12,m13,m14,m15);
	}

	static void* operator new( size_t size )		{ return ALIGNED_MALLOC(size,16); }
    static void* operator new[]( size_t size )		{ return ALIGNED_MALLOC(size,16); }
	static void  operator delete( void* pointer )	{ ALIGNED_FREE(pointer); }
    static void  operator delete[]( void* pointer )	{ ALIGNED_FREE(pointer); }
};
/*
// http://cs.helsinki.fi/u/ilmarihe/mmul.cpp

inline vec4 operator*( const mat4& m, const vec4& v )
{
	vec4 v0, v1, v2, v3, temp;
	v0.m = _mm_set_ps( v.x, v.x, v.x, v.x );
	v1.m = _mm_set_ps( v.y, v.y, v.y, v.y );
	v2.m = _mm_set_ps( v.z, v.z, v.z, v.z );
	v3.m = _mm_set_ps( v.w, v.w, v.w, v.w );
	v0.m = _mm_mul_ps( v0.m, m.m[0] );
	v1.m = _mm_mul_ps( v1.m, m.m[1] );
	v2.m = _mm_mul_ps( v2.m, m.m[2] );
	v3.m = _mm_mul_ps( v3.m, m.m[3] );
	temp.m = _mm_add_ps( v0.m, v1.m );
	temp.m = _mm_add_ps( temp.m, v2.m );
	temp.m = _mm_add_ps( temp.m, v3.m );
	return temp;
}
*/
/*
inline mat4 transpose( const mat4& a )
{
	mat4 temp(a);
	_MM_TRANSPOSE4_PS( temp.r0.m, temp.r1.m, temp.r2.m, temp.r3.m );
	return temp;
}
*/
// TODO more functions


#endif // __SSE__

inline mat4 projection( vec2 viewport, double fov, double near_plane, double far_plane )
{
	double yf = 1 / tan( fov * PI/360 );
    double xf = yf * viewport.width / viewport.height;
	double f0 = (near_plane + far_plane) / (near_plane - far_plane);
	double f1 = (2*near_plane*far_plane) / (near_plane - far_plane);

	return mat4(   yf,  0,  0,  0,
					0, xf,  0,  0,
					0,  0, f0, -1,
					0,  0, f1,  0   );
}

inline mat4 inverseProjection( mat4& p )
{
	return mat4(   1/p[0],        0,      0,            0,
					    0,   1/p[5],      0,  		    0,
					    0,        0,      0,      1/p[14],
					    0,        0,     -1,  p[10]/p[14] );
}

// TODO det( const mat4 )
// TODO inverse( const mat4 )

inline mat4 transform( const vec3& rotation, const vec3& scale, const vec3& translation )
{
	mat3 r = rotate( rotation );

/*	return mat4(r[0], r[1], r[2], translation.x,
				r[3], r[4], r[5], translation.y,
				r[6], r[7], r[8], translation.z,
				   0,    0,    0, 1);
*/
	return mat4(r[0]*scale.x, r[3]*scale.x, r[6]*scale.x, 0,
				r[1]*scale.y, r[4]*scale.y, r[7]*scale.y, 0,
				r[2]*scale.z, r[5]*scale.z, r[8]*scale.z, 0,
				translation.x, translation.y, translation.z, 1);
}

// RotationAroundAxis
// TODO optimize
inline mat4 rotate4( float angle, float x, float y, float z )
{
	mat4 m;
	
	// normalize
	float mag = float(sqrt( x*x + y*y + z*z ));
	if( mag == 0.0f ) 
	{
		return mat4(1);
	}
	x /= mag;
	y /= mag;
	z /= mag;

	angle *= PI/180.0f;
	float s = float(sin(angle));
	float c = float(cos(angle));

	float xx = x * x;
	float yy = y * y;
	float zz = z * z;
	float xy = x * y;
	float yz = y * z;
	float zx = z * x;
	float xs = x * s;
	float ys = y * s;
	float zs = z * s;
	float c1 = 1.0f - c;

	#define M(row,col)  m[col*4+row]

	M(0,0) = (c1 * xx) + c;
	M(0,1) = (c1 * xy) - zs;
	M(0,2) = (c1 * zx) + ys;
	M(0,3) = 0.0f;

	M(1,0) = (c1 * xy) + zs;
	M(1,1) = (c1 * yy) + c;
	M(1,2) = (c1 * yz) - xs;
	M(1,3) = 0.0f;

	M(2,0) = (c1 * zx) - ys;
	M(2,1) = (c1 * yz) + xs;
	M(2,2) = (c1 * zz) + c;
	M(2,3) = 0.0f;

	M(3,0) = 0.0f;
	M(3,1) = 0.0f;
	M(3,2) = 0.0f;
	M(3,3) = 1.0f;

    #undef M
	
	return m;
}

// EulerToRotation
inline mat4 RotationMatrix( vec3 rotation )
{
	mat4 m;

	float A  = cos( rotation.pitch * PI/180.0f );
	float B  = sin( rotation.pitch * PI/180.0f );
	float C  = cos( rotation.yaw * PI/180.0f );
	float D  = sin( rotation.yaw * PI/180.0f );
	float E  = cos( rotation.roll * PI/180.0f );
	float F  = sin( rotation.roll * PI/180.0f );
	float AD = A * D;
	float BD = B * D;

	m[0]  =   C * E;
	m[1]  = -BD * E + A * F;
	m[2]  =  AD * E + B * F;
	m[3]  = 0.0;
	m[4]  =  -C * F;
	m[5]  =  BD * F + A * E;
	m[6]  = -AD * F + B * E;
	m[7]  = 0.0;
	m[8]  =  -D;
	m[9]  =  -B * C;
	m[10] =   A * C;
	m[11] = 0.0;
	m[12] = 0.0;
	m[13] = 0.0;
	m[14] = 0.0;
	m[15] = 1.0;
	
	return m;
}

// (angles,scale,translate) => matrix
inline mat4 TransformationMatrix( const mat3& transform )
{ 
	mat4 t = mat4(transform.scale);
	t.translation.xyz = transform.scale * transform.position;
	return RotationMatrix(transform.rotation) * t;
}


// RGBA 
inline rgba average( const rgba& color1, const rgba& color2 )
{ 
	rgba temp;
	temp.r = (int(color1.r) + int(color2.r))/int(2);
	temp.g = (int(color1.g) + int(color2.g))/int(2);
	temp.b = (int(color1.b) + int(color2.b))/int(2);
	temp.a = (int(color1.a) + int(color2.a))/int(2);
	return temp;
}

inline rgba average( const rgba& color1, const rgba& color2, const rgba& color3, const rgba& color4 )
{ 
	rgba temp;
	temp.r = (int(color1.r) + int(color2.r) + int(color3.r) + int(color4.r))/int(4);
	temp.g = (int(color1.g) + int(color2.g) + int(color3.g) + int(color4.g))/int(4);
	temp.b = (int(color1.b) + int(color2.b) + int(color3.b) + int(color4.b))/int(4);
	temp.a = (int(color1.a) + int(color2.a) + int(color3.a) + int(color4.a))/int(4);
	return temp;
}

//
// Packing
//

/*
//http://www.gamedev.net/topic/442138-packing-a-float-into-a-a8r8g8b8-texture-shader/
float pack( const vec4& c )
{
	const float B = 64.0;
	const vec4 shift = vec4(1/(B*B*B), 1/(B*B), 1/B, 1);
	float x = dot( c, shift );
	return x;
}
float pack( const rgba& c )
{
	vec4 t = vec4( float(c.r), float(c.g), float(c.b), float(c.a) );///255.0f;
	return pack( t ); 
}
vec4 unpack( const float x )
{
	const float B = 64.0;
	const vec4 shift = vec4(B*B*B, B*B, 1/B, 1);
	const vec4 mask  = vec4(0, 1/B, 1/B, 1/B);
	vec4 c = x * shift;
	c = fract(c);
	c -= vec4(c.x, c.x, c.y, c.z) * mask;  // c.xxyz * mask;
	return c;
}
*/

// TODO use exponent bits too
float pack2( const vec4& c )
{
	const float B = 64.0;
	vec4 t = floor( c * (B-1) );
	return dot(t, vec4(B*B*B, B*B, B, 1));
}
float pack2( const rgba& c )
{
	vec4 t = vec4( float(c.r), float(c.g), float(c.b), float(c.a) ) / 255.0f;
	return pack2( t ); 
}
vec4 unpack2( float x )
{
	const float B = 64.0;
	vec4 t = fract( x * vec4(1/(B*B*B*B), 1/(B*B*B), 1/(B*B), 1/B) );
	return t * B/(B-1);
}

float packNormal( const vec3& n )
{
	float BASE = 256.0f;
	return dot(vec4(floor((n + 1.0f)*0.5f*(BASE-1)),0.0f), vec4(1.0f, BASE, BASE*BASE, BASE*BASE*BASE))/(BASE*BASE*BASE*BASE);;
}

vec4 unpackNormal( float f )
{
	float BASE = 256.0f;
	return (mod(f*(BASE*BASE*BASE*BASE)/vec4(1.0f, BASE, BASE*BASE, BASE*BASE*BASE), BASE)/(BASE-1))*2.0f - 1.0f;
}



// Clipmap stuff
inline float scrollValue( double x, double y )
{
	return float(fract(x/y)*y);
}

// TODO shall return int
inline float tileValue( double x, double y )
{
	return float(floor(x/y));
}



//
// Misc
//
//#include <stdio.h>

template<int N> inline char* str(vec<float,N> v)
{
    static char buffer[256];
    buffer[0] = 0;
    char* p = buffer;
	for( int i = 0; i < N; i++ )
		p += sprintf( p, "%+.2f ", v[i] );
	return buffer;
}

template<int N> inline char* str(vec<double,N> v)
{
    static char buffer[256];
    buffer[0] = 0;
    char* p = buffer;
	for(int i = 0; i < N; i++)
		p += sprintf(p, "%+.2f ", v[i]);
	return buffer;
}


// TEMP just to prove the side effect between calls
template<int N> inline void STR(vec<double,N> v, char* buffer)
{
    buffer[0] = 0;
    char* p = buffer;
	for(int i = 0; i < N; i++)
		p += sprintf(p, "%+.2f ", v[i]);
}

template<int N> inline char* str(vec<int,N> v)
{
    static char buffer[256];
    buffer[0] = 0;
    char* p = buffer;
	for( int i = 0; i < N; i++ )
		p += sprintf(p, "%+i ", v[i]);
	return buffer;
}

template<int N,int M> inline char* str( mat<float,N,M> a )
{
    static char buffer[256];
    buffer[0] = 0;
    char* p = buffer;
	for( int j=0; j<M; j++ )
	{
		for( int i=0; i<N; i++ )
			p += sprintf( p, "%+.2f ", a[j*N+i] );
		//p += sprintf( p, "\n" );
	}
	return buffer;
}


template<class T, int N> inline T max( vec<T,N>& x )
{
	T t = T(-1E37);
	for( int i=0; i<N; i++ ) if( x[i]>t ) t = x[i];
	return t;
}

/*
void testPacking()
{
	
//	vec4 x = vec4( 1.10f, 1.17f, -0.79f, 0.52f );
//	printf( "abs(x) = %s\n", str( abs(x) ) );
//	printf( "max(abs(x)) = %f\n", max(abs(x)) );
//	exit_on_error("...");


	float tolerance = 0.016;
	float BASE = 16.0f;
	
	for( float i=0.0f; i<BASE*BASE*BASE*BASE; i++ )
	{
		vec4 c1 = floor(mod(i/vec4(1.0f, BASE, BASE*BASE, BASE*BASE*BASE), BASE)) / (BASE-1);

		rgba c;
		c.r = c1.r * 255;
		c.g = c1.g * 255;
		c.b = c1.b * 255;
		c.a = c1.a * 255;
		

		float f = pack2( c );
		vec4 c2 = unpack2( f );

		printf( "[%s] [%.2f] [%s] ", str(c1), f, str(c2) );

		vec4 n1 = c1*2.0f - 1.0f;
		f = packNormal( n1 );
		vec4 n2 = unpackNormal( f );

		
	
		if( max(abs(c1-c2))>tolerance || max(abs(n1-n2).xyz)>tolerance )
		{
			exit_on_error("failed");
		}

		printf("\n");
	}
	exit_on_error("done");
}
*/


/*
//#include <stdio.h>
//#include <assert.h>

void test_glsl()
{
	vec4 diffuse = vec4( 0.1, 0.2, 0.3, 0.4 );
    ASSERT( diffuse.x==0.1f && diffuse.y==0.2f && diffuse.z==0.3f && diffuse.w==0.4f );
    ASSERT( diffuse.r==0.1f && diffuse.g==0.2f && diffuse.b==0.3f && diffuse.a==0.4f );    
    ASSERT( diffuse.s==0.1f && diffuse.t==0.2f && diffuse.p==0.3f && diffuse.q==0.4f );

//#ifdef _MSC_VER 
//	ASSERT( diffuse.xy==vec2(0.1,0.2) && diffuse.zw==vec2(0.3,0.4) );
//	ASSERT( diffuse.rgb==vec3(0.1,0.2,0.3) && diffuse.tpq==vec3(0.2,0.3,0.4) );
//	diffuse.yz += 1.0f;
//	ASSERT( diffuse==vec4(0.1, 1.2, 1.3, 0.4) );
//#endif

	vec3 normal = vec3( -0.7, 0.3, 1.1 );
	vec3 light = vec3( 12.5, 6.4, 9.3 );

	mat3 m = rot(30,45,80);
// TODO verify that Rt*R = I
//	mat3 i = transpose(m) * m;
//	ASSERT( nearEqual(i,mat3(1)) );

	normal = m * normal;
	light = m * light;

	ASSERT( equal(normal,vec3(-0.654860, 1.092511, 0.409364)) );
	ASSERT( equal(light, vec3(-0.584514, -3.414564, 16.483297)) );

	normal = normalize(normal);
	light = normalize(light);

	ASSERT( equal(length(normal),1) );
	ASSERT( equal(length(light),1) );

	float NdotL = dot(light, normal);
	ASSERT( equal(NdotL,0.150877) );

	vec4 color = diffuse * NdotL;
	ASSERT( equal(color, vec4(0.015088, 0.181052, 0.196140, 0.060351)) );

	mat4 a = mat4( 1.1, 1.2, 1.3, 1.4, 
				   2.1, 2.2, 2.3, 2.4,
				   3.1, 3.2, 3.3, 3.4,
				   4.1, 4.2, 4.3, 4.4);
	vec4 t = vec4( 6.1, -3.2, 0.1, -2.0 );
	vec4 result = a * t;
	result = -result;
	ASSERT( equal( result, vec4(7.9,7.8,7.7,7.6)) );
}
*/
/*
#include "os/timer.h"

void test_sse()
{
	svec4 v;
	ASSERT_ALIGNED( v.array, 16 );

const int size = 1024;

	for(int i=0; i<size; i++)
	{
		svec4* t = new svec4[11];
		ASSERT_ALIGNED( t[0].array, 16 );
		delete[] t;
	}

	svec4 v1[size];
	ASSERT_ALIGNED( v1[0].array, 16 );

	svec4* v2 = new svec4[size];
	svec4* v3 = new svec4[size];

	Timer timer;

for(int t=0; t<1000; t++)
{
	svec4* v1_ptr = v1;
	svec4* v2_ptr = v2;
	for(int i=0; i<size; i++) 
	{
		*v1_ptr++ = svec4(1,2,3,4);
		*v2_ptr++ = svec4(5,6,7,8);
	}

	v1_ptr = v1;
	v2_ptr = v2;
	svec4* v3_ptr = v3;
	for(int i=0; i<size; i++)
	{
		*v3_ptr++ = (*v1_ptr++) + (*v2_ptr++);
	}

	v3_ptr = v3;
	for(int i=0; i<size; i++)
	{
//		ASSERT( (*v3_ptr++) == svec4(6,8,10,12) );
	}
}

	printf("elapsed %f\n", timer.delta() );

	delete[] v2;
	delete[] v3;
}
*/
/*
void test_sum_the_other_way()
{
	printf( "test_sum_the_other_way\n" );
	const int n = 100000000;
	
	vec4 a(1),b(2),c(3),d(4);

	//ASSERT_ALIGNED(&a.m,16);
	//ASSERT_ALIGNED(&b.m,16);
	//ASSERT_ALIGNED(&c.m,16);
	//ASSERT_ALIGNED(&d.m,16);

	Timer timer;
	for(int i=0; i<n/10; i++)
	{
		a = fract( atan(b,c) );
		b = clamp( a,b,c );
	}
	printf( "fract atan clamp: elapsed %f\n", timer.delta() );

	for(int i=0; i<n; i++)
	{
		a = b - c;
		b = c - d;
		c = d + a;
		d = a + b;
	}
	printf( "[+,-]: elapsed %f\n", timer.delta() );

	for(int i=0; i<n; i++)
	{
		a = b + c;
		b = c + d;
		c = d + a;
		d = a + b;
	}
	printf( "[+]: elapsed %f\n", timer.delta() );

	for(int i=0; i<n; i++)
	{
		a += b;
		b += c;
		c += d;
		d += a;
	}
	printf( "[+=]: elapsed %f\n", timer.delta() );

	//printf(" results: %s, %s, %s, %s\n", str(x), str(y), str(z), str(w) );

	
	getc(stdin);
	exit(0);
}
*/

