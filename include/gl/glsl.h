#pragma once

#include <math.h>
#undef min
#undef max

#if defined(_MSC_VER) 
    // Disable warnings for possible loss of data. 
	#pragma warning(disable: 4244)
#endif



// ======================================
// BASE TEMPLATE FORWARD DEFINITIONS
// ======================================

template<class T, int N> union vec;
template<class T, int N, int M=N> union mat;


// ======================================
// GLSL TYPES
// ======================================

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



// ======================================
// SCALARS
// ======================================

const double PI			= 3.141592653589793238462;
const double LN2		= 0.693147180559945309417;
const double EPSILON	= 0.0002;//0.000001;

// Functions
template<class T> inline T radians(T x)		{ return (PI/180)*x; }
template<class T> inline T degrees(T x)		{ return (180/PI)*x; }
template<class T> inline T sin(T x)			{ return T(::sin(x)); }
template<class T> inline T cos (T x)		{ return T(::cos(x)); }
template<class T> inline T tan (T x)		{ return T(::tan(x)); }
template<class T> inline T asin(T x)		{ return T(::asin(x)); }
template<class T> inline T acos(T x)		{ return T(::acos(x)); }
template<class T> inline T atan(T x)		{ return T(::atan(x)); }
template<class T> inline T atan(T x, T y)	{ return T(::atan(x,y)); }
template<class T> inline T sinh(T x, T y)	{ return T(::sinh(x,y)); }
template<class T> inline T cosh(T x, T y)	{ return T(::cosh(x,y)); }
template<class T> inline T tanh(T x, T y)	{ return T(::tanh(x,y)); }
template<class T> inline T asinh(T x, T y)	{ return T(::asinh(x,y)); }
template<class T> inline T acosh(T x, T y)	{ return T(::acosh(x,y)); }
template<class T> inline T atanh(T x, T y)	{ return T(::atanh(x,y)); }
template<class T> inline T pow (T x, T y)	{ return T(::pow(x,y)); }
template<class T> inline T exp (T x)		    { return T(::exp(x)); }
template<class T> inline T log (T x)		    { return T(::log(x)); }
template<class T> inline T exp2(T x)		    { return T(::exp(x*LN2)); }
template<class T> inline T log2(T x)		    { return T(::log(x)/LN2); }
template<class T> inline T sqrt(T x)		    { return T(::sqrt(x)); }
template<class T> inline T inversesqrt(T x)	{ return T(1/::sqrt(x)); }
template<class T> inline T abs(T x)			{ return x<0 ? -x : x; }
template<class T> inline T sign(T x)		    { return x>0 ? 1.0 : x<0 ? -1.0 : 0.0; }
template<class T> inline T floor(T x)		{ return T(::floor(x)); }
template<class T> inline T ceil(T x)		    { return T(::ceil(x)); }
template<class T> inline T fract(T x)        { return T(x - ::floor(x)); }
template<class T> inline T mod(T x, T y)     { return T(::fmod(x,y)); }
template<class T> inline T min(T x, T y)	    { return x < y ? x : y; }
template<class T> inline T max(T x, T y)	    { return x > y ? x : y; }
template<class T> inline T clamp(T x, T x0 = T(0), T x1 = T(1)) { return x<x0 ? x0 : x>x1 ? x1 : x; }
template<class T> inline T mix(T x0, T x1, float a) { return x0*(1-a) + x1*a; }
template<class T> inline T step(T x0, T x)   { return x<=x0 ? 0.0f : 1.0f; }
template<class T> inline T smoothstep(T x0, T x1, T x)  { T t = clamp((x-x0)/(x1-x0), T(0), T(1)); return t*t*(3 - 2*t); }
template<class T> inline bool equal(T a, T b, T epsilon=T(EPSILON)) { return abs(a-b)<epsilon; }

// Geometric functions
template<class T> inline T length(T x)       { return ::sqrt(x*x); }
template<class T> inline T distance(T p0, T p1) { return length(p0-p1); }
template<class T> inline T dot(T x, T y)     { return x*y; }
template<class T> inline T normalize(T)      { return 1; }
template<class T> inline T faceforward(T n, T i, T nref) { return nref*i < 0 ? n : -n; }
template<class T> inline T reflect(T i, T n) { return i - 2*dot(n,i)*n; }
/*
template<class T> inline T refract(T i, T n, double eta) { 
    float k = 1.0 - eta*eta*(1.0 - dot(n,i)*dot(n,i)); 
    return k < 0.0 ? 0.0 : eta*i - (eta*dot(n,i) + sqrt(k))*n; 
}
*/ /*
float refract(float i, float n, double eta) { 
    float k = 1.0 - eta*eta*(1.0 - dot(n,i)*dot(n,i)); 
    return k < 0.0 ? 0.0 : eta*i - (eta*dot(n,i) + sqrt(k))*n; 
}*/





// ======================================
// VECTORS
// ======================================

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

// Component-wise operations
#define scalar const T& 
#define var(type) type<T,N>
#define out(type) template<class T, int N> inline type<T,N>
#define in(type) const type<T,N>&
#define foreach(f) {var(vec) t; for(int i=0; i<N; i++) (t)[i] = (f); return t;}

out(vec) operator+(in(vec) a)                   { return a; }
out(vec) operator-(in(vec) a)                   foreach(-a[i])
out(vec) operator+(in(vec) a, in(vec) b)        foreach(a[i]+b[i])
out(vec) operator-(in(vec) a, in(vec) b)        foreach(a[i]-b[i])
out(vec) operator*(in(vec) a, in(vec) b)        foreach(a[i]*b[i])
out(vec) operator/(in(vec) a, in(vec) b)        foreach(a[i]/b[i])

out(vec) operator+(in(vec) a, scalar s)         foreach(a[i]+s)
out(vec) operator-(in(vec) a, scalar s)         foreach(a[i]-s)
out(vec) operator*(in(vec) a, scalar s)         foreach(a[i]*s)
out(vec) operator/(in(vec) a, scalar s)         foreach(a[i]/s)

out(vec) operator+(scalar s, in(vec) a)         foreach(s+a[i])
out(vec) operator-(scalar s, in(vec) a)         foreach(s-a[i])
out(vec) operator*(scalar s, in(vec) a)         foreach(s*a[i])
out(vec) operator/(scalar s, in(vec) a)         foreach(s/a[i])

out(vec)& operator+=(var(vec)& a, in(vec) b)    {a = a + b; return a;}
out(vec)& operator-=(var(vec)& a, in(vec) b)    {a = a - b; return a;}
out(vec)& operator*=(var(vec)& a, in(vec) b)    {a = a * b; return a;}
out(vec)& operator/=(var(vec)& a, in(vec) b)    {a = a / b; return a;}

out(vec)& operator+=(var(vec)& a, scalar s)   {a = a + s; return a;}
out(vec)& operator-=(var(vec)& a, scalar s)   {a = a - s; return a;}
out(vec)& operator*=(var(vec)& a, scalar s)   {a = a * s; return a;}
out(vec)& operator/=(var(vec)& a, scalar s)   {a = a / s; return a;}

out(vec) radians(in(vec) a)               foreach(radians(a[i]))
out(vec) degrees(in(vec) a)               foreach(degrees(a[i]))
out(vec) sin(in(vec) a)                   foreach(sin(a[i]))
out(vec) cos(in(vec) a)                   foreach(cos(a[i]))
out(vec) tan(in(vec) a)                   foreach(tan(a[i]))
out(vec) asin(in(vec) a)                  foreach(asin(a[i]))
out(vec) acos(in(vec) a)                  foreach(acos(a[i]))
out(vec) atan(in(vec) a)                  foreach(atan(a[i]))
out(vec) atan(in(vec) y, in(vec) x)       foreach(atan2(y[i],x[i]))
out(vec) atan(in(vec) y, scalar s)        foreach(atan2(y[i],s))
out(vec) pow(in(vec) y, in(vec) x)        foreach(pow(y[i],x[i]))
out(vec) pow(in(vec) y, scalar s)         foreach(pow(y[i],s))
out(vec) exp(in(vec) a)                   foreach(exp(a[i]))
out(vec) log(in(vec) a)                   foreach(log(a[i]))
out(vec) exp2(in(vec) a)                  foreach(exp2(a[i]))
out(vec) log2(in(vec) a)                  foreach(log2(a[i]))
out(vec) sqrt(in(vec) a)                  foreach(sqrt(a[i]))
out(vec) inversesqrt(in(vec) a)           foreach(inversesqrt(a[i]))
out(vec) abs(in(vec) a)                   foreach(abs(a[i]))
out(vec) sign(in(vec) a)                  foreach(sign(a[i]))
out(vec) floor(in(vec) a)                 foreach(floor(a[i]))
out(vec) ceil(in(vec) a)                  foreach(ceil(a[i]))
out(vec) fract(in(vec) a)                 foreach(fract(a[i]))
out(vec) mod(in(vec) a, scalar b)         foreach(mod(a[i],b))
out(vec) mod(scalar a, in(vec) b)         foreach(mod(a,b[i]))
out(vec) mod(in(vec) a, in(vec) b)        foreach(mod(a[i],b[i]))
out(vec) min(in(vec) a, scalar b)         foreach(min(a[i],b))
out(vec) min(scalar a, in(vec) b)         foreach(min(a,b[i]))
out(vec) min(in(vec) a, in(vec) b)        foreach(min(a[i],b[i]))
out(vec) max(in(vec) a, scalar b)         foreach(max(a[i],b))
out(vec) max(scalar a, in(vec) b)         foreach(max(a,b[i]))
out(vec) max(in(vec) a, in(vec) b)        foreach(max(a[i],b[i]))
out(vec) clamp(in(vec) a, scalar b, scalar c) foreach(clamp(a[i],b,c))
out(vec) clamp(in(vec) a, in(vec) b, in(vec) c) foreach(clamp(a[i],b[i],c[i]))
out(vec) mix(in(vec) a, in(vec) b, scalar c) foreach(mix(a[i],b[i],c))
out(vec) mix(in(vec) a, in(vec) b, in(vec) c) foreach(mix(a[i],b[i],c[i]))
out(vec) step(in(vec) a, in(vec) b)        foreach(step(a[i],b[i]))
out(vec) smoothstep(in(vec) a, in(vec) b, in(vec) c) foreach(smoothstep(a[i],b[i],c[i]))

// Geometric functions
template<class T, int N> T dot(in(vec) a, in(vec) b)
{
	T t = 0;
	for( int i=0; i<N; i++ ) t += a[i]*b[i];
	return t;
}
template<class T, int N> T length(in(vec) a)                { return T(sqrt(dot(a,a))); }
template<class T, int N> T distance(in(vec) a, in(vec) b)   { return length(a-b); }

out(vec) normalize(in(vec) a)
{
	T len = length(a);
	if (len  == 0) return a;
	return a/len;
}
out(vec) faceforward(in(vec) a, in(vec) b, in(vec) c) { return dot(c,b) < 0 ? a : -a; }
out(vec) reflect(in(vec) a, in(vec) b) { return a - 2*dot(b,a)*b; }
out(vec) refract(in(vec) a, in(vec) b, const float& c)
{
	T d = dot(b,a);
	T k = 1.0f - c*c*( 1 - d*d );
	return k<0.0f ? T(0.0f) : c*a - (c*d + sqrt(k))*b;
}
template<class T, int N> bool operator==(in(vec) a, in(vec) b) { return distance(a,b) < EPSILON; }
template<class T, int N> bool operator!=(in(vec) a, in(vec) b) { return distance(a,b) > EPSILON; }
template<class T, int N> bool equal(in(vec) a, in(vec) b, const T& epsilon = T(EPSILON)) { return distance(a,b) < epsilon; }

// Partial specializations for T==bool
template<class T, int N> inline vec<bool,N> lessThan( in(vec) a, in(vec) b )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = a[i] < b[i];
	return t;
}
template<class T, int N> inline vec<bool,N> lessThanEqual( in(vec) a, in(vec) b )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = a[i] <= b[i];
	return t;
}
template<class T, int N> inline vec<bool,N> greaterThan( in(vec) a, in(vec) b )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = a[i] > b[i];
	return t;
}
template<class T, int N> inline vec<bool,N> greaterThanEqual( in(vec) a, in(vec) b )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = a[i] >= b[i];
	return t;
}
template<class T, int N> inline vec<bool,N> equal( in(vec) a, in(vec) b )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = equal(a[i],b[i]);
	return t;
}
template<class T, int N> inline vec<bool,N> notEqual( in(vec) a, in(vec) b )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = !equal(a[i],b[i]);
	return t;
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
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = !v[i];
	return t;
}
#endif
// this operator!() can be used instead of the above not()
template<int N> inline vec<bool,N> operator!( const vec<bool,N>& v )
{
	vec<bool,N> t;
	for( int i=0; i<N; i++ ) t[i] = !v[i];
	return t;
}

out(vec) operator*( const mat<T,N>& m, in(vec) v )
{
	vec<T,N> t;
	for(int j=0; j<N; j++)
	{
		t[j] = 0;
		for(int i=0; i<N; i++)
			t[j] += m[j+N*i] * v[i];
	}
	return t;
}

out(vec) operator*( in(vec) v, const mat<T,N>& m )
{
	vec<T,N> t;
	for(int j=0; j<N; j++) 
	{
		t[j] = 0;
		for(int i=0; i<N; i++)
			t[j] += m[i+N*j] * v[i];
	}
	return t;
}
out(vec)& operator*=( vec<T,N>& a, const mat<T,N>& m )
{
	a = a * m;
	return a;
}






// ======================================
// MATRICES
// ======================================

template<class T, int N, int M> union mat
{
	T array[N*M];

	inline T& operator[](int i)	{ return array[i]; }
	inline const T& operator[](int i) const { return array[i]; }

	inline mat() {}
	inline mat( const T& s )		{ diag( vec<T,N>(s) ); }
	inline mat( in(vec) v ) { diag(v); }

	inline mat( const mat<T,N-1>& m ) 
	{ 
		for( int j=0; j<N-1; j++ )
		for( int i=0; i<N-1; i++ ) array[j*N+i] = m[j*N+i];

		for( int j=0; j<N-1; j++ ) array[j*N + N-1] = 0;
		for( int i=0; i<N-1; i++ ) array[(N-1)*N + i] = 0; 
		array[N*N-1] = T(1);
	} 
};

#undef out
#undef in
#undef foreach
#define out(mat) template<class T, int N, int M> inline mat<T,N,M>
#define in(mat) const mat<T,N,M>&
#define scalar const T&
#define foreach(f) {mat<T,N,M> t; for(int i=0; i<N*M; i++) (t)[i] = (f); return t;}
out(mat) operator-(in(mat) a) foreach(-a[i])
out(mat) operator+(in(mat) a, in(mat) b) foreach(a[i]+b[i])
out(mat) operator-(in(mat) a, in(mat) b) foreach(a[i]-b[i])
out(mat) matrixCompMult(in(mat) a, in(mat) b) foreach(a[i]*b[i])
out(mat) operator/(in(mat) a, in(mat) b) foreach(a[i]/b[i])
out(mat) operator+(in(mat) a, scalar s) foreach(a[i]+s)
out(mat) operator-(in(mat) a, scalar s) foreach(a[i]-s)
out(mat) operator*(in(mat) a, scalar s) foreach(a[i]*s)
out(mat) operator/(in(mat) a, scalar s) foreach(a[i]/s)
out(mat) operator+(scalar s, in(mat) a) foreach(s+a[i])
out(mat) operator-(scalar s, in(mat) a) foreach(s-a[i])
out(mat) operator*(scalar s, in(mat) a) foreach(s*a[i])
out(mat) operator/(scalar s, in(mat) a) foreach(s/a[i])
out(mat)& operator+=(mat<T,N,M>& a, in(mat) b)   {a = a + b; return a;}
out(mat)& operator-=(mat<T,N,M>& a, in(mat) b)   {a = a - b; return a;}
out(mat)& operator*=(mat<T,N,M>& a, in(mat) b)   {a = a * b; return a;}
out(mat)& operator/=(mat<T,N,M>& a, in(mat) b)   {a = a / b; return a;}
out(mat)& operator+=(mat<T,N,M>& a, scalar s)   {a = a + s; return a;}
out(mat)& operator-=(mat<T,N,M>& a, scalar s)   {a = a - s; return a;}
out(mat)& operator*=(mat<T,N,M>& a, scalar s)   {a = a * s; return a;}
out(mat)& operator/=(mat<T,N,M>& a, scalar s)   {a = a / s; return a;}

// Geometric operations

template<class T, int N, int M> inline T distance(const mat<T,N,M>& A, const mat<T,N,M>& B)
{
    T result = T(0.0);
    for(int i=0; i<N*M; i++)
    {
        result += (A[i]-B[i])*(A[i]-B[i]);
    }
    return sqrt(result);
}
/*
template<class T, int N, int M> inline T distance2(const mat<T,N,M>& A, const mat<T,N,M>& B)
{
    T result = T(0.0);
    for(int i=0; i<N*M; i++)
    {
        result += abs(A[i]-B[i]);
    }

    return result;
}
*/



// Geometric matrix operations

// matrix-matrix multiplication
// TODO N*M version
template<class T, int N> inline mat<T,N> operator*( const mat<T,N>& a, const mat<T,N>& b )
{
	mat<T,N> t;
	for(int i=0; i<N; i++)
	for(int j=0; j<N; j++)
	{
		t[j+i*N] = 0;
		for(int k=0; k<N; k++)
			t[j+i*N] += a[j+k*N] * b[k+i*N];
	}
	return t;
}


template<class T, int N> inline mat<T,N>& operator*=( mat<T,N>& a, const mat<T,N>& b )
{
	a = a * b;
	return a;
}



template<class T, int N, int M> inline mat<T,M,N> outerProduct( const vec<T,M>& a, const vec<T,N>& b )
{
	mat<T,N,M> t;
	for(int j=0; j<M; j++)
	for(int i=0; i<N; i++)
		t[j*N+i] = a[j]*b[i];
	return t;
}

template<class T, int N, int M> inline mat<T,M,N> transpose( const mat<T,N,M>& a )
{
	mat<T,M,N> t;
	for(int j=0; j<N; j++ )
	for(int i=0; i<M; i++ )
		t[j*M+i] = a[i*N+j];
	return t;
}
/*
template<class T, int N, int M> inline mat<T,M,N> determinant( const mat<T,N,M>& a )
{
	T t;
	// TODO
	return t;
}

// https://chi3x10.wordpress.com/2008/05/28/calculate-matrix-inversion-in-c/
// https://www.tutorialspoint.com/cplusplus-program-to-find-inverse-of-a-graph-matrix
template<class T, int N, int M> inline mat<T,M,N> inverse( const mat<T,N,M>& a )
{
	mat<T,M,N> t;
	// TODO
	return t;
}
*/
template<class T, int N, int M> inline bool equal( const mat<T,N,M>& a, const mat<T,N,M>& b, const T& epsilon=T(EPSILON) )
{
	for( int i=0; i<N*M; i++ ) if( !equal(a[i],b[i],epsilon) ) return false;
	return true;
}




// ======================================
// SPECIALIZATIONS
// ======================================


//
// Partial specializations for N=2
//
template<class T> union vec<T,2>
{
    T array[2];
    struct { T x,y; };
    struct { T r,g; };
    struct { T s,t; };

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
    struct { vec<T,2> rg; T b; };
    struct { vec<T,2> st; T p; };
    struct { T x; vec<T,2> yz; };
    struct { T r; vec<T,2> gb; };
    struct { T s; vec<T,2> tp; };

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
	vec<T,3> t;
	t[0] = a[1]*b[2] - a[2]*b[1];
	t[1] = a[2]*b[0] - a[0]*b[2];
	t[2] = a[0]*b[1] - a[1]*b[0];
	return t;
}

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
	inline vec( const vec<T,3>& v, const T s1 ) { x = v[0]; y = v[1]; z = v[2]; w = s1; }
	inline vec( const vec<T,2>& v, const T s1, const T s2 ) { x = v[0]; y = v[1]; z = s1; w = s2; }
	inline vec( const T s, const vec<T,3>& v  ) { x = s; y = v[0]; z = v[1]; w = v[2]; }
};

template<class T> inline vec<T,4> cross( const vec<T,4>& a, const vec<T,4>& b )
{
	vec<T,4> t;
	t[0] = a[1]*b[2] - a[2]*b[1];
	t[1] = a[2]*b[0] - a[0]*b[2];
	t[2] = a[0]*b[1] - a[1]*b[0];
	t[3] = 1;
	return t;
}




// TODO lookat
// TODO slerp


// ======================================
// MATRIX PARTIAL SPECIALIZATIONS
// ======================================

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






//
// Partial specializations for T=float
//

/* TODO verify 
template<> inline float sqrt(float x)
{
	float t;
	_asm
	{
		movupd xmm0, x;
		sqrtss xmm0, xmm0;
		movupd temp, xmm0;
	}
	return t;
}
*/



//
// Partial specializations for T=float and N=2
//

#if 0
template<> inline vec<float,2> operator*( const vec<float,2>& a, const vec<float,2>& b )
{ 
	vec<float,2> t; 
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
	return t;
}
#endif


// Geometric operations
/*
float distance(mat2 A, mat2 B)
{
    float result = 0.0;
    for(int i=0; i<=3; i++)
    {
        result += (A[i]-B[i])*(A[i]-B[i]);
    }
    return sqrt(result);
}
*/

// TODO transformation factory: translation matrix, rotation matrix, scale matrix

// TODO does it work with both mat2 and dmat2 ?
template<class T> inline mat<T, 2> rotate(T angle)
{
	mat<T, 2> t;

	T A  = (T) cos( angle * PI / (-180.0) );
	T B  = (T) sin( angle * PI / (-180.0) );

	t[0] =  A;
	t[1] = -B;
	t[2] =  B;
	t[3] =  A;

	return t;
}

template<class T> inline T determinant(const mat<T, 2>& m)
{
	return m[0]*m[3] - m[1]*m[2];
}

template<class T> inline mat<T, 2> inverse(const mat<T, 2>& m)
{
	mat<T, 2> t = mat<T, 2>(m[3], -m[1], -m[2], m[0]);
	t /= determinant(m);
	return t;
}



//
// Partial specializations for T=float and N=3
//

/* TODO debug
inline vec3 normalize( const vec3& v )
{
	vec3 t;
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
	return t;
}
*/


template<class T> inline mat<T, 3> rotate( T pitch, T yaw, T roll )
{
	mat<T, 3> t;
	T Cx, Sx, Cy, Sy, Cz, Sz, CxSy, SxSy;

	Cx  = (T) cos( pitch * PI / (-180.0) );
	Sx  = (T) sin( pitch * PI / (-180.0) );
	Cy  = (T) cos(   yaw * PI / (-180.0) );
	Sy  = (T) sin(   yaw * PI / (-180.0) );
	Cz  = (T) cos(  roll * PI / (-180.0) );
	Sz  = (T) sin(  roll * PI / (-180.0) );
	CxSy = Cx * Sy;
	SxSy = Sx * Sy;

	t[0] =    Cy * Cz;
	t[1] =  SxSy * Cz + Cx * Sz;
	t[2] = -CxSy * Cz + Sx * Sz;
	t[3] =   -Cy * Sz;
	t[4] = -SxSy * Sz + Cx * Cz;
	t[5] =  CxSy * Sz + Sx * Cz;
	t[6] =    Sy;
	t[7] =   -Sx * Cy;
	t[8] =    Cx * Cy;

	return t;
}

template<class T> inline mat<T, 3> rotate( const vec<T, 3>& rotation )
{
	return rotate( rotation.x, rotation.y, rotation.z );
}

template<class T> inline mat<T, 3> transform( T rotation, T scale, T translation )
{
	return mat<T, 3>( mat<T, 2>(scale)*rotate(rotation), vec2(translation) );
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
	mat3 t;
	t[0] = + (m[4] * m[8] - m[5] * m[7]);
	t[1] = - (m[1] * m[8] - m[2] * m[7]);
	t[2] = + (m[1] * m[5] - m[2] * m[4]);
	t[3] = - (m[3] * m[8] - m[5] * m[6]);
	t[4] = + (m[0] * m[8] - m[2] * m[6]);
	t[5] = - (m[0] * m[5] - m[2] * m[3]);
	t[6] = + (m[3] * m[7] - m[4] * m[6]);
	t[7] = - (m[0] * m[7] - m[1] * m[6]);
	t[8] = + (m[0] * m[4] - m[1] * m[3]);
	t /= determinant(m);

	return t;
}

inline mat3 inverseTranspose( const mat3& m )
{
	mat3 t;
	t[0] = + (m[4] * m[8] - m[5] * m[7]);
	t[1] = - (m[3] * m[8] - m[5] * m[6]);
	t[2] = + (m[3] * m[7] - m[4] * m[6]);
	t[3] = - (m[1] * m[8] - m[2] * m[7]);
	t[4] = + (m[0] * m[8] - m[2] * m[6]);
	t[5] = - (m[0] * m[7] - m[1] * m[6]);
	t[6] = + (m[1] * m[5] - m[2] * m[4]);
	t[7] = - (m[0] * m[5] - m[2] * m[3]);
	t[8] = + (m[0] * m[4] - m[1] * m[3]);
	t /= determinant(m);

	return t;
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
	vec4 v0, v1, v2, v3, t;
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
	return t;
}
*/
/*
inline mat4 transpose( const mat4& a )
{
	mat4 temp(a);
	_MM_TRANSPOSE4_PS( temp.r0.m, temp.r1.m, temp.r2.m, temp.r3.m );
	return t;
}
*/
// TODO more functions


#endif // __SSE__


// ======================================
// NON STANDARD FUNCTIONS
// ======================================


inline mat4 projection( vec2 viewport, double fov, double near_plane, double far_plane )
{
	double yf = 1 / tan( fov * PI/360 );
    double xf = yf * viewport.x / viewport.y;
	double f0 = (near_plane + far_plane) / (near_plane - far_plane);
	double f1 = (2*near_plane*far_plane) / (near_plane - far_plane);

	return mat4( yf,  0,  0,  0,
			      0, xf,  0,  0,
				  0,  0, f0, -1,
				  0,  0, f1,  0  );
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

	float A  = cos( rotation.x * PI/180.0f );
	float B  = sin( rotation.x * PI/180.0f );
	float C  = cos( rotation.y * PI/180.0f );
	float D  = sin( rotation.y * PI/180.0f );
	float E  = cos( rotation.z * PI/180.0f );
	float F  = sin( rotation.z * PI/180.0f );
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

//
    // Quaternion
    //

    // https://github.com/mattatz/ShibuyaCrowd/blob/master/source/shaders/common/quaternion.glsl
    // https://www.geeks3d.com/20141201/how-to-rotate-a-vertex-by-a-quaternion-in-glsl/

    vec4 qmul(vec4 q1, vec4 q2) {
	    return vec4(
		    q2.xyz * q1.w + q1.xyz * q2.w + cross(q1.xyz, q2.xyz),
		    q1.w * q2.w - dot(q1.xyz, q2.xyz)
	    );
    }
/*    
    vec3 qtransform(vec4 q, vec3 v) { 
	    return v + 2.0*cross(cross(v, q.xyz ) + q.w*v, q.xyz);
	} 
*/

    vec3 rotate_vector(vec3 v, vec4 r) {
	    vec4 r_c = r * vec4(-1, -1, -1, 1);
	    return qmul(r, qmul(vec4(v, 0), r_c)).xyz;
    }

    vec3 rotate_vector_at(vec3 v, vec3 center, vec4 r) {
	    vec3 dir = v - center;
	    return center + rotate_vector(dir, r);
    }

    vec4 rotate_angle_axis(float angle, vec3 axis) {
	    float sn = sin(angle * 0.5);
	    float cs = cos(angle * 0.5);
	    return vec4(axis * sn, cs);
    }

    vec4 conjugate(vec4 q) {
	    return vec4(-q.xyz, q.w);
    }

    vec4 invert(vec4 q) {
        return vec4(-q.xyz, q.w) / sqrt(dot(q, q));
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




// Restore warning for possible loss of data
#if defined(_MSC_VER) 
    //#pragma warning(pop) 
#endif