#pragma once

#include "glsl.h"

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
//typedef vec<float,8>		vec8;		// TEMP terrain
/*
#define PRECISION 16
typedef long long			fixed;
typedef vec<fixed,2>		fpvec2;
typedef vec<fixed,3>		fpvec3;
*/


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

template<class T> inline T amod(T x, T y)		{ return (x/y - floor(x/y))*y; }

//_OUT(vec) amod( _IN(vec) v, T& s )
template<class T, int N> inline vec<T,N> amod( const vec<T,N>& v, T& s )
{
	vec<T,N> t;
	for( int i=0; i<N; i++ ) t[i] = amod(v[i],s);
	return t;
}

//_OUT_S(T) swap(T& a, T& b)  { T t; t=a; a=b; b=t; return t;}
template<class T, int N> inline T swap(T& a, T& b)  { T t; t=a; a=b; b=t; return t;}




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