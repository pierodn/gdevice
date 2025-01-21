#pragma once

#include "type/glsl.h"

#define QUINTIC_INTERPOLATION


//
// Pseudo Random Number Generators
// 

#if 0 // ORIGINAL
inline static float rand( int s ) 
{
	return ((s^37457)*(s^92837)*(s^47563) & 0x0fffffff)/(float)0x10000000;
}
inline static float rand( int x, int y ) 
{
	return rand( y*1234567 + x^7654321 );
}
#elif 1
inline static float rand( int x ) 
{
	return ((x^1234567897)*(x^9876543211) & 0x0fffffff)/(float)0x10000000;
}
inline float rand( int x, int y ) 
{
	return ( ((x^1234567897)*(y^9876543211)) & 0x0fffffff)/(float)0x10000000;
}
inline float rand( vec2 p ) 
{
	return ( ((int(p.x)^1234567897)*(int(p.y)^9876543211)) & 0x0fffffff)/(float)0x10000000;
}
float hash2( vec2 p )
{
	return fract(sin(p.x+p.y)*43758.5453);
}
#else // GLSL
inline static float rand( int x, int y ) 
{
	const float seed = 12345666;
	return 1.0 - 2.0 * fract(sin(dot((cos(vec2(x, y)) * sin(seed)), vec2(15.12345, 91.98765))) * 115309.76543);
}
#endif
/*
inline float rand( int x, int y, int z ) 
{
	return ( ((x^123456789)*(y^876543210)*(z^543212345)) & 0x0fffffff)/(float)0x10000000;
}
*/

//
// Coherent noise functions WITHOUT derivatives
// 

// TODO float zvalue
// TODO float zgradient (perlin or simplex)
// TODO float zcell (voronoi)





//
// Coherent noise functions WITH derivatives
//
#if 0
inline vec3 value( vec2 p ) 
{
	// integer part
	int ix = p.x>0 ? int(p.x) : int(p.x)-1;
	int iy = p.y>0 ? int(p.y) : int(p.y)-1;

	// fractional part
	float fx = p.x - ix;
	float fy = p.y - iy;

#ifdef QUINTIC_INTERPOLATION
	// interpolation
	float u = fx * fx * fx * (fx * (fx * 6.0f - 15.0f) + 10.0f);
	float v = fy * fy * fy * (fy * (fy * 6.0f - 15.0f) + 10.0f);
	// derivatives
	float du = 30.0f * fx * fx * (fx * (fx - 2.0f) + 1.0f);
	float dv = 30.0f * fy * fy * (fy * (fy - 2.0f) + 1.0f);
#else
	// interpolation
	float u = (1-cos(PI*fx))/2;
	float v = (1-cos(PI*fy))/2;
	// derivatives
	float du = (PI*sin(PI*fx))/2;
	float dv = (PI*sin(PI*fy))/2;
#endif

	// lattice points
	float a = rand(ix+0,iy+0);
	float b = rand(ix+1,iy+0);
	float c = rand(ix+0,iy+1);
	float d = rand(ix+1,iy+1);

	// bilinear interpolation
	float k0 = a;
	float k1 = b - a;
	float k2 = c - a;
	float k3 = a - b - c + d;

	vec3 n;
	n.x = du*(k1 + k3*v);					// dNdx
	n.y = dv*(k2 + k3*u);					// dNdy
	n.z = k0 + k1*u + k2*v + k3*u*v;		// N(x,y)

	return n;
}
#else
inline vec3 value( vec2 p ) 
{
	vec2 i = floor(p);
	vec2 f = fract(p);

	float a = rand(i+vec2(0,0));
	float b = rand(i+vec2(1,0));
	float c = rand(i+vec2(0,1));
	float d = rand(i+vec2(1,1));

	vec2 u = f*f*f*(f*(f*6.0f-15.0f)+10.0f);
	vec2 du = f*f*(f*(30.0f*f-60.0f)+30.0f);

	return vec3(
		du * (vec2(b-a,c-a)+(a-b-c+d)*vec2(u.y,u.x)),	
		a + (b-a)*u.x + (c-a)*u.y + (a-b-c+d)*u.x*u.y );
}
#endif
/*
inline vec4 value( const vec3& p ) 
{
	// integer part
	int ix = p.x>0 ? (int)p.x : int(p.x)-1;
	int iy = p.y>0 ? (int)p.y : int(p.y)-1;
	int iz = p.z>0 ? (int)p.z : int(p.z)-1;

	// fractional part
	float fx = p.x - ix;
	float fy = p.y - iy;
	float fz = p.z - iz;

	// quintic interpolation
	float u = fx * fx * fx * (fx * (fx * 6.0f - 15.0f) + 10.0f);
	float v = fy * fy * fy * (fy * (fy * 6.0f - 15.0f) + 10.0f);
	float w = fz * fz * fz * (fz * (fz * 6.0f - 15.0f) + 10.0f);

	// derivatives
	float du = 30.0f * fx * fx * (fx * (fx - 2.0f) + 1.0f);
	float dv = 30.0f * fy * fy * (fy * (fy - 2.0f) + 1.0f);
	float dw = 30.0f * fz * fz * (fz * (fz - 2.0f) + 1.0f);

	// lattice points
	float a = rand(ix+0,iy+0,iz+0);
	float b = rand(ix+1,iy+0,iz+0);
	float c = rand(ix+0,iy+1,iz+0);
	float d = rand(ix+1,iy+1,iz+0);
	float e = rand(ix+0,iy+0,iz+1);
	float f = rand(ix+1,iy+0,iz+1);
	float g = rand(ix+0,iy+1,iz+1);
	float h = rand(ix+1,iy+1,iz+1);

	// bilinear interpolation
	float k0 =  a;
	float k1 =  b - a;
	float k2 =  c - a;
	float k3 =  e - a;
	float k4 =  a - b - c + d;
	float k5 =  a - c - e + g;
	float k6 =  a - b - e + f;
	float k7 = -a + b + c - d + e - f - g + h;

	vec4 n;
	n.x = du * (k1 + k4*v + k6*w + k7*v*w); // dNdx
	n.y = dv * (k2 + k5*w + k4*u + k7*w*u); // dNdy
	n.z = dw * (k3 + k6*u + k5*v + k7*u*v); // dNdz
	n.w = k0 + k1*u + k2*v + k3*w + k4*u*v + k5*v*w + k6*w*u + k7*u*v*w;

	return n;
}
*/

/*
// GLSL version: kept to test performances of the BLAS
// Reference: http://www.gamedev.net/community/forums/topic.asp?topic_id=502913
inline vec3 noise_glsl( const vec2& p ) 
{
	vec2 i = floor(p);
	vec2 f = fract(p);
	vec2 q = f * f * f * (f * (f * 6.0f - 15.0f) + 10.0f);
	vec2 dq = 30.0f * f * f * (f * (f - 2.0f) + 1.0f);

	float a = rand(i.x+0,i.y+0);
	float b = rand(i.x+1,i.y+0);
	float c = rand(i.x+0,i.y+1);
	float d = rand(i.x+1,i.y+1);

	float k0 = a;
	vec2 k12 = vec2(b - a, c - a);
	float k3 = a - b - c + d;

	vec3 n;
	n.xy = dq * (k12 + k3*vec2(q.v,q.u) );
	n.z = k0 + dot(k12,q) + k3*q.u*q.v;
	return n;
}
*/

// TODO gradient (perlin)
// TODO cell     (voronoi)
// TODO metaball (http://www.geisswerks.com/ryan/BLOBS/blobs.html)

/*
	void voronoi(float3 position, out float f1, out float3 pos1, out float f2, out float3 pos2, float jitter=.9, bool manhattanDistance = false )
	{
	  float3 thiscell = floor(position)+.5;
	  f1 = f2 = 1000;
	  float i, j, k;
	  
	  float3 c;
	  for(i = -1; i <= 1; i += 1)
	  {
		for(j = -1; j <= 1; j += 1)
		{
		  for(k = -1; k <= 1; k += 1)
		  {
			float3 testcell = thiscell  + float3(i,j,k);
			float3 randomUVW = testcell * float3(0.037, 0.119, .093);
			float3 cellnoise = perm(perm2d(randomUVW.xy)+randomUVW.z);
			float3 pos = testcell + jitter*(cellnoise-.5);
			float3 offset = pos - position;
			float dist;
			if(manhattanDistance)
			  dist = abs(offset.x)+abs(offset.y) + abs(offset.z); 
			else
			  dist = dot(offset, offset);
			if(dist < f1)
			{
			  f2 = f1;
			  pos2 = pos1;
			  f1 = dist;
			  pos1 = pos;
			}
			else if(dist < f2)
			{
			  f2 = dist;
			  pos2 = pos;
			}
		  }
		}
	  }
	  if(!manhattanDistance)
	  {
		f1 = sqrt(f1);
		f2 = sqrt(f2);
	  }
	}
*/





//
// Domain/amplitude modifiers (Worley, Musgrave)
//
/*
template<vec3 (*noise)(vec2)>
inline static vec3 warp( const vec2& point, const float& scale )
{ 
	vec3 n = noise( scale*point );
	n.xy *= scale;
	return n;
}
*/
template<vec3 (*noise)(vec2)>
inline static vec3 warp( const vec2& point, const mat2& scale )
{ 
	vec3 n = noise( scale*point );
	n.xy *= scale;
	return n;
}
/*
template<vec3 (*noise)(vec2)>
inline static vec3 warp( const vec2& point, const mat3& scale )
{ 
	vec3 n = noise( (scale*vec3(point,1)).xy );
	n.xy = (scale*vec3(n.xy,1)).xy;
	return n;
}
*/

template<vec3 (*noise)(vec2)>
inline static vec3 warp( const vec2& point, const mat2 scale, float ridging ) // sharpness=1
{
	vec3 n = warp<noise>(point, scale); 
	
	n.z -= ridging;
	if( n.z<0 ) n.xy = -n.xy;

	float m = max( ridging, 1-ridging );
	n.z = 1-abs(n.z)/m;
	n.xy = (-1/m) * n.xy;

	return n;
}
/*
template<vec3 (*noise)(vec2)>
inline static vec3 warp( const vec2& point, const mat2& scale, float ridging, float sharpness )
{
	vec3 n = warp<noise>(point,scale); 

	n.z -= ridging;
	if( n.z<0 ) n.xy = -n.xy;

	float m = max( ridging, 1-ridging );
	n.z = 1-abs(n.z)/m;
	n.xy = (-sharpness/m * pow(n.z, sharpness-1)) * n.xy;
	n.z = pow(n.z, sharpness);

	return n;
}
*/







//
// Operations
// 

inline vec3 bias( float offset, const vec3& n )
{ 	
	return vec3( n.x, n.y, n.z + offset );
}

inline vec3 flip( const vec3& n )
{ 
	return vec3(0,0,1) - n;
}

inline vec3 mul( const vec3& n1, const vec3& n2 )
{
	return vec3( n1.xy * n2.z + n2.xy * n1.z, n1.z * n2.z );
}

// TEST
inline vec3 div( const vec3& n1, const vec3& n2 )
{
	return vec3( (n1.xy * n2.z - n2.xy * n1.z)/(n2.z * n2.z) , n1.z / n2.z );
}

inline vec3 sqr( const vec3& n )
{ 
	return mul(n,n);
}

inline vec3 lerp( const vec3& n, const vec3& n1, const vec3& n2 )
{ 
	return mul(flip(n),n1) + mul(n,n2);
}

inline vec3 clamp( const vec3& n, const float z1, const float z2=1.0f )
{ 
	return n.z < z1 ? vec3(0,0,z1) : 
		   n.z > z2 ? vec3(0,0,z2) : 
		   n;
}

/*
template<class T, int N> inline vec<T,N> operator*( const vec<T,N>& v, const mat<T,N>& m )
inline vec<T,N> clamp( const vec<T,N>& n, const T z1, const T z2=T(1) )
{ 
	return n.z < z1 ? vec<T,N>(0,0,z1) : 
		   n.z > z2 ? vec<T,N>(0,0,z2) : 
		   n;
}
*/
inline vec3 crop( const vec3& n, const float z1=0.0f, const float z2=1.0f )
{
	return bias(-z1,clamp(n,z1,z2)) / (z2-z1);
}

// PROP
inline float slope( const vec3& n)
{
	return 1-normalize(vec3(n.xy,1)).z;
}

inline vec3 erode( const vec3& n, const vec3& amount )
{
	//float s = normalize( vec3(n.xy,1) ).z;
	float d = 1/(1+dot(n.xy,n.xy));
	return mul( n, d*amount );
}

inline vec3 modulate( const vec3& n1, const vec3& n2, const float average=0.5f )
{
	float average2 = average*average;
	vec3 n = mul(n1,n2);
	return clamp(n,0,average2)/average + bias(-average2,clamp(n,average2,1))/(1+average);
}



//
// Fractional Brownian Motion
//
/*
// ORIGINAL VERSIONS

// first version
template<vec3 (*noise)(vec2)>
inline vec3 fbm2d( vec2 point, int octaves, float persistence = 0.931f, float lacunarity = 1.717f )
{
	vec3  result    = 0.0f;
	vec3  magnitude = 0.0f;
	float weight    = 1.0f;
	float frequency = 1.0f;
	float scale     = 1.0f;
	
    for( int i=0; i<octaves ; i++ )
    {
		float w = 1/( weight + frequency ); 
		scale *= frequency;

		//vec3 signal = w * noise(point,scale);
		vec2 p = point*scale;		//
		vec3 signal = w * noise(p);	//
		signal.xy *= scale;			//
		
		result     += signal;
		magnitude += w;
        weight    *= persistence;
		frequency *= lacunarity;
    }	

    return result/magnitude;
}


// 2nd version (corrected and warped)
template<vec3 (*noise)(vec2)>
inline vec3 wfbm2d( vec2 point, int octaves, float persistence = 0.931f, mat2 lacunarity = mat2(1.31,0.0,0.0,1.31) )
{
	vec3  result     = 0.0f;
	vec3  magnitude = 0.0f;
	float weight    = 1.0f;
	mat2  frequency = mat2(0.37) + rot( 0.707 );
	mat2  scale     = 1; //mat2(0.73) + rotationMatrix(-0.303 );
	
    for( int i=0; i<octaves ; i++ )
    {
		float w = 1/(weight + frequency[0] ); //+ rand(weight)/10 ) ; 
		scale *= frequency;

		vec3 signal = w * noise(point,scale);
		
		result     += signal;
		magnitude += w;

        weight    *= persistence;
		frequency *= lacunarity;
    }	

    return result/magnitude;
}

// erosive proposal that spoils normals but looks nice
template<vec3 (*noise)(vec2)>
inline vec3 pfbm( vec2 point, int octaves, 
	float weight = 1,  float persistence = 0.5, 
	mat2  scale  = 1,  mat2  lacunarity  = mat2(1.6,-1.2,1.2,1.6),
	float ridging=0.0)
{
	vec3  result     = 0.0f;
	vec3  magnitude = 0.0f;
	
    for( int i=0; i<octaves; i++ )
    {
		vec3 signal = weight * noise(point,scale,ridging);
		
		result     += signal;
		magnitude += weight;

		float k = normalize(vec3(signal.xy,1)).z;
		k = 1 - 0.1*k*k*k - 0.05*k*k;

		weight    *= persistence * k;
		scale     *= lacunarity;
    }	

    return result/magnitude;
}

// a turbulence with 2 noises
template<vec3 (*noise)(vec2)>
inline vec3 pfbm2( vec2 point, int octaves, 
	float weight = 1,  float persistence = 0.5, 
	mat2  scale1 = 1,  mat2  lacunarity1  = mat2(1.6,-1.2,1.2,1.6), 
	mat2  scale2 = 1,  mat2  lacunarity2  = rot(66)*mat2(1.666), 
	float ridging = 0.0 )
{
	vec3  result     = 0.0f;
	vec3  magnitude = 0.0f;
	
    for( int i=0; i<octaves; i++ )
    {
		vec3 n1 = noise(point, scale1, ridging);
		vec3 n2 = noise(point, scale2, ridging);
		vec3 signal = weight * mul(n1,n2);
		
		result     += signal;
		magnitude += weight;

		float k = normalize(vec3(signal.xy,1)).z;
		k = 1 ;//- 0.1*k*k*k - 0.05*k*k;

		weight    *= persistence * k;
		scale1    *= lacunarity1;
		scale2    *= lacunarity2;
    }	

    return result/magnitude;
}

// a turbulence with 3 noises
template<vec3 (*noise)(vec2)>
inline vec3 pfbm3( vec2 point, int octaves, 
	float weight = 1,  float persistence = 0.5, 
	mat2  scale1 = 1,  mat2  lacunarity1  = mat2(1.6,-1.2,1.2,1.6), 
	mat2  scale2 = 1,  mat2  lacunarity2  = rot(66)*mat2(1.666),
	mat2  scale3 = 1,  mat2  lacunarity3  = rot(11)*mat2(1.111),
	float ridging = 0.0)
{
	vec3  result     = 0.0f;
	vec3  magnitude = 0.0f;
	
    for( int i=0; i<octaves; i++ )
    {
		vec3 n1 = noise(point, scale1, ridging);
		vec3 n2 = noise(point, scale2, ridging);
		vec3 n3 = noise(point, scale3, ridging);
		vec3 signal = weight * mul(n1,mul(n2,n3));
		
		result     += signal;
		magnitude += weight;

		float k = normalize(vec3(signal.xy,1)).z;
		k = 1 ;//- 0.1*k*k*k - 0.05*k*k;

		weight    *= persistence * k;
		scale1    *= lacunarity1;
		scale2    *= lacunarity2;
		scale3    *= lacunarity3;
    }	

    return result/magnitude;
}
*/

/*
// TEMP alternative version to be tested
inline static vec3 __fbm( const vec2& point, int octaves, float weight, mat2 scale )
{
	vec3  result     = 0;
	float magnitude = 0;

    while( octaves>0 )
    {
		result     += weight * noise( point, scale );
		magnitude += weight;

		weight *= weight;
		scale *= scale;
		octaves--;
    }	

    return result/magnitude;
}

inline static vec3 fbm( const vec2& point, int octaves, float persistence, mat2 frequency, mat2 lacunarity )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight    *= persistence;
		frequency *= lacunarity;
		scale      = (scale + frequency)*frequency;	

		result     += weight * noise( point, scale );
		magnitude += weight;	
    }	

    return result/magnitude;
}



inline static vec3 ___fbm( const vec2& point, int octaves, float weight, mat2 scale, float ridging )
{

	vec3  result     = 0;
	float magnitude = 0.01;

    for( int i=0; i<octaves; i++ )
    {
		result     += weight * noise( point, scale, ridging );
		magnitude += weight;	

		weight    *= weight;
		scale	  *= scale;
    }	

    return result/magnitude;
}
*/
/*
template<vec3 (*noise)(vec2)>
inline static vec3 fbm( const vec2& point, int octaves, 
		float persistence, mat2 frequency, mat2 lacunarity )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight    *= persistence;
		frequency *= lacunarity;
		scale      = (scale + frequency)*frequency;	

		result     += weight * noise( point, scale );
		magnitude += weight;	
    }	

    return result/magnitude;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm( const vec2& point, int octaves, 
		float persistence, mat2 frequency, mat2 lacunarity, float ridging )
{

	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight    *= persistence;
		frequency *= lacunarity;
		scale      = (scale + frequency)*frequency;	

		result     += weight * noise( point, scale, ridging );
		magnitude += weight;	
    }	

    return result/magnitude;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm2( const vec2& point, int octaves, 
		float persistence, mat2 frequency, mat2 lacunarity)
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight    *= persistence;
		frequency *= lacunarity;
		scale      = (scale + frequency)*frequency;

		vec3 n1 = noise(point, scale    );
		vec3 n2 = noise(point, frequency);

		result     += weight * mul( n1, n2 );
		magnitude += weight;			
    }	

    return result/magnitude;
}


template<vec3 (*noise)(vec2)>
inline static vec3 __fbm2( const vec2& point, int octaves, 
		float weight, mat2 scale1, mat2 scale2, float ridging)
{
	vec3  result     = 0;
	float magnitude = 0;

    for( int i=0; i<octaves; i++ )
    {
		vec3 n1 = noise(point, scale1, ridging);
		vec3 n2 = noise(point, scale2, ridging);

		result     += weight * mul( n1, n2 );
		magnitude += weight;	

		weight *= weight;
		scale1  *= scale2;
		scale2  *= scale1;
    }	

    return result/magnitude;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm2( vec2 point, int octaves, 
		float persistence, mat2 frequency, mat2 lacunarity, float ridging)
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight    *= persistence;
		frequency *= lacunarity;
		scale      = (scale + frequency)*frequency;

		vec3 n1 = noise(point, scale,     ridging);
		vec3 n2 = noise(point, frequency, ridging);

		result     += weight * mul( n1, n2 );
		magnitude += weight;	
    }	

    return result/magnitude;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm3( vec2 point, int octaves, 
		float persistence, mat2 scale1, mat2 scale2, mat2 scale3, float ridging)
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight *= persistence;
		scale1	= (scale1 + scale2)*scale2;
		scale2	= (scale2 + scale3)*scale3;
		scale3	= (scale3 + scale1)*scale1;

		vec3 n1 = warp<noise>(point, scale1, ridging);
		vec3 n2 = warp<noise>(point, scale2, ridging);
		vec3 n3 = warp<noise>(point, scale3, ridging);

		result     += weight * mul( mul( n1, n2 ), n3);
		magnitude += weight;	
    }	

    return result/magnitude;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm3( vec2 point, int octaves, 
		float persistence, mat2 scale1, mat2 scale2, mat2 scale3)
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight *= persistence;
		scale1	= (scale1 + scale2)*scale2;
		scale2	= (scale2 + scale3)*scale3;
		scale3	= (scale3 + scale1)*scale1;

		vec3 n1 = warp<noise>(point, scale1);
		vec3 n2 = warp<noise>(point, scale2);
		vec3 n3 = warp<noise>(point, scale3);

		result    += weight * mul( mul( n1, n2 ), n3);
		magnitude += weight;	
    }	

    return result/magnitude;
}


// TODO versions with sharpness != 1
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2( vec2 point, int octaves, 
		float persistence, mat2 frequency, mat2 lacunarity, float ridging, float sharpness)
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		weight    *= persistence;
		frequency *= lacunarity;
		scale      = (scale + frequency)*frequency;

		vec3 n1 = noise(point, scale,     ridging, sharpness);
		vec3 n2 = noise(point, frequency, ridging, sharpness);

		result     += weight * mul( n1, n2 );
		magnitude += weight;	
    }	

    return result/magnitude;
}


	template<vec3 (*dnoise)(vec2)>
	inline static vec3 fPm_optimized( vec2 point, int octaves, 
			float persistence = 0.212179,
			mat2 frequency = mat2(0.37) + rot(0.707),
			mat2 lacunarity = mat2(0.07) + mat2(1.37,0.0,0.0,1.31) )
	{
		vec3  value     = 0;
		vec3  magnitude = 0;
		float weight	= 1;
		mat2  scale		= mat2(1);

	    for( int i=0; i<octaves; i++ )
	    {
			frequency *= lacunarity;
			scale      = (scale + frequency)*frequency;	

			vec3 signal = weight * dnoise( scale*point );
			signal.xy *= scale;

			value     += signal;
			magnitude += weight;	

			//signal.z = 1;
			//float k = normalize(signal).z;
			weight    *= persistence;// * (1-0.987*k*k);
			//frequency = frequency*(1-0.00073f*k*k);	
	    }	

	    return value/magnitude;
	}


	template<vec3 (*dnoise)(vec2)>
	inline static vec3 fbm1_scale_after( vec2 point, int octaves, 
			float persistence = 0.212179,
			mat2 frequency = mat2(0.37) + rot(0.707),
			mat2 lacunarity = mat2(0.07) + mat2(1.37,0.0,0.0,1.31) )
	{
		mat2  scale		= mat2(1);
		float weight	= 1;
		vec3  magnitude = 0;
		vec3  value     = 0;

	    for( int i=0;; i++ )
	    {
			magnitude += weight;
			value     += weight * noise( point, scale );
			if( i>=octaves ) break;
			weight    *= persistence;
			frequency *= lacunarity;
			scale      = (scale + frequency)*frequency;	
	    }	

	    return value/magnitude;
	}
*/
/*
	inline static vec3 fpm19_ridged( vec2 point, int octaves, 
			float persistence = 0.212179,
			mat2 frequency = mat2(0.37) + rot(0.707),
			mat2 lacunarity = mat2(0.07) + mat2(1.37,0.0,0.0,1.31), float ridging )
	{
		mat2  scale		= mat2(1);
		float weight	= 1;
		vec3  magnitude = 0;
		vec3  value     = 0;

	    for( int i=0; i<octaves; i++ )
	    {
			weight    *= persistence;
			frequency *= lacunarity;
			scale      = (scale + frequency)*frequency;	
			magnitude += weight;
			value     += weight * noise( point, scale, ridging );
	    }	

	    return value/magnitude;
	}

	*/




// ------------------------------
//
// RESEARCH section
//
// ------------------------------



// NOTE: First fBm presented by Musgrave in "Procedural Fractal Terrains", 
// named as "archetypal fBm" that is statistically homogeneous and isotropic,
// where: 
//   H ("fractal increment") has to do with smoothness (H->0 => rougher, H->1 => smoother)
//   Lacunarity is the gap between frequencies and has to do with bandwidth (L->2.0)
//   Octaves has to do with LOD and should be Log2(resolution)-2 

template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_archetypal( vec2 point, float H, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	float wi		= pow(lacunarity.a,-H); // TODO lacunarity.a -> pow(det(lacunarity),1/2)
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		result += weight * noise(point, scale);	// TODO noise2
		magnitude += weight;
		scale *= lacunarity;
		weight *= wi;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * noise(point,scale);
	magnitude += weight;

    return result/magnitude;
}

// [MUSGRAVE] Scaling higher frequencies by the previous frequency amplitude
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_hybrid1( vec2 point, float H, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	vec3  previous  = vec3(0,0,1);
	float wi		= pow(lacunarity.a,-H); // TODO lacunarity.a -> pow(det(lacunarity),1/2)
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		vec3 signal = noise(point, scale);
		result += weight * mul(previous, signal);
		magnitude += weight;// * previous.z;
		scale *= lacunarity;
		weight *= wi;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, noise(point,scale) );
	magnitude += weight;// * previous.z;

    return result/magnitude;
}
/*
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_hybrid( vec2 point, float harshness, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	vec3  previous  = noise(point);
	float wi		= pow(lacunarity.a,-harshness); // TODO lacunarity.a -> pow(det(lacunarity),-Harshness/N)
	mat2  scale		= mat2(lacunarity);

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point, scale);
		result += weight * mul(previous, signal);
		magnitude += weight;// * previous.z;
		scale *= lacunarity;
		weight *= wi;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale) );
	magnitude += weight;// * previous.z;

    return result/magnitude;
}

// PROP Failed
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_hybrid_eroded( vec2 point, float H, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	vec3  previous  = noise(point);
	float wi		= pow(lacunarity.a,-H);
	mat2  scale		= mat2(lacunarity);

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = noise(point, scale);
//float a = atan2(signal.y,signal.x); a = (sin(a*a*5)+1)/2;

		result += weight * mul(previous, signal);
		magnitude += weight;// * previous.z;
		scale *= lacunarity;

		float a = atan2(signal.y,signal.x);
		a = (sin(a*a*5)+1)/2;
		float e = (1-a*normalize(vec3(signal.xy,1)).z);
		
		weight *= wi * e;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, noise(point,scale) );
	magnitude += weight;// * previous.z;

    return result/magnitude;
}



// [DENICOLA] Scaling higher frequencies by the previous frequency slope
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_archetypal_eroded( vec2 point, float H, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	vec3  previous  = vec3(0,0,1);
	float wi		= pow(lacunarity.a,-H);
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		vec3 signal = noise(point, scale);
		result += weight * normalize(vec3(previous.xy,1)).z * signal;
		magnitude += weight; // * normalize(vec3(previous.xy,1)).z;
		scale *= lacunarity;
		weight *= wi;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * noise(point,scale);
	magnitude += weight;

    return result/magnitude;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_archetypal_eroded_2( vec2 point, float H, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	vec3  previous  = vec3(0,0,1);
	float wi		= pow(lacunarity.a,-H);
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		vec3 signal = noise(point, scale);
		result += weight * signal;
		magnitude += weight;
		scale *= lacunarity;
		weight *= wi * normalize(vec3(signal.xy,1)).z;;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * noise(point,scale);
	magnitude += weight;

    return result/magnitude;
}

// [ALTAVILLA]
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_archetypal_enriched( vec2 point, float H, mat2 lacunarity, float octaves )
{
	vec3  result     = 0;
	float magnitude = 0;
	float weight	= 1;
	float wi		= pow(lacunarity.a,-H); // TODO lacunarity.a -> pow(det(lacunarity),1/2)
	mat2  scale		= mat2(1);

	//for( int i=0; i<octaves; i++ )
	float i=0;
	while(i<octaves)
    {
		vec3 signal = noise(point, scale);
		result += weight * signal;
		magnitude += weight;
		scale *= lacunarity;
		weight *= wi;
		octaves += normalize(vec3(signal.xy,1)).z*0.8;
		i+=1;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * noise(point,scale);
	magnitude += weight;

    return result/magnitude;
}
*/


//
//
// APPROVING
//
//
/*
// archetypal
template<vec3 (*noise)(vec2)>
inline static vec3 fbm1( vec2 point, float octaves, mat2 lacunarity, float roughness )
{
	vec3  result    = 0;
	float magnitude = 0;
	float weight	= 1;
	float weighti	= pow(det(lacunarity),-1/(2*roughness));
	mat2  scale		= mat2(1);

    for( int i=0; i<octaves; i++ )
    {
		result += weight * warp<noise>(point, scale);
		magnitude += weight;
		scale *= lacunarity;
		weight *= weighti;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * warp<noise>(point, scale);
	magnitude += weight;

    return result/magnitude;
}

// hybrid
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2( vec2 point, float octaves, mat2 scale, mat2 lacunarity, float roughness )
{
	vec3  result    = 0;
	float magnitude = 0;
	float weight	= 1;
	float weighti	= pow(det(lacunarity),-1/(2*roughness));
	vec3  previous  = warp<noise>(point,scale);
	scale *= lacunarity;

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale);
		result += weight * mul(previous, signal);
		magnitude += weight;
		scale *= lacunarity;
		weight *= weighti;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale) );
	magnitude += weight;

    return result/magnitude;
}

// working bad
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_eroded( vec2 point, float octaves, mat2 scale, mat2 lacunarity, float roughness, float erosion )
{
	vec3  result    = 0;
	float magnitude = 0;
	float weight	= 1;
	float weighti	= pow(det(lacunarity),-1/(2*roughness));
	vec3  previous  = warp<noise>(point,scale);
	scale *= lacunarity;

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale);
		result += weight * mul(previous, signal);
		magnitude += weight;
		scale *= lacunarity;
			float n = 1/(1+dot(result.xy,result.xy));
			float h = result.z;
			float k = 1 - (1-0.0f*h*h)*n*n*n*n*erosion;
		weight *= weighti * k;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale) );
	magnitude += weight;

    return result/magnitude;
}


//
// 2011-04-07 new proposal
//
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_eroded_ridged_20110407( vec2 point, float octaves, 
		mat2 scale, mat2 lacunarity, float roughness, float erosion, float ridging, vec2 slope )
{
	vec3  result		= 0;
	float weight		= 1;
	float persistence	= pow(det(lacunarity),-1/(2*roughness));

	// octave 0
	vec3  previous		= warp<noise>(point,scale,ridging);
	scale *= lacunarity;

	// octaves 1..n
    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale,ridging);
		result += weight * mul(previous, signal);
		scale *= lacunarity;
			float d = dot(result.xy,result.xy);
			float n = 1/(1+d);
			float k = 1-0.0435287*n;//*n*n*n*n*n*n*n;
			//float k = 1-0.15287*smoothstep(0.0f,2.0f,k*k*k*k);
			weight    *= persistence * k;
		previous = signal;
    }	

	// octave rest
	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale,ridging) );

    return result;
}
*/

// TODO return a 4th component as flatness of the first octave

template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_eroded_s2( vec2 point, float octaves, mat2 scale, mat2 lacunarity, float roughness, float erosion, float scraping )
{
	vec3  result		= 0;
	float weight		= 1;
	float persistence	= pow(det(lacunarity),-1/roughness);

	// octave 0
	vec3  previous		= warp<noise>(point,scale);
	scale *= lacunarity;

	// octave 1..n
    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale);
		scale *= lacunarity;

		result += weight * mul(previous, signal);
		previous = signal;

			float d = dot(result.xy,result.xy);
			float N = 1/(1 + d);
			float e = N*(0.99f+erosion*0.01f) - scraping*(signal.x-signal.y);
			float k = 1 - smoothstep(0.99f, 1.00f, e );

		weight  *= persistence * k;
    }	

	// remaining octave
	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale) );

    return result;
}

// attempt to correct normals
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_eroded3( vec2 point, float octaves, mat2 scale, mat2 lacunarity, float roughness, float erosion, float scraping )
{
	vec3  result		= 0;
	float weight		= 1;
	float persistence	= pow(det(lacunarity),-1/roughness);

	// octave 0
	vec3  previous		= warp<noise>(point,scale);
	scale *= lacunarity;

	// octave 1..n
    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale);
		scale *= lacunarity;

		result += weight * mul(previous, signal);
		previous = signal;

			float d = dot(result.xy,result.xy);
			float N = 1/(1+d*d);
			//float e = N*(0.99f+erosion*0.01f) ;//- scraping*(signal.x-signal.y);
			//float k = 1 - smoothstep(0.99f, 1.00f, e );

		weight  *= persistence * //k;
					 //smoothstep(0.01f, 0.02f, N ) * 0.919;
					 N * 0.937;
    }	

	// remaining octave
	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale) );

    return result;
}



/*
// without magnitude
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_eroded_s( vec2 point, float octaves, mat2 scale, mat2 lacunarity, float roughness, float erosion, vec2 slope )
{
	vec3  result    = 0;
	float weight	= 1;
	float weighti	= pow(det(lacunarity),-1/(2*roughness));
	vec3  previous  = warp<noise>(point,scale);
	scale *= lacunarity;

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale);
		result += weight * mul(previous, signal);
		scale *= lacunarity;
			float d = dot(result.xy,result.xy);
			float n = 1/(1+d*d); // as approximation of normal length 
			//float h = result.z;
			//float k = 1 - ((1-0.0f*h*h)*n)*erosion;
				float k = 1-n*n*n*erosion;
				k = k*k*k;
		weight *= weighti * k;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale) );

    return result;
}

template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_eroded_ridged( vec2 point, float octaves, mat2 scale, mat2 lacunarity, float roughness, float erosion, float ridging )
{
	vec3  result    = 0;
	float magnitude = 0;
	float weight	= 1;
	float weighti	= pow(det(lacunarity),-1/(2*roughness));
	vec3  previous  = warp<noise>(point,scale,ridging);
	scale *= lacunarity;

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale);
		result += weight * mul(previous, signal);
		magnitude += weight;
		scale *= lacunarity;
			float n = 1/(1+dot(result.xy,result.xy));
			float h = result.z;
			float k = 1 - (1-0.0f*h*h)*n*n*n*erosion;
		weight *= weighti * k;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale,ridging) );
	magnitude += weight;

    return result/magnitude;
}


// heavier and not versatile 
template<vec3 (*noise)(vec2)>
inline static vec3 fbm2_ridged( vec2 point, float octaves, mat2 lacunarity, float roughness, float ridging )
{
	vec3  result    = 0;
	float magnitude = 0;
	float weight	= 1;
	float weighti	= pow(det(lacunarity),-1/(2*roughness));
	mat2  scale		= mat2(lacunarity);
	vec3  previous  = noise(point);

    for( int i=1; i<octaves; i++ )
    {
		vec3 signal = warp<noise>(point,scale,ridging);
		result += weight * mul(previous, signal);
		magnitude += weight;
		scale *= lacunarity;
		weight *= weighti;
		previous = signal;
    }	

	weight *= (octaves - (int)octaves);
	result += weight * mul(previous, warp<noise>(point,scale,ridging) );
	magnitude += weight;

    return result/magnitude;
}
*/


// http://lists.apple.com/archives/perfoptimization-dev/2005/Jan/msg00051.html

float satan2( float x, float y )
{
	const float ONEQTR_PI = PI / 4.0;
	const float THRQTR_PI = 3.0 * PI / 4.0;
	float r, angle;
	float abs_y = fabs(y) + 1e-10f;      // kludge to prevent 0/0 condition
	if ( x < 0.0f )
	{
		r = (x + abs_y) / (abs_y - x);
		angle = THRQTR_PI;
	}
	else
	{
		r = (x - abs_y) / (x + abs_y);
		angle = ONEQTR_PI;
	}
	angle += (0.1963f * r * r - 0.9817f) * r;
	//if ( y < 0.0f ) angle = -angle;
	return angle;
}

float trimod( float x, float y )
{
	x = mod(x,2*y);
	if( x>y ) x = y-x;
	return x;
}