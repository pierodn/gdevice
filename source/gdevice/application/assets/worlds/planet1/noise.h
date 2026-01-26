
#include "gpu/program.h"

GLSL_(

// https://www.shadertoy.com/view/WlfSzX    
vec4 noise(vec2 point) 
{		
    vec2 i = floor(point);
    vec2 f = fract(point);

	vec2 u = f*f*(3.0-2.0*f);  //f*f*f*(f*(6.0*f-15.0)+10.0);
	vec2 du = 6.0*f*(1.0-f);   //30.0*f*f*(f*(f-2.0)+1.0);
	
	const float WIDTH = 0.3137; // 13.0
    vec4 L = vec4(0.0, 1.0, WIDTH, WIDTH + 1.0);	

    L = fract((L + dot(i, L.yz))*0.1731);
    L += L*(L + 1.0);
    L = fract(7.7771*L*L) - 0.5;

	L = vec4(L.x, L.y-L.x, L.z-L.x, L.x-L.y-L.z+L.w);
	return vec4(du*(L.yz + L.w*u.yx), L.x + L.y*u.x + L.z*u.y + L.w*u.x*u.y, 0.0 );
}

// Noise algebra
const vec4 unit = vec4(0.0, 0.0, 1.0, 0.0);
const mat2 identity = mat2(1.0,0.0,0.0,1.0);
vec4 noise(vec2 point, mat2 scale)  { vec4 n = noise(scale*point); n.xy *= scale; return n; }
vec4 noise(vec2 point, float scale) { vec4 n = noise(scale*point); n.xy *= scale; return n; }
vec4 bias(vec4 n, float offset )    { return n + offset * unit; }
vec4 multiply(vec4 n1, vec4 n2 )    { return vec4(n1.xy * n2.z + n2.xy * n1.z, n1.z * n2.z, n1.z * n2.w + n2.z*n1.w); }
//vec4 quotient(vec4 n1, vec4 n2)		{ return vec4((n1.xy * n2.z - n2.xy * n1.z)/(n2.z*n2.z),  n1.z/n2.z, 0.0); }	
vec4 lerp(vec4 a, vec4 n1, vec4 n2 ){ return multiply(vec4(0.0, 0.0, 1.0, 0.0)-a,n1) + multiply(a,n2); }
//vec4 power(vec4 n, float k)         { return vec4(k*pow(n.z, k-1.0)*n.xy, pow(n.z, k), k*pow(n.w, k+0.1)); }
vec4 power(vec4 n, float k)         { return vec4(k*pow(n.z, k-1.0)*n.xy, pow(n.z, k), n.w); }
vec4 maximum(vec4 n1, vec4 n2)      { return n1.z > n2.z ? n1 : n2; }
vec4 minimum(vec4 n1, vec4 n2)      { return n1.z < n2.z ? n1 : n2; }
vec4 saturate(vec4 n, float z1, float z2)   { return n.z < z1 ? z1*unit : n.z > z2 ? z2*unit : n; }
vec4 saturate(vec4 n)               { return saturate(n, 0.0, 1.0); }
//float converge(float x0, float x1, float x) { return 1.0/(x + 1.0/(x0-x1)) + x1; } 
vec4 saturate(vec4 n, vec4 a, vec4 b) { return n.z < a.z ? a : n.z >= b.z ? vec4(-b.xy, b.zw) : n; }

vec4 smoothStep(float h0, float h1, vec4 n) { n = vec4(n.xy, n.z - h0, 0.0)/(h1 - h0); n = saturate(n, 0.0, 1.0); return multiply(multiply(n,n), bias(-2.0 * n, +3.0)); }
vec4 absolute(vec4 n, float zero)   { return n.z >= zero ? n : -n; }    
vec4 invert(vec4 n)                 { float z = 1.0/(n.z + 1.0); return vec4(-z*z*n.xy, z, 0.0); }
vec4 invert2(vec4 n)			    { return vec4(-n.xy, 1.0 - n.z, 0.0); }
vec4 minus(float t, vec4 n)         { return vec4(-n.xy, t - n.z, 0.0); }
vec4 tone(vec4 n, float t)			{ return smoothStep(0.0, 1.0, multiply(n, invert(bias(n,t))))*(1.0 + t); }
vec4 harmonic(vec2 p)               { return bias(0.5*vec4(cos(p.x)*sin(p.y), sin(p.x)*cos(p.y), sin(p.x)*sin(p.y), 0.0), 0.5); }

vec4 exponent2(vec4 a)
{
   // TODO
   return vec4( exp2(a.z), 1.0, 1.0, 0.0 );
}

vec4 harmonicVoronoi(vec2 point, float scale)
{
    point *= scale;
    
    const mat2 R1 = mat2( -0.7373, -0.6755, +0.6755, -0.7373 ); // 137.5077 (Golden Angle)
    //const mat2 R2 = mat2( +0.3623, -0.9320, +0.9320, +0.3623 ); // 19.64390 (Ga/7)
    //const mat2 R3 = mat2( -0.7159, -0.6981, +0.6981, -0.7159 ); // 27.50154 (Ga/5)
    vec4 F1 = harmonic(point); 
    vec4 F2 = harmonic(R1*point + R1[0]); F2.xy *= R1;
    //vec4 F3 = harmonic(R2*point + R2[0]); F3.xy *= R2;
    //vec4 F4 = harmonic(R3*point + R3[0]); F4.xy *= R3;

    const float k = 2.0;
    F1 = power(F1, -k) + power(F2, -k);
    F1 = power(F1, -1.0/k);

    vec4 result = vec4(F1.xyz, (F2-F1).z);

	result.xy *= scale;
    return result;
}

vec4 fbm(vec2 p, int octaves, 
    float amplitude, // = 1.0,
    float gain, // = 0.5,
    mat2 frequency, // = mat2(1.0,0.0,0.0,1.0),
    mat2 lacunarity // = 2.03*mat2(0.8,-0.6,0.6,0.8) 
    )
{
    vec4 signal = vec4(0.0);
    for( int i=0; i<octaves; i++ ) {
        signal += amplitude*noise(p, frequency);
		amplitude *= gain;
		frequency *= lacunarity;
    }
    return signal;
}

vec4 fbm(vec2 p, int octaves)
{
	return fbm(p, octaves, 1.0, 0.5, mat2(1.0,0.0,0.0,1.0), 2.03*mat2(0.8,-0.6,0.6,0.8));
}

vec4 fbm(vec2 p, float scale)
{
    const mat2 M2 = mat2(0.8,-0.6,0.6,0.8);

    vec4 f = vec4(0.0);
    f += 0.5000*noise(p*SCALE); p = M2*p*2.01;
    f += 0.2500*noise(p*SCALE); p = M2*p*2.02;
    f += 0.1250*noise(p*SCALE); p = M2*p*2.03;
    f += 0.0625*noise(p*SCALE);
    return f/0.9375;
}
    
vec4 hybrid(vec2 point, int octaves, 
    inout vec4 signal, // = vec4(0.0), 
    inout vec4 weight, // = vec4(0.0, 0.0, 0.17, 0.0),
    inout float scale,
    inout float frequency)
{  
    float H = 0.11;             // 0.13 - The lower, the rougher
    float offset = 0.69;	    // The heigher, the rougher (**)
    float lacunarity = 2.03;    // 2.00 - Decomposes bodies but increases spikes (?)
    float maxFrequency = 17.0;  // The lower, the rougher
    float maxWeight = 0.05;     // 1.00 - Low values reduce the spikes, reduce height, introduce artifacts
 				
    for(int i=0; i<octaves; i++)
    {
        vec4 noise = noise(point, scale);
        noise = bias(noise, offset) * pow(frequency, -H); 
        signal += weight; 
        weight = multiply(weight, noise);
        weight = minimum(weight, maxWeight*unit); 
        weight.w = dot(weight.xy, weight.xy);
        frequency = min(frequency*lacunarity, maxFrequency);
	    scale *= lacunarity;
    }
    return signal;
}
            
)