#pragma once

#include <stdio.h>
#include <stdarg.h>

#include "type/glsl.h"

#if defined(_MSC_VER) 
	#pragma warning(disable: 4172) // Returning local variable address
    #pragma warning(disable: 4996) // Unsafe functions like sprintf
#endif

// TODO Make multi-threaded
char* format(const char* str, ...) 
{
    const int bufferSize = 1024;
    static char buffer[bufferSize];

    va_list args;
    va_start(args, str);
    vsnprintf(buffer, bufferSize-1, str, args);
    va_end(args);

    return buffer;
}

char* ToString(double a){return format("%.8g", a);}

template<class T, int N> inline
char* ToString(const vec<T,N>& v)
{
    const int bufferSize = 1024;
    static char buffer[bufferSize];

    char* p = buffer;
    p += sprintf(p, "vec%d(", N );
    for(int i=0; i<N; i++) p += sprintf(p, "%.8g, ", v[i]);
    sprintf(p-2, ")", 0);

    return buffer;
}

template<class T, int N, int M> inline
char* ToString(const mat<T,N,M>& m)
{
    const int bufferSize = 4096;
    static char buffer[bufferSize];

    char* p = buffer;
    p += sprintf(p, "mat%d", N );
    if(N != M) p += sprintf(p, "x%d", M);
    p += sprintf(p, "%s", "(");
    for(int i=0; i<N*M; i++) p += sprintf(p, "%.8g, ", m[i]);
    sprintf(p-2, ")", 0);

    return buffer;
}

#define PRINTLN(x) printf("%s = %s\n", #x, ToString(x));

#define ASSERT_EQUAL(x,y) \
    if(distance((x),(y)) > EPSILON) { \
    printf("\nAssert failed at line %d\n\nExpression: " #x "\n     Value: %s\n  Expected: " #y "\n\n", __LINE__, ToString(x)); }





void test_all()
{
    printf("Executing %s in %s\n\n", __FUNCTION__, __FILE__);

    //
    // Construction
    //
    const float  pi         = float(PI);
    const double EULER		= 2.718281828459045235360;
    const float  euler      = float(EULER);

    float s = 0.5;

    vec2 u = vec2(1.0, 2.0);
    vec2 v = vec2(3.0, 4.0);
    vec3 a = vec3(1.0, 2.0, 3.0);
    vec3 b = vec3(4.0, 5.0, 6.0);
    vec4 p = vec4(1.0, 2.0, 3.0, 4.0);
    vec4 q = vec4(5.0, 6.0, 7.0, 8.0);
    mat2 M = mat2(-3, 1, 5, -2);
    mat2 N = inverse(M);
    mat3 A = mat3(0, -3, -2, 1, -4, -2, -3, 4, 1);
    mat3 B = inverse(A);
    mat4 T = mat4(-2, -2, 3, 0, 0, 1, -1, 0, 0.5, 0, -0.5, 0, 0, 0, 0, 1);
    mat4 Q = inverse(T);

    //
    // Component access
    //
    ASSERT_EQUAL(u.x, 1.0f);
    ASSERT_EQUAL(u.y, 2.0f);
    ASSERT_EQUAL(a.x, 1.0f);
    ASSERT_EQUAL(a.y, 2.0f);
    ASSERT_EQUAL(a.z, 3.0f);
    ASSERT_EQUAL(a.xy, vec2(1.0, 2.0));
    ASSERT_EQUAL(a.yz, vec2(2.0, 3.0));
    ASSERT_EQUAL(p.x, 1.0f);
    ASSERT_EQUAL(p.y, 2.0f);
    ASSERT_EQUAL(p.z, 3.0f);
    ASSERT_EQUAL(p.w, 4.0f);
    ASSERT_EQUAL(p.xy, vec2(1.0, 2.0));
    ASSERT_EQUAL(p.yz, vec2(2.0, 3.0));
    ASSERT_EQUAL(p.zw, vec2(3.0, 4.0));
    ASSERT_EQUAL(p.xyz, vec3(1.0, 2.0, 3.0));
    ASSERT_EQUAL(p.yzw, vec3(2.0, 3.0, 4.0));

    //
    // Component-wise operations
    //
    ASSERT_EQUAL(1.0 + 1.0, 2.0);

    ASSERT_EQUAL(   -v, vec2(-3.0, -4.0));
    ASSERT_EQUAL(u + v, vec2(4.0, 6.0));
    ASSERT_EQUAL(u - v, vec2(-2.0));
    ASSERT_EQUAL(u * v, vec2(3.0, 8.0));
    ASSERT_EQUAL(u / v, vec2(0.333333, 0.5));
    ASSERT_EQUAL(u + s, vec2(1.5, 2.5));
    ASSERT_EQUAL(u - s, vec2(0.5, 1.5));
    ASSERT_EQUAL(u * s, vec2(0.5, 1.0));
    ASSERT_EQUAL(u / s, vec2(2.0, 4.0));
    ASSERT_EQUAL(s + v, vec2(3.5, 4.5));
    ASSERT_EQUAL(s - v, vec2(-2.5, -3.5));
    ASSERT_EQUAL(s * v, vec2(1.5, 2.0));
    ASSERT_EQUAL(s / v, vec2(0.1666666, 0.125));

    ASSERT_EQUAL(   -b, vec3(-4.0, -5.0, -6.0));
    ASSERT_EQUAL(a + b, vec3(5.0, 7.0, 9.0));
    ASSERT_EQUAL(a - b, vec3(-3.0));
    ASSERT_EQUAL(a * b, vec3(4.0, 10.0, 18.0));
    ASSERT_EQUAL(a / b, vec3(0.25, 0.4, 0.5));
    ASSERT_EQUAL(a + s, vec3(1.5, 2.5, 3.5));
    ASSERT_EQUAL(a - s, vec3(0.5, 1.5, 2.5));
    ASSERT_EQUAL(a * s, vec3(0.5, 1.0, 1.5));
    ASSERT_EQUAL(a / s, vec3(2.0, 4.0, 6.0));
    ASSERT_EQUAL(s + b, vec3(4.5, 5.5, 6.5));
    ASSERT_EQUAL(s - b, vec3(-3.5, -4.5, -5.5));
    ASSERT_EQUAL(s * b, vec3(2.0, 2.5, 3.0));
    ASSERT_EQUAL(s / b, vec3(0.125, 0.1, 0.0833333));

    ASSERT_EQUAL(   -p, vec4(-1.0, -2.0, -3.0, -4.0));
    ASSERT_EQUAL(p + q, vec4(6.0, 8.0, 10.0, 12.0));
    ASSERT_EQUAL(p - q, vec4(-4.0));
    ASSERT_EQUAL(p * q, vec4(5.0, 12.0, 21.0, 32.0));
    ASSERT_EQUAL(p / q, vec4(0.2, 0.333333, 0.428571, 0.5));
    ASSERT_EQUAL(p + s, vec4(1.5, 2.5, 3.5, 4.5));
    ASSERT_EQUAL(p - s, vec4(0.5, 1.5, 2.5, 3.5));
    ASSERT_EQUAL(p * s, vec4(0.5, 1.0, 1.5, 2.0));
    ASSERT_EQUAL(p / s, vec4(2.0, 4.0, 6.0, 8.0));
    ASSERT_EQUAL(s + q, vec4(5.5, 6.5, 7.5, 8.5));
    ASSERT_EQUAL(s - q, vec4(-4.5, -5.5, -6.5, -7.5));
    ASSERT_EQUAL(s * q, vec4(2.5, 3.0, 3.5, 4.0));
    ASSERT_EQUAL(s / q, vec4(0.1, 0.0833333, 0.071428, 0.0625));

    ASSERT_EQUAL(   -N, mat2(2, 1, 5, 3));
    ASSERT_EQUAL(M + N, mat2(-5, 0, 0, -5));
    ASSERT_EQUAL(M - N, mat2(-1, 2, 10, 1));
    ASSERT_EQUAL(matrixCompMult(M, N), mat2(6, -1, -25, 6)); // M * N is a geometric function
    ASSERT_EQUAL(M / N, mat2(1.5, -1, -1, 0.666666));
    ASSERT_EQUAL(M + s, mat2(-2.5, 1.5, 5.5, -1.5));
    ASSERT_EQUAL(M - s, mat2(-3.5, 0.5, 4.5, -2.5));
    ASSERT_EQUAL(M * s, mat2(-1.5, 0.5, 2.5, -1));
    ASSERT_EQUAL(M / s, mat2(-6, 2, 10, -4));
    ASSERT_EQUAL(s + N, mat2(-1.5, -0.5, -4.5, -2.5));
    ASSERT_EQUAL(s - N, mat2(2.5, 1.5, 5.5, 3.5));
    ASSERT_EQUAL(s * N, mat2(-1, -0.5, -2.5, -1.5));
    ASSERT_EQUAL(s / N, mat2(-0.25, -0.5, -0.1, -0.166666));

    ASSERT_EQUAL(   -B, mat3(-4, 5, 2, -5, 6, 2, 8, -9, -3));
    ASSERT_EQUAL(A + B, mat3(4, -8, -4, 6, -10, -4, -11, 13, 4));
    ASSERT_EQUAL(A - B, mat3(-4, 2, 0, -4, 2, 0, 5, -5, -2));
    ASSERT_EQUAL(matrixCompMult(A, B), mat3(0, 15, 4, 5, 24, 4, 24, 36, 3));
    ASSERT_EQUAL(A / B, mat3(0, 0.6, 1, 0.2, 0.666666, 1, 0.375, 0.444444, 0.333333));
    ASSERT_EQUAL(A + s, mat3(0.5, -2.5, -1.5, 1.5, -3.5, -1.5, -2.5, 4.5, 1.5));
    ASSERT_EQUAL(A - s, mat3(-0.5, -3.5, -2.5, 0.5, -4.5, -2.5, -3.5, 3.5, 0.5));
    ASSERT_EQUAL(A * s, mat3(0, -1.5, -1, 0.5, -2, -1, -1.5, 2, 0.5));
    ASSERT_EQUAL(A / s, mat3(0, -6, -4, 2, -8, -4, -6, 8, 2));
    ASSERT_EQUAL(s + B, mat3(4.5, -4.5, -1.5, 5.5, -5.5, -1.5, -7.5, 9.5, 3.5));
    ASSERT_EQUAL(s - B, mat3(-3.5, 5.5, 2.5, -4.5, 6.5, 2.5, 8.5, -8.5, -2.5));
    ASSERT_EQUAL(s * B, mat3(2, -2.5, -1, 2.5, -3, -1, -4, 4.5, 1.5));
    ASSERT_EQUAL(s / B, mat3(0.125, -0.1, -0.25, 0.1, -0.0833333, -0.25, -0.0625, 0.0555555, 0.166666));

    ASSERT_EQUAL(   -Q, mat4(1, 2, 2, 0, 1, 1, 4, 0, 1, 2, 4, 0, 0, 0, 0, -1)); 
    ASSERT_EQUAL(T + Q, mat4(-3, -4, 1, 0, -1, 0, -5, 0, -0.5, -2, -4.5, 0, 0, 0, 0, 2));
    ASSERT_EQUAL(T - Q, mat4(-1, 0, 5, 0, 1, 2, 3, 0, 1.5, 2, 3.5, 0, 0, 0, 0, 0));
    ASSERT_EQUAL(matrixCompMult(T, Q), mat4(2, 4, -6, 0, -0, -1, 4, 0, -0.5, -0, 2, 0, 0, 0, 0, 1));
    ASSERT_EQUAL(T / (Q + s), mat4(4, 4/3.0, -2, 0, -0, -2, 0.285714, 0, -1, -0, 0.142857, 0, 0, 0, 0, 0.666666));
    ASSERT_EQUAL(T + s, mat4(-1.5, -1.5, 3.5, 0.5, 0.5, 1.5, -0.5, 0.5, 1, 0.5, 0, 0.5, 0.5, 0.5, 0.5, 1.5));
    ASSERT_EQUAL(T - s, mat4(-2.5, -2.5, 2.5, -0.5, -0.5, 0.5, -1.5, -0.5, 0, -0.5, -1, -0.5, -0.5, -0.5, -0.5, 0.5));
    ASSERT_EQUAL(T * s, mat4(-1, -1, 1.5, 0, 0, 0.5, -0.5, 0, 0.25, 0, -0.25, 0, 0, 0, 0, 0.5));
    ASSERT_EQUAL(T / s, mat4(-4, -4, 6, 0, 0, 2, -2, 0, 1, 0, -1, 0, 0, 0, 0, 2));
    ASSERT_EQUAL(s + Q, mat4(-0.5, -1.5, -1.5, 0.5, -0.5, -0.5, -3.5, 0.5, -0.5, -1.5, -3.5, 0.5, 0.5, 0.5, 0.5, 1.5));
    ASSERT_EQUAL(s - Q, mat4(1.5, 2.5, 2.5, 0.5, 1.5, 1.5, 4.5, 0.5, 1.5, 2.5, 4.5, 0.5, 0.5, 0.5, 0.5, -0.5));
    ASSERT_EQUAL(s * Q, mat4(-0.5, -1, -1, 0, -0.5, -0.5, -2, 0, -0.5, -1, -2, 0, 0, 0, 0, 0.5));
    ASSERT_EQUAL(s / (Q + s), mat4(-1, -1/3.0, -1/3.0, 1, -1, -1, -0.142857, 1, -1, -1/3.0, -0.142857, 1, 1, 1, 1, 1/3.0));

    //
    // Functions
    //
    ASSERT_EQUAL(radians(180.0), PI);
    ASSERT_EQUAL(degrees(PI), 180.0);
    ASSERT_EQUAL(sin(PI), 0.0);
    ASSERT_EQUAL(cos(0.0), 1.0);
    ASSERT_EQUAL(tan(PI/4.0), 1.0);
    ASSERT_EQUAL(tan(PI/6.0), sqrt(3.0)/3.0);
    ASSERT_EQUAL(atan(sqrt(3.0)/3.0), PI/6.0);
    ASSERT_EQUAL(pow(3.0, 0.5), sqrt(3.0));
    ASSERT_EQUAL(pow(PI, 2.0), PI*PI);
    ASSERT_EQUAL(exp(1.0), EULER);
    ASSERT_EQUAL(exp(0.0), 1.0);
    ASSERT_EQUAL(exp2(PI), exp(PI*LN2));
    ASSERT_EQUAL(log(EULER), 1.0);
    ASSERT_EQUAL(log2(256.0), 8.0);
    ASSERT_EQUAL(inversesqrt(PI*PI), 1.0/PI);
    ASSERT_EQUAL(abs(+PI), PI);
    ASSERT_EQUAL(abs(-PI), PI);
    ASSERT_EQUAL(sign(+PI), +1.0);
    ASSERT_EQUAL(sign(-PI), -1.0);
    ASSERT_EQUAL(floor(PI), 3.0);
    ASSERT_EQUAL(ceil(PI), 4.0);
    ASSERT_EQUAL(fract(PI), PI - floor(PI));
    ASSERT_EQUAL(mod(PI, 3.0), PI - 3.0);
    ASSERT_EQUAL(min(PI, 3.0), 3.0);
    ASSERT_EQUAL(max(PI, 3.0), PI);
    ASSERT_EQUAL(clamp(PI, 2.0, 3.0), 3.0);
    ASSERT_EQUAL(mix(PI, LN2, 0.5), (PI+LN2)*0.5);
    ASSERT_EQUAL(step(3.0, PI), 1.0);
    ASSERT_EQUAL(step(4.0, PI), 0.0);
    ASSERT_EQUAL(smoothstep(0.3, 0.6, 0.4), 0.259259);

    ASSERT_EQUAL(radians(180.0f*u), pi*u);
    ASSERT_EQUAL(degrees(v*pi),     v*180.0f);
    ASSERT_EQUAL(sin(v*pi),         vec2(0));
    ASSERT_EQUAL(cos(vec2(0)),      vec2(1));
    ASSERT_EQUAL(tan(vec2(PI, 0)/4.0f), vec2(1, 0));
    ASSERT_EQUAL(atan(vec2(1, 0)),  vec2(PI/4, 0));
    ASSERT_EQUAL(pow(u, s),         vec2(1.0, sqrt(2.0)));
    ASSERT_EQUAL(exp(u),            vec2(euler, 7.38905609));
    ASSERT_EQUAL(exp2(u),           vec2(2,4));
    ASSERT_EQUAL(log(u),            vec2(0, 0.693147));
    ASSERT_EQUAL(log2(u),           vec2(0,1));
    ASSERT_EQUAL(inversesqrt(u*u),  1.0f/u);
    ASSERT_EQUAL(abs(u-v),          vec2(+2.0));
    ASSERT_EQUAL(abs(v-u),          vec2(+2.0));
    ASSERT_EQUAL(sign(v-u),         vec2(+1.0));
    ASSERT_EQUAL(sign(u-v),         vec2(-1.0));
    ASSERT_EQUAL(floor(u*pi),       3.0f*u);
    ASSERT_EQUAL(ceil(u*pi),        vec2(4,7));
    ASSERT_EQUAL(fract(u*pi),       u*pi - floor(u*pi));
    ASSERT_EQUAL(mod(u*pi, v),      vec2(pi-3.0f, 2.0f*pi - 4.0f));
    ASSERT_EQUAL(min(u*pi, 3.0f),   vec2(3));
    ASSERT_EQUAL(max(u*pi, 4.0f),   vec2(4, 2*pi));
    ASSERT_EQUAL(clamp(0.5f*u, 0.5f, 2.0f), vec2(0.5, 1.0));
    ASSERT_EQUAL(mix(u, v, 0.5),    (u+v)*0.5f);
    ASSERT_EQUAL(step(u, v/2.0f),   vec2(1, 0));
    ASSERT_EQUAL(smoothstep(u, v, v/2.0f), vec2(0.15625, 0));

    ASSERT_EQUAL(radians(180.0f*a), pi*a);
    ASSERT_EQUAL(degrees(b*pi),     vec3(720, 900.00006, 1080)); //b*180.0f);
    ASSERT_EQUAL(sin(b*pi),         vec3(0));
    ASSERT_EQUAL(cos(vec3(0)),      vec3(1));
    ASSERT_EQUAL(tan(vec3(PI, 0)/4.0f), vec3(1, 0, 0));
    ASSERT_EQUAL(atan(vec3(1, 0)),  vec3(PI/4, 0, 0));
    ASSERT_EQUAL(pow(a, s),         vec3(1.0, sqrt(2.0), 1.73205));
    ASSERT_EQUAL(exp(a),            vec3(euler, 7.3890561, 20.085537));
    ASSERT_EQUAL(exp2(a),           vec3(2, 4, 8));
    ASSERT_EQUAL(log(a),            vec3(0, 0.69314718, 1.0986123));
    ASSERT_EQUAL(log2(a),           vec3(0, 1, 1.5849625));
    ASSERT_EQUAL(inversesqrt(a*a),  1.0f/a);
    ASSERT_EQUAL(abs(a-b),          vec3(3));
    ASSERT_EQUAL(abs(b-a),          vec3(3));
    ASSERT_EQUAL(sign(b-a),         vec3(+1.0));
    ASSERT_EQUAL(sign(a-b),         vec3(-1.0));
    ASSERT_EQUAL(floor(a*pi),       3.0f*a);
    ASSERT_EQUAL(ceil(a*pi),        vec3(4, 7, 10));
    ASSERT_EQUAL(fract(a*pi),       a*pi - floor(a*pi));
    ASSERT_EQUAL(mod(a*pi, b),      vec3(pi, 2.0f*pi - 5.0f, 3.424778));
    ASSERT_EQUAL(min(a*pi, 3.0f),   vec3(3));
    ASSERT_EQUAL(max(a*pi, 4.0f),   vec3(4, 2*pi, 9.424778));
    ASSERT_EQUAL(clamp(0.5f*a, 0.5f, 2.0f), vec3(0.5, 1.0, 1.5));
    ASSERT_EQUAL(mix(a, b, 0.5),    (a+b)*0.5f);
    ASSERT_EQUAL(step(a, b/2.0f),   vec3(1, 1, 0));
    ASSERT_EQUAL(smoothstep(a, b, b/2.0f), vec3(0.25925928, 0.074074082, 0));

    //vec4 r = pow(vec4(2.0), p);
    ASSERT_EQUAL(radians(180.0f*p), pi*p);
    ASSERT_EQUAL(degrees(q*pi),     q*180.0f);
    ASSERT_EQUAL(sin(q*pi),         vec4(0.0));
    ASSERT_EQUAL(cos(vec4(0.0)),    vec4(1.0));
    ASSERT_EQUAL(tan(vec4(pi/4)),   vec4(1.0))
    ASSERT_EQUAL(atan(vec4(1.0)),   vec4(pi/4.0));
    //ASSERT_EQUAL(pow(p, s),         vec4(1.0, sqrt(2.0), 1.7320508, 2));
    //ASSERT_EQUAL(exp(p),            vec4(euler, 7.3890561, 20.085537, 0)); // SKIP
    ASSERT_EQUAL(exp2(p),           pow(vec4(2.0), p)); //vec4(2, 4, 8, 0));
    //ASSERT_EQUAL(log(p),            vec4(0, 0.69314718, 1.0986123, 0));
    ASSERT_EQUAL(log2(pow(vec4(2.0), p)),           p);
    ASSERT_EQUAL(inversesqrt(p*p),  1.0f/p);
    ASSERT_EQUAL(abs(p-q),          vec4(4));
    ASSERT_EQUAL(abs(q-p),          vec4(4));
    ASSERT_EQUAL(sign(q-p),         vec4(+1.0));
    ASSERT_EQUAL(sign(p-q),         vec4(-1.0));
    ASSERT_EQUAL(floor(p*pi),       3.0f*p);
    ASSERT_EQUAL(ceil(p*pi),        3.0f*p + 1.0f);
    ASSERT_EQUAL(fract(p*pi),       p*pi - floor(p*pi));
    ASSERT_EQUAL(mod(p*pi, euler),  p*pi - floor(p*pi/euler)*euler);
    ASSERT_EQUAL(min(p*pi, pi),     vec4(pi));
    ASSERT_EQUAL(max(p*pi, p.y*pi), vec4(p.y, p.yzw)*pi);
    ASSERT_EQUAL(clamp(0.5f*p, 0.5f, 2.0f), vec4(0.5, 1.0, 1.5, 2.0));
    ASSERT_EQUAL(mix(p, q, 0.5),    (p + q)*0.5f);
    ASSERT_EQUAL(step(p, q*0.5f),   vec4(1, 1, 1, 0));
    ASSERT_EQUAL(smoothstep(p, q, q/2.0f), vec4(0.31640625, 0.15625, 0.04296875, 0));

    //
    // Geometric functions
    //
    ASSERT_EQUAL(length(PI), PI);
    ASSERT_EQUAL(length(u), 2.236068f);
    ASSERT_EQUAL(length(a), 3.7416575f);
    ASSERT_EQUAL(length(p), 5.4772258f);

    ASSERT_EQUAL(distance(PI, LN2), PI - LN2);
    ASSERT_EQUAL(distance(u, v), 2.8284271f);
    ASSERT_EQUAL(distance(a, b), 5.1961522f);
    ASSERT_EQUAL(distance(p, q), 8.0f);

    ASSERT_EQUAL(dot(PI, LN2), PI*LN2);
    ASSERT_EQUAL(dot(u, v), 11.0f);
    ASSERT_EQUAL(dot(a, b), 32.0f);
    ASSERT_EQUAL(dot(p, q), 70.0f);

    ASSERT_EQUAL(cross(a, b), vec3(-3, 6, -3));

    ASSERT_EQUAL(normalize(PI), 1.0);
    ASSERT_EQUAL(normalize(u), vec2(0.44721359, 0.89442718));
    ASSERT_EQUAL(normalize(a), vec3(0.26726124, 0.53452247, 0.80178368));
    ASSERT_EQUAL(normalize(p), vec4(0.18257418, 0.36514837, 0.54772258, 0.73029673));

    ASSERT_EQUAL(faceforward(PI, -1.0, LN2), +PI );
    ASSERT_EQUAL(faceforward(v, u, -u), vec2(3, 4) );
    ASSERT_EQUAL(faceforward(a, vec3(+1.0), b), -a );
    ASSERT_EQUAL(faceforward(p, -vec4(1.0), q), +p );

    ASSERT_EQUAL(reflect(PI, LN2), 0.122817);
    ASSERT_EQUAL(reflect(u, v), vec2(-65, -86));
    ASSERT_EQUAL(reflect(a, b), vec3(-255, -318, -381));
    ASSERT_EQUAL(reflect(p, q), vec4(-699, -838, -977, -1116));

//    ASSERT_EQUAL(refract(PI, LN2, 0.5), -0.1482120);
    ASSERT_EQUAL(refract(u, v, 0.5), vec2(-32.703293, -43.271057));
    ASSERT_EQUAL(refract(a, b, 0.5), vec3(-127.59369, -159.11711, -190.64053));
    ASSERT_EQUAL(refract(p, q, 0.5), vec4(-349.55356, -419.06427, -488.57498, -558.08569));

    ASSERT_EQUAL(M * N, mat2(1));
    ASSERT_EQUAL(A * B, mat3(1));
    ASSERT_EQUAL(T * Q, mat4(1));

    // Quaternions



    // Vector Pack/unpack




    // TODO more examples: 
    // https://gist.github.com/patricknelson/f4dcaedda9eea5f5cf2c359f68aa35fd
    // https://gist.github.com/rlane/1223480
    // https://en.wikibooks.org/wiki/GLSL_Programming/Vector_and_Matrix_Operations

    // TODO performance tests
}