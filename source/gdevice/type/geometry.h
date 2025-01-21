#pragma once

#include "type/glsl.h"

inline float triangleArea( const vec3& v1, const vec3& v2, const vec3& v3 )
{
	float a = distance(v1,v2);
	float b = distance(v2,v3);
	float c = distance(v3,v1);
	float p = (a+b+c)/2;
	return sqrt(p*(p-a)*(p-b)*(p-c)); // Thank you Heron
}

inline float quadrilateralArea( const vec3& v1, const vec3& v2, const vec3& v3, const vec3& v4 )
{
	return triangleArea(v1,v2,v3) + triangleArea(v3,v4,v1);
}



