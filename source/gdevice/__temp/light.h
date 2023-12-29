#pragma once

#include "type/gpu.h"

struct Light
{
	int id;

	vec4 ambient;
	vec4 diffuse;
	vec4 specular;

	union {
		struct { vec3 spot_direction; float positional; };
		struct { vec4 position; };
	};
	int   spot_exponent;
	int   spot_cutoff;
	float constant_attenuation;
	float linear_attenuation;
	float quadratic_attenuation;

	// precomputations (multiply by current material)
	vec4 FrontLightProduct_ambient; //KaLa;
	vec4 FrontLightProduct_specular; //KsLs;
};