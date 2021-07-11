#pragma once

#include "libraries/math/algebra.h"

// TODO: Use the code from asset.h
// TODO: Make a "procedural material" that works with 4 patterns, mixmap and color.
// TODO: Make a simple OpenGL material? (that acts as a base class?)

struct Material
{
	vec4 ambient;
	vec4 diffuse;	// Not used at the moment // TEMP1
	vec4 emission;
	vec4 specular;
	float shininess;

	// textures
	// uniforms
	vec4 specular_weights;
};