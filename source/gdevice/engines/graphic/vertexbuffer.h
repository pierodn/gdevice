#pragma once

#include "libraries/math/algebra.h"

#include "cacheable.h"
#include "texture.h"

struct VertexBuffer : Cacheable // this will be used as either VBO or FBO (not both!).
{
	Texture<vec4> quartets;
	Texture<vec4> gradients;
	Texture<vec4> colors;
	Texture<vec4> mixmaps;

	void init(int tileRes)
	{
		quartets .resample<0>(tileRes, tileRes);
		gradients.resample<0>(tileRes, tileRes);
		colors	 .resample<0>(tileRes, tileRes);
		mixmaps  .resample<0>(tileRes, tileRes);
    }
	
	int bytes()
	{
		return quartets.bytes() + gradients.bytes() + colors.bytes() + mixmaps.bytes(); 
	}
};