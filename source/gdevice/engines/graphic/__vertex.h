#pragma once

#include "libraries/math/algebra.h"


// DONE interleave for best performance 
// TODO also pad to 32 bytes ( 12+12+4+4 = 32 )
// http://www.opengl.org/wiki/Vertex_Arrays

// TODO position as vec4(x,y,z,zc)
// Zc is the coarser clipmap Z of the vertex. 
// It is used only by programmed function to smooth Z and avoid either poppings and cracks.
// Fixed function will only rely on a flat ground-colored panel and show-backfaces mode.

// TODO normal as vec3(dx,dy,1)
// So that will normalized when in fixed function 
// or passed as just vec2(dx,dy) in case of programmed function


struct Vertex 
{
	vec3 position;		
	vec3 normal;
	rgba color;
	rgba mixmap; 
	vec4 quartet;	// used by GL4 function (the only attribute used)
};
