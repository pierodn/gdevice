#pragma once

#include "libraries/math/algebra.h"



vec3 polar_to_cartesian( const vec3& p )
{
	 return p.radius * vec3( 
		cos(p.theta) * sin(p.phi), 
					  -cos(p.phi), 
	   -sin(p.theta) * sin(p.phi) ); 
}

vec3 cartesian_to_polar( const vec3& p )
{
	float radius = length(p);
	float theta  = acos( p.y/radius );
	float phi    = atan2( p.z, p.x );
	return vec3( theta, phi, radius );
}

vec3 cylindrical_to_cartesian( const vec3& p )
{ 
	return vec3( 
		 p.radius * cos(p.theta),
		 p.height - 0.5,
		-p.radius * sin(p.theta));
}

vec3 cartesian_to_cylindrical( const vec3& p )
{
	float theta  = atan2( -p.z, p.x );
	float radius = sqrt( p.x*p.x + p.z*p.z );
	float height = p.y + 0.5;
	return vec3( theta, height, radius );
}

vec3 conical_to_cartesian( const vec3& p )
{
	return vec3(
		p.radius*(1-p.y)*cos(p.theta),
		p.height - 0.5,
		p.radius*(p.y-1)*sin(p.theta));
}

vec3 cartesian_to_conical( const vec3& p )
{
	// TODO
	return vec3(0);
}

vec3 toroidal_to_cartesian( const vec3& p, float hole = 0.5f )
{
	return vec3(
		cos(p.theta) * (cos(p.phi)*p.radius + hole),
		sin(p.theta) * (cos(p.phi)*p.radius + hole),
		sin(p.phi)*p.radius);
}

vec3 cartesian_to_toroidal( const vec3& p, float hole = 0.5f )
{
	// TODO
	return vec3(0);
}

