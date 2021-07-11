#pragma once

// PROP

struct Shaper
{
	enum Primitive
	{
		FlatPrimitive,
		TubularPrimitive,
		ConicalPrimitive,
		SpheroidalPrimitive,
		ToroidalPrimitive
	};

	static inline Texture<vec3> genWall( int width, int height=0, float xsmooth = 0, float ysmooth = 0 )
	{
		return genPrimitive( width, height, FlatPrimitive, xsmooth, ysmooth );
	}

	static inline Texture<vec3> genTube( int width, int height=0, float xsmooth = 0, float ysmooth = 0 )
	{
		return genPrimitive( width, height, TubularPrimitive, xsmooth, ysmooth );
	}

	static inline Texture<vec3> genCone( int width, int height=0, float xsmooth = 0, float ysmooth = 0 )
	{
		return genPrimitive( width, height, ConicalPrimitive, xsmooth, ysmooth );
	}

	static inline Texture<vec3> genBall( int width, int height=0, float xsmooth = 0, float ysmooth = 0 )
	{
		return genPrimitive( width, height, SpheroidalPrimitive, xsmooth, ysmooth );
	}

	static inline Texture<vec3> genRing( int width, int height=0, float xsmooth = 0, float ysmooth = 0 )
	{
		return genPrimitive( width, height, ToroidalPrimitive, xsmooth, ysmooth );
	}

	static inline Texture<vec3> genPrism( int width, int height=0, float xsmooth = 0, float ysmooth = 0 )
	{
		return genPrimitive( width, height, TubularPrimitive, xsmooth, ysmooth );
	}

	static inline Texture<vec3> genCube( float xsmooth = 0, float ysmooth = 0 )
	{
//		return genPrimitive( 4, 2, SpheroidalShape, xsmooth, ysmooth ).scale( vec3(1.4,1.0,1.4) );
	}

	static Texture<vec3> genPrimitive( int width, int height, Primitive shape = FlatPrimitive, 
		float xsmooth = 0, float ysmooth = 0, float radius = 1 )
	{
		Texture<vec3> positions;

		width = xsmooth==0 ? width : width*2;
		height = ysmooth==0 ? height : height*2;
		height = height>0 ? height : width<=128 ? 256/width : 65536/width;

		positions.resample<0>( width, height );
		positions.wrap_s = (primitive != FlatPrimitive); 
		positions.wrap_t = (primitive == ToroidalPrimitive);

		for( int j=0; j<positions.height; j++ )
		for( int i=0; i<positions.width;  i++ )
		{
			float x = xsmooth==0 ? (float)i/positions.width		  : (float)(i/2+(i%2)*xsmooth)/positions.width*2;
			float y = ysmooth==0 ? ((float)j+0.5)/positions.height : (float)(j/2+(0.5-ysmooth)+(j%2)*ysmooth)/positions.height*2;
			float z = 0.0; // TODO use a fz(x,y) as parameter, to map a texture-heights onto (like a bump)

			vec3 vertex = 
				primitive==FlatPrimitive       ?							vec3(x-0.5,  y-0.5,  0) : 
				primitive==TubularPrimitive    ? cylindrical_to_cartesian(	vec3(2*PI*x, y,      0.5)) :
				primitive==ConicalPrimitive    ? conical_to_cartesian(		vec3(2*PI*x, y,      0.5)) :
				primitive==SpheroidalPrimitive ? polar_to_cartesian(		vec3(2*PI*x, PI*y,   0.5)) :
				primitive==ToroidalPrimitive   ? toroidal_to_cartesian(		vec3(2*PI*x, 2*PI*y, 0.15), 0.35) :
				0; // TODO add primitives

			vertex *= radius;

			positions.at(i,j) = vertex;
		}
		return positions;
	}
};