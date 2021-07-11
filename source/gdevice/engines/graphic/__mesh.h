#pragma once


// TODO use this as remind for the mesh algebra (which is now a plain texture algebra)

struct Mesh : public Cacheable 
{
	//
	// Geometric operations
	//
/*
	Mesh& scale( const vec3& a )
	{
		for( int y=0; y<positions.height; y++ )
		for( int x=0; x<positions.width;  x++ ) 
		{
			positions.at(x,y) *= a;
		}

		return *this;
	}

	Mesh& rotate( float ax, float ay, float az )
	{
		for( int j=0; j<positions.height; j++ )
		for( int i=0; i<positions.width;  i++ )
		{
			positions.at(i,j) = rot(ax,ay,az) * positions.at(i,j);
		}
		return *this;
	}

	Mesh& translate( const vec3& a )
	{
		for( int y=0; y<positions.height; y++ )
		for( int x=0; x<positions.width;  x++ ) 
			positions.at(x,y) += a;

		return *this;
	}

	Mesh& mix( Mesh& mesh, float alpha )
	{
		for(int x=0; x<positions.width;  x++)
		for(int y=0; y<positions.height; y++)
		{
			float u = (float)x*mesh.positions.width/positions.width;
			float v = (float)y*mesh.positions.height/positions.height;
			vec3 v1 = positions[y][x];
			vec3 v2 = mesh.positions.interpolate<coserp>( u, v);
			positions.at(x,y) = lerp( alpha, v1, v2 );
		}
		return *this;
	}

	Mesh& tubemuly( float (*ft)(float), float (*fh)(float)=unit, float (*fr)(float)=unit )
	{
		for( int j=0; j<positions.height; j++ )
		for( int i=0; i<positions.width;  i++ )
		{
			float x = (float)i/positions.width;
			float y = (float)j/positions.height;
			
			vec3 p = cartesian_to_cylindrical( positions.at(i,j) );
			p.theta  *= ft(y);
			p.height *= fh(y);
			p.radius *= fr(y);
			vec3 q = cylindrical_to_cartesian( p );
			positions.at(i,j) = q;
		}
		return *this;
	}

	Mesh& tubemul( float (*ft)(float,float), float (*fh)(float,float), float (*fr)(float,float) )
	{
		for( int j=0; j<positions.height; j++ )
		for( int i=0; i<positions.width;  i++ )
		{
			float x = (float)i/positions.width;
			float y = (float)j/positions.height;
			
			vec3 p = cartesian_to_cylindrical( positions.at(i,j) );
			p.theta  *= ft(x,y);
			p.height *= fh(x,y);
			p.radius *= fr(x,y);
			vec3 q = cylindrical_to_cartesian( p );
			positions.at(i,j) = q;	
		}
		return *this;
	}


	template<float* profile>
	Mesh& mould_y()
	{
		return tubemul( unit, unit, add<zero,cubicinterpolate<profile>> );
	}



///
	Mesh& test_args( double x0, double y0, ... )
	{
		va_list argptr;
		double f;

		va_start(argptr, y0);
		do {
			f = va_arg(argptr, double);
			printf("arg = %f\n",f);
		} while (f<1.0);
		va_end(argptr);

		return *this;
	}

	Mesh& test_mould_y( double x0, double y0, ... )
	{
		// ---x0---x1-x-x2---x3---

		va_list argptr;
		va_start(argptr, y0);

		double x1 = x0;
		double y1 = y0;
		double x2 = va_arg(argptr, double);
		double y2 = va_arg(argptr, double);
		double x3 = va_arg(argptr, double);
		double y3 = va_arg(argptr, double);

		for( int j=0; j<positions.height; j++ )
		{
			float y = (float)j/positions.height;   // the 'x' to be interpolated 

			while( y>x2 && x2<1.0 )
			{
				x0 = x1; 
				x1 = x2; 
				x2 = x3;
				y0 = y1; 
				y1 = y2; 
				y2 = y3;
				if( x3<1.0 )
				{
					x3 = va_arg(argptr, double);
					y3 = va_arg(argptr, double);
				}
			}

			double fy = cubic( (y-x1)/(x2-x1), y0, y1, y2, y3 );

			for( int i=0; i<positions.width;  i++ )
			{
				float x = (float)i/positions.width;
								
				vec3 p = cartesian_to_cylindrical( positions.at(i,j) );
				//p.theta  *= ft(x,y);
				//p.height *= fh(x,y);
				p.radius *= fy;
				vec3 q = cylindrical_to_cartesian( p );
				positions.at(i,j) = q;	
			}
		}

		return *this;
	}

	Mesh& test2_mould_y( float* s )
	{
		float x0 = *s++;
		float y0 = *s++;
		float x1 = x0;
		float y1 = y0;
		float x2 = *s++;
		float y2 = *s++;
		float x3 = *s++;
		float y3 = *s++;

		for( int j=0; j<positions.height; j++ )
		{
			float y = (float)j/positions.height;   // the 'x' to be interpolated 

			while( y>x2 && x2<1.0 )
			{
				x0 = x1; 
				x1 = x2; 
				x2 = x3;
				y0 = y1; 
				y1 = y2; 
				y2 = y3;
				if( x3<1.0 )
				{
					x3 = *s++;
					y3 = *s++;
				}
			}

			double fy = cubic( (y-x1)/(x2-x1), y0, y1, y2, y3 );

			for( int i=0; i<positions.width;  i++ )
			{
				float x = (float)i/positions.width;
								
				vec3 p = cartesian_to_cylindrical( positions.at(i,j) );
				//p.theta  *= ft(x,y);
				//p.height *= fh(x,y);
				p.radius *= fy;
				vec3 q = cylindrical_to_cartesian( p );
				positions.at(i,j) = q;	
			}
		}

		return *this;
	}
////

	Mesh& twist_y( float (*fx)(float), float (*fy)(float), float (*fz)(float) )
	{
		for( int j=0; j<positions.height; j++ )
		for( int i=0; i<positions.width;  i++ )
		{
			float y = (float)j/positions.height;
			positions.at(i,j) = rot(fx(y),fy(y),fz(y)) * positions.at(i,j);
		}

		return *this;
	}

	// TEMP to test
	template< vec3 interpolator(float,float,const vec3&,const vec3&,const vec3&,const vec3&) >
	Mesh& resample( int width, int height )
	{
		positions.resample<interpolator>( width, height );
		return *this;
	}
*/



	//

	Mesh& join( Mesh& mesh, float glue=0.0 );

	//
	// mapping operations
	// 

	Mesh& add( Mesh& mesh );
	Mesh& sub( Mesh& mesh );
	Mesh& mul( Mesh& mesh );
	Mesh& div( Mesh& mesh );
	Mesh& rotate( int dx, int dy );


};
