#pragma once

// ---- Premise ----
// gdevice is not meant to be a generic 3D engine.
// It heavily addresses procedural content and comes with
// specific object types and algebras designed to promote
// orders of manipulation.
// However, conversion facilities will allow to inject
// generic content from outside, that thus can be profitably
// used for further manipulation enabling composition of
// non-generic content with generated content.


// ---- Object semantics ----
// A Texture<T> object holds an array of T objects
// and can profitably embody graphic objects such as:
//	* image objects     (textures, PBOs, FBOs..)  
//	* geometric objects (2D VBOs, generic VBOs..)
//
// Examples: 
//		texture<texel>  texture;
//		texture<pixel>  fbo; 
//		texture<vertex> vbo;
//		texture<float>  terrain;
//
// Note: the one who couldn't be (yet) fitted into a texture
// is the index buffer in combination with its lod segments
// because it is required to be like an "array of arrays"
// within a monolitic buffer.
// Actually a texture can be seen as a array of arrays, but
// that implementation would bring more complicated handling
// and probably worse performances.


// ---- Dimensions ----
// Textures have deferred constructors that give
// a dimensional feature still retaining a promiscuous access.
// Examples: 
//		vbo.resample(20,30);	// shapes vbo as a 2D image
//		x = vbo.at(3);			// accesses vbo as a linear buffer


// ---- Algebra semantics ----
// Textures are subject to a uniform algebra (add,lerp,etc.)
// that enables to compose and manipulate images, geometries 
// and virtually every kind of T objects who are compliant to 
// the algebra.
// Since this algebra has a preeminent role for manipulation,
// those objects who are compliant to the algebra are thus
// considered (from a procedural standpoint) as 1st class objects
// and their names are not capitalized.
// Examples: pixel, texel, vertex, voxel..


// ---- Details ----
// Texture 1D has width>0, height=0, depth=0
// Texture 2D has width>0, height>0, depth=0
// Texture 3D has width>0, height>0, depth>0


#include "gpu/cacheable.h"


template<class T> 
struct Texture : public Cacheable
{
	// TODO when this leaves
	// Cacheable vbo;	// when binded with GL::Texturing::bind() (possibily PBO-copying from TBO counterpart, if exists)
	// Cacheable tbo;	// when binded with VBO::bind() (possibily PBO-copying from VBO counterpart, if exists)

	T*		array;		
	ivec3	origin;		// TODO Toroidal access
	ivec3	size;
	vec3	scale;
	bvec3	wrap;

    Texture() : array(NULL), origin(0), size(0), scale(1), wrap(true) {
	}

    Texture(char* filename) : array(NULL), origin(0), size(0), scale(1), wrap(true) {
        IO::load<T>(*this, filename);
	}

    void store(char* filename) {
        //IO::saveRGB<T>(*this, filename);
        IO::save<T>(*this, filename);
    }

    ~Texture() {
        resample<0>(0,0);
    }

	int elements() {
		return (size.x>0 ? size.x : 1) * (size.y>0 ? size.y : 1) * (size.z>0 ? size.z : 1);
	}

	int bytes() {
		return elements() * sizeof(T);
	}	

	bool empty() {
		return elements()>0;
	}

	//
	// Deferred constructors
	//
	template<T interpolator(float,float,const T&,const T&,const T&,const T&)>
	Texture<T>& resample( int width, int height )
	{
        if( height==0 ) {
            height = width;
        }

		if( size.x == width && size.y == height ) {
			return *this;
		}

		T* temp = NULL;
		if( width*height>0 ) {
			temp = new T[width*height];
            if( interpolator!=NULL ) {
				for(int x=0; x<width;  x++)
                for(int y=0; y<height; y++) {
					temp[ y*width+x ] = interpolate<interpolator>( (float)x*size.x/width, (float)y*size.y/height );
                }
            }
		}	

		if( array ) delete array;
		array = temp;
		origin = 0;
		size.x = width;
		size.y = height;
		return *this;
	}

	// copy constructor
	template<T interpolator(float,float,const T&,const T&,const T&,const T&)>
	void copy( Texture<T>& t )
	{
		//printf("--- copy invoked\n"); //
		float kx = t.size.width/size.width;
		float ky = t.size.height/size.height;

		for( int x=0; x<size.width;  x++ )
		for( int y=0; y<size.height; y++ )
			at(x,y) = t.interpolate<interpolator>(x*kx, y*ky);
	}

	Texture<T>& operator=(const Texture<T>& t) 
	{
		if( this!= &t ) copy(t);
		return *this;
	}

	template< T interpolator(float,float,const T&,const T&,const T&,const T&) >
	inline T interpolate( float x, float y )
	{
		int i = (int)x;
		int j = (int)y;
		float ax = x - (float)i;
		float ay = y - (float)j;
		return interpolator( ax, ay, at(i,j), at(i,j+1), at(i+1,j), at(i+1,j+1) );
	}
/*
	// TEMP
	float pick( vec2 point )
	{
		int res = size.height - 1;
		float u = point.x*res, v = point.y*res;
		//vec3 ve = vertex_buffer->positions.interpolate<lerp>( u, v );
		//return ve.z;
		
		int i = (int)u;
		int j = (int)v;
		float ax = u - (float)i;
		float ay = v - (float)j;
		float z0 = at(i,j).position.z;
		float z1 = at(i,j+1).position.z;
		float z2 = at(i+1,j).position.z;
		float z3 = at(i+1,j+1).position.z;
		return lerp( ax, ay, z0, z1, z2, z3 );
	}
*/

	inline int index( int x )
	{
		// BUG this line makes normals not be shown (F2)
		//x = x % size.x  + (x<0 ? size.x : 0); 
		return x;
	}

	inline int index( int x, int y )
	{
		x = x%size.x + (x<0 ? size.x : 0); 
		y = y%size.y + (y<0 ? size.y : 0); 
		return x + y*size.x;
	}

	inline int index(int x, int y, int z)
	{
		x = x%size.x  + (x<0 ? size.x : 0); 
		y = y%size.y  + (y<0 ? size.y : 0); 
		z = z%size.z  + (z<0 ? size.z : 0); 
		return x + y*size.x + z*size.x*size.y;
	}

	inline T& at(int x)				    { return array[index(x)]; }
	inline T& at(int x, int y)		    { return array[index(x,y)]; }
	inline T& at(int x, int y, int z)	{ return array[index(x,y,z)]; }

    // Derivative operators
    inline T forwardDx(int x, int y)    { return at(x+1,y) - at(x, y); }
	inline T forwardDy(int x, int y)    { return at(x,y+1) - at(x, y); }
	inline T backwardDx(int x, int y)   { return at(x,y) - at(x-1, y); }
	inline T backwardDy(int x, int y)   { return at(x,y) - at(x, y-1); }
    inline T centralDx(int x, int y)    { return 0.5f*(at(x+1,y) - at(x-1, y)); }
	inline T centralDy(int x, int y)    { return 0.5f*(at(x,y+1) - at(x, y-1)); }

    // Derivative operators that integrate mild blurring:
    // Prewitt  => average blur
    // Sobel    => gaussian blur
    // https://theailearner.com/tag/scharr-operator/

    #define SAMPLE_NEIGHBORS  \
    	T C00 = at(x-1, y-1); \
		T C01 = at(x-1, y+0); \
		T C02 = at(x-1, y+1); \
		T C10 = at(x+0, y-1); \
		T C12 = at(x+0, y+1); \
		T C20 = at(x+1, y-1); \
		T C21 = at(x+1, y+0); \
		T C22 = at(x+1, y+1)

    inline T prewittDx(int x, int y) {
        SAMPLE_NEIGHBORS;
        return 0.333334f*(C20 + C21 + C22 - C00 - C01 - C02);
	}

	inline T prewittDy(int x, int y) {
        SAMPLE_NEIGHBORS;
		return 0.333334f*(C02 + C12 + C22 - C00 - C10 - C20);
	}

    inline T sobelDx( int x, int y ) {
		SAMPLE_NEIGHBORS;
        return 0.5f*(C21 - C01) + 0.25f*(C20 + C22 - C00 - C02); 
	}

    inline T sobelDy( int x, int y ) {
		SAMPLE_NEIGHBORS;
		return 0.5f*(C12 - C10) + 0.25f*(C02 + C22 - C00 - C20);  
	}

    
/*
    // Only for T == float

	inline vec2 sobelGradient( int x, int y )
	{
        SAMPLE_NEIGHBORS;
		float dx = 0.5f*(C21 - C01) + 0.25f*(C20 + C22 - C00 - C02);
		float dy = 0.5f*(C12 - C10) + 0.25f*(C02 + C22 - C00 - C20);  
		return vec2(dx, dy);
	}

	inline vec2 gradient(int x, int y) {
		return vec2(dx, dy);
	}

	inline vec3 normal( int x, int y ) {
		return normalize(vec3(gradient(x,y),1)); 
	}
*/

    //
    // Blur filters
    // 

    inline T box2x2(int x, int y) 
    {
        return 0.25f*(  at(x+0, y+0) + 
                        at(x+1, y+0) + 
                        at(x+0, y+1) + 
                        at(x+1, y+1) );  
    }

    // TODO: gaussian
};

//
// algebra
//

// PROP texture1 op texture2 => apply foreach item: texture1.item op texture2.item


/*
Texture<T>& lerp( float alpha, const Texture<T>& t1, const Texture<T>& t2 )
{
	Texture<T> t;

	if( t1.width*t1.height >  t2.width*t2.height )
		t.resample(t1.width,t1.height);
	else t.resample(t2.width,t2.height);

	for( int x=0; x<t.width;  x++ )
	for( int y=0; y<t.height; y++ )
	{
		float x2 = (float)x * t2.width/width;
		float y2 = (float)y * t2.height/height;
		T v1 = at(x,y);
		T v2 = t2.interpolate<coserp>(x2,y2);
		//at(x,y) = lerp( alpha, v1, v2 );
		at(x,y) = t2.interpolate<coserp>(x2,y2);
	}
}
*/

/*
template<class T> inline Texture<T> operator=( const Texture<T>& a, const Texture<T>& b )
{
	Texture<T> temp;
	// TODO
	return temp;
}
*/

/*
template<class T, int N> inline vec<T,N> operator+( const vec<T,N>& a, const vec<T,N>& b )
{ 
	vec<T,N> temp;
	for( int i=0; i<N; i++ ) temp[i] = a[i]+b[i];
	return temp;
}
*/