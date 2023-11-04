#pragma once

#include "type/cpu.h"
#include "type/array.h"

#include "__temp/texture.h"


// PROP 
// struct IndexBuffer : Array<ivec2>, Cacheable


struct IndexBuffer : public Cacheable
{
	int index_size;
	uint8* array;
	int length;
	int allocated;
	int granularity;

	int lods_count;
	struct LOD 
	{
		int start;
		int count;
	} *lods;




	int elements()
	{
		return lods[lods_count-1].start + lods[lods_count-1].count;
	}

	int bytes()
	{
		return elements() * index_size;
	}

	
	IndexBuffer( int allocated=0, int granularity=1 ) :
		index_size(0), array(NULL), length(0), allocated(allocated), granularity(granularity),
		lods_count(0), lods(NULL)
	{
		if( allocated ) resize( allocated );
	}

	~IndexBuffer() 
	{
		delete[] array;
		delete[] lods;
	}

	void clear()
	{
		//index_size = 0;
		delete[] array;
		array = NULL;
		length = 0;
		allocated = 0;
		delete[] lods;
		lods = NULL;
	}
 
	int indexSize()
	{
		return index_size;
	}

	int size()
	{
		return length;
	}

	IndexBuffer& resize( int size, int granularity=0 )
	{
		if( allocated!=size )
		{
			array = (uint8*) realloc( array, size*index_size );
			allocated = size;
		}
		if( granularity>0 )
		{
			this->granularity = granularity;
		}
		return *this;
	}

	inline unsigned int get( int i )
	{
		return index_size==1 ? ((uint8*)array)[i] : ((uint16*)array)[i];
	}

	inline void set( unsigned int i, unsigned int value )
	{
		if( index_size==1 )
		{
			((uint8*)array)[i] = value;
		}
		else
		{
			((uint16*)array)[i] = value;
		}
	}

	int operator[]( int i ) 
	{ 
		return get(i);
	}

	void push( int value )
	{
		if( length+1>allocated ) resize( allocated+granularity );
		set( length++, value );
	}
/*
	const T& pop()
	{
		const T& value = array[--length];
		if( length<=allocated-granularity ) resize( allocated-granularity );
		return value;
	}
*/	
	IndexBuffer& fit()
	{
		return resize(length);
	}


	//
	// quads for GL4Profile
	//
	void init( int tileRes )
	{
		clear();
		index_size = tileRes<=16 ? 1 : 2;
		lods_count = 1;
		lods = new IndexBuffer::LOD[1];
		lods->start = 0;
		lods->count = (tileRes-1)*(tileRes-1)*4; // no strips, no stickers here, just quads
		resize( lods->count );

		for( int j=0; j<tileRes-1; j++ )
		for( int i=0; i<tileRes-1; i++ )
		{
			push( (i+1) + (j+0)*tileRes );
			push( (i+0) + (j+0)*tileRes );
			push( (i+0) + (j+1)*tileRes );
			push( (i+1) + (j+1)*tileRes );	
		}
		videomem_invalidated = true;
		hostmem_invalidated = false;
	}


	

};
