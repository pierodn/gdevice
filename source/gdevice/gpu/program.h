#pragma once

#include "type/array.h"
#include "gpu/cacheable.h"


struct Program : public Array<Cacheable>, Cacheable
{
	char* source;

	Program( char* source = NULL )
	{
		this->source = source;
	}

	~Program() 
	{ 
		for( int i=0; i<length; i++ )
		{
			array[i].deallocate();
		}
	}
};
