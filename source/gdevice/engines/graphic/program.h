#pragma once

#include "engines/graphic/array.h"
#include "engines/graphic/cacheable.h"


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
