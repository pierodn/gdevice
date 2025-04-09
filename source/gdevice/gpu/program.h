#pragma once

#include "type/array.h"
#include "gpu/opengl/gl.h"
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

    void Bind()
    {
        // TODO GL::GLSL::bind(*this);
    }

    template<class T> 
    void Set(char* name, T value)
    {
        GL::GLSL::set(*this, name, value);
    }
};
