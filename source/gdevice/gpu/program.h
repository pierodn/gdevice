#pragma once


#include "type/array.h"
#include "gpu/cacheable.h"
#include "gpu/opengl/gl.h"

struct Program : public Array<Cacheable>, Cacheable
{
    Program()
	{
        id = 0;
	}

	~Program() 
	{ 
		for( int i=0; i<length; i++ )
		{
			array[i].deallocate();
		}
	}

    // private
    static void deallocateProgram( uint& id )
    {
        if( id ) glDeleteProgram(id);
        id = 0;
    }
    static void deallocateShader( uint& id )
    {
        if( id ) glDeleteShader(id);
        id = 0;
    }

    void Build(char* source)
    {
        ASSERT(source != NULL);
        ASSERT(strlen(source) < MAX_PROGRAM_SOURCE_LENGTH);
        ASSERT(id == 0);

        id = glCreateProgram();  // TODO GL::CreateProgram();
        deallocator = &Program::deallocateProgram;

        const int ShaderTypesCount = 10;
        int shaders[ShaderTypesCount];
        int shaderCount = GL::BuildShaders(source, shaders);
        
        for(int i=0; i<shaderCount; i++)
        {
            ASSERT(shaders[i]);

            Cacheable& shader = push().tail();
            shader.id = shaders[i];
            shader.deallocator = &Program::deallocateShader;        
            GL::AttachShaderToProgram(id, shader.id ); 
        }

        GL::LinkProgram(id);
    }
    
    template<class T> 
    void SetUniform(char* name, T value)
    {
        ASSERT(id);
        // TODO ASSERT program is bound
        //ASSERT(GetCurrentProgramId() == id); // FIX
        GL::SetUniform(id, name, value);
    }

    void Use()
    {
        GL::UseProgram(id);
    }

    void Run()
    {
        ASSERT(id);

/*
        glPatchParameteri(GL_PATCH_VERTICES, 4);		
	    glDrawElements(
		    GL_PATCHES, 
		    ibo.lods[lod].count, // It is 1024 for tiles 17x17
		    GL_UNSIGNED_SHORT, 
		    0
	    ); 
*/
    }
};
