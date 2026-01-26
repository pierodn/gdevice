#pragma once

//
// GLoom - GL Object Oriented Middleware
// 

#include "gpu/opengl/gluon.h"

#include "gpu/cacheable.h"
#include "__temp/Texture.h"
#include "__temp/VertexBuffer.h"
#include "gpu/IndexBuffer.h"

namespace GL
{ 
	inline char* vendor()
	{
		return (char*)glGetString(GL_VENDOR);
	}

    inline char* renderer()
	{
        return (char*)glGetString(GL_RENDERER);
    }

    inline char* error()
	{
		GLenum error = glGetError();
		return 
            error == GL_NO_ERROR                        ? NULL :
			error == GL_INVALID_ENUM		            ? "Invalid enum" :
			error == GL_INVALID_VALUE		            ? "Invalid value" :
			error == GL_INVALID_OPERATION               ? "Invalid operation" :
			error == GL_STACK_OVERFLOW                  ? "Stack overflow" :
			error == GL_STACK_UNDERFLOW                 ? "Stack underflow" :
			error == GL_OUT_OF_MEMORY                   ? "Out of memory" : 
            error == GL_INVALID_FRAMEBUFFER_OPERATION   ? "Invalid framebuffer operation" :
                                                          "Unknown error";
	}
}


namespace GL
{ 
    namespace VBO
    {

        void create( Cacheable& buffer )
        {
	        glGenBuffers( 1, &buffer.id );
	        buffer.deallocator = deallocate;
        }

        void bind( Cacheable& buffer, int type )//=GL_ARRAY_BUFFER )
        {
	        glBindBuffer( type, buffer.id );
        }

        void unbind( int type = GL_ARRAY_BUFFER)
        {
	        glBindBuffer( type, 0 );
        }

        void assert_allocated_is( int size, int type = GL_ARRAY_BUFFER)
        {
	        int allocated_size;
	        glGetBufferParameteriv( type, GL_BUFFER_SIZE, &allocated_size );
	        if( allocated_size != size ) 
	        {
		        //buffer.deallocate(); // OUT_OF_MEMORY ??
		        DEBUG_CRITICAL("Allocation fault");
	        }
        }

        void allocate( Cacheable& buffer, int size, void* data = NULL, int type = GL_ARRAY_BUFFER, int usage = GL_STATIC_DRAW )
        {
	        glBindBuffer( type, buffer.id );
	        glBufferData( type, size, data, usage );
	        assert_allocated_is( size, type );
	        if( data ) 
	        {
		        buffer.videomem_invalidated = false;
		        buffer.hostmem_invalidated = false; // it's supposed to be its data, already in the array 
	        }
        }

        void update( Cacheable& buffer, int size, void* data = NULL, int type = GL_ARRAY_BUFFER, int usage = GL_STATIC_DRAW )
        {
	        if( buffer.videomem_invalidated )
	        {
		        buffer.deallocate();
		        VBO::create( buffer );
		        VBO::allocate( buffer, size, data, type, usage );
	        }
        }

        template <class T>
        void bind( Texture<T>& vbo )
        {
	        update( vbo, vbo.bytes(), vbo.array, GL_ARRAY_BUFFER );
	        bind( vbo, GL_ARRAY_BUFFER );
        }


        // TODO: Handle abstract attributes...
        void bind( VertexBuffer& vbo )
        {
	        if( vbo.videomem_invalidated ) 
	        {
                if( !vbo.id )
                {
                    create( vbo );
		            allocate( vbo, vbo.bytes() );
                }
		        glBufferSubData( GL_ARRAY_BUFFER, 0, vbo.quartets.bytes(), vbo.quartets.array);
		        glBufferSubData( GL_ARRAY_BUFFER, vbo.quartets.bytes(), vbo.gradients.bytes(), vbo.gradients.array);
		        glBufferSubData( GL_ARRAY_BUFFER, vbo.quartets.bytes()+ vbo.gradients.bytes(), vbo.colors.bytes(), vbo.colors.array);
		        glBufferSubData( GL_ARRAY_BUFFER, vbo.quartets.bytes()+ vbo.gradients.bytes()+ vbo.colors.bytes(), vbo.mixmaps.bytes(), vbo.mixmaps.array);
                
                vbo.videomem_invalidated = false;
	        }

	        VBO::bind( vbo, GL_ARRAY_BUFFER );
        }

        void bind( IndexBuffer& ibo )
        {
            VBO::update( ibo, ibo.bytes(), ibo.array, GL_ELEMENT_ARRAY_BUFFER );
            VBO::bind( ibo, GL_ELEMENT_ARRAY_BUFFER );
        }

        void upload( Cacheable& buffer, int offset, int size, void* data, int type = GL_ARRAY_BUFFER )
        {
	        glBindBuffer( type, buffer.id );
	        glBufferSubData( type, offset, size, data );
        }
    }



    namespace Texturing
    {
        /////// GLOOM
		template<typename texel>
		inline void upload( Texture<texel>& texture )
		{
			DEBUG_ASSERT( texture.id>0 && textureBound() == texture.id );

			if( texture.videomem_invalidated )
			{
				int internalFormat, format, type;
				getFormat<texel>( internalFormat, format, type );

				glBindTexture( GL_TEXTURE_2D, texture.id );
				glTexImage2D( GL_TEXTURE_2D, 0, internalFormat, texture.size.width, texture.size.height, 
											 0, format, type, texture.array );

				texture.videomem_invalidated = false;
			}
		}

		template<typename texel>
		inline void download( Texture<texel>& texture )
		{
			DEBUG_ASSERT( texture.id>0 );//&& textureBound() == texture.id );

			if( texture.hostmem_invalidated )
			{
				int internalFormat, format, type;
				getFormat<texel>( internalFormat, format, type );

				glBindTexture( GL_TEXTURE_2D, texture.id );
				glGetTexImage( GL_TEXTURE_2D, 0, format, type, texture.array );

				texture.hostmem_invalidated = false;
			}
		}

        // TODO: Making 1 bind function out of the two? Is it possible?

		template<typename texel>
		void bind( int slot, Texture<texel>& texture )
		{	
            glEnable( GL_TEXTURE_2D );
			glActiveTexture( GL_TEXTURE0 + slot );
				
			if( texture.id==0 )
			{
				glGenTextures( 1, &texture.id );
				texture.deallocator = deallocate;
            }
			
            glBindTexture(GL_TEXTURE_2D, texture.id);
	
			if( texture.videomem_invalidated )
			{
				int internalFormat, format, type;
                getFormat<texel>( internalFormat, format, type );

                glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); 
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR); 
				glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, texture.size.x, texture.size.y, 
											 0, format, type, texture.array);
                texture.videomem_invalidated = false;
			}
		}

        // Bind and update a texture together with its given mipmaps (used for normalmap).
        template<typename texel>
		void bind( int slot, Texture<texel>* texture )
		{	
            glEnable( GL_TEXTURE_2D );
			glActiveTexture( GL_TEXTURE0 + slot );
				
			if( texture[0].id==0 )
			{
				glGenTextures( 1, &texture[0].id );
				texture[0].deallocator = deallocate;
            }
			
            glBindTexture(GL_TEXTURE_2D, texture[0].id);
	
			if( texture[0].videomem_invalidated )
			{
				int internalFormat, format, type;
                getFormat<texel>(internalFormat, format, type);

                glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); 
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR); 

                // Loading all mipmaps.
                // https://books.google.dk/books?id=E7eVY78Jo5YC&pg=PA167&lpg=PA167&dq=glTexImage2D+mipmapping&source=bl&ots=D7-QE78mql&sig=gSu2suthOqW5F3oIXYz8D5Xlj7Q&hl=en&sa=X&ved=0ahUKEwiEnb7tosLaAhWGJlAKHQ41DKAQ6AEIiAEwBw#v=onepage&q=glTexImage2D%20mipmapping&f=false                           
                int lods = log2(float(texture[0].size.x));
                for(int i = 0; i<lods; i++) {
                    int s = 1<<(lods-i);
				    glTexImage2D( GL_TEXTURE_2D, i, internalFormat, s, s, 0, format, type, texture[i].array );
                }
                texture[0].videomem_invalidated = false;
			}
		}
/////// GLOOM //

    }
} // GL



namespace GL
{
    // private
    void printSource(char* source, int length, int focusLine = -1)
	{
        ASSERT( length <= MAX_PROGRAM_SOURCE_LENGTH );

        char buffer[MAX_PROGRAM_SOURCE_LENGTH];
        strncpy(buffer, source, length);

        char* p = buffer;
        char lineText[256];
        for(int line = 2; *p; line++) 
		{
	        char* nIndex = strstr2(p, "\n");
	        if(nIndex == NULL) break;

            if((focusLine < 0) || (focusLine >= 0) && (focusLine-6 <= line) && (line <= focusLine+3)) 
            { 
	            int length = nIndex - p;			
	            strncpy(lineText, p, length);
	            lineText[length] = 0;

                if(line == focusLine) {
                    color(CMD_WHITE, CMD_RED); 
                } else {
                    color(CMD_LIGHTGRAY, 0); 
                }
	            printf("%i: %s\n", line, lineText);
            }

	        p = nIndex + 1;
        }
    }

    int CompileShader(char* source, int len, int type)
    {
        for(; *source>0 && *source<=' '; source++) len--;

        ASSERT(*source);
        ASSERT(len >= 0);

        GLint shaderId = glCreateShader( type );
        glShaderSource(shaderId, 1, (const char**)&source, &len );
        glCompileShader(shaderId);

        int successful = 0;
        glGetShaderiv( shaderId, GL_COMPILE_STATUS, &successful );

        if( !successful )
        {
            color(CMD_LIGHTRED,0); DEBUG_PRINT("\n\n");

            char message[512];
            glGetShaderInfoLog( shaderId, sizeof(message), NULL, message );
			printf("%s\n", message);
            //DEBUG_PRINT(message);

			// Parse to get the error line.
			// NOTE: Different drivers give different error messages
			// https://gamedev.stackexchange.com/questions/38685/is-this-a-reliable-method-of-parsing-glgetshaderinfolog
			int line = -1;
			if(line <= 0)
            {
				// NVIDIA error message
				char *p1 = 0, *p2 = 0;
				p1 = strchr(message, '(');
				if( p1>0 ) p2 = strchr(p1+1, ')');
				if( p1>0 && p2>0 ) line = strtol(p1+1, &p2, 10) + 1; // it was +0 
			}
			if(line <= 0)
            {
				// ATI/Intel error message
				char *p0 = 0, *p1 = 0, *p2 = 0;
				p0 = strchr(message, ':');
				if( p0>0 ) p1 = strchr(p0+1, ':');
				if( p1>0 ) p2 = strchr(p1+1, ':');
				if(p0>0 && p1>0 && p2>0) line = strtol(p1+1, &p2, 10) + 1;
			}
			ASSERT(line > 0);

            #if defined(DEBUG_SHOW_GLSL_SOURCE)
	            printSource( source, len, line ); 
                printf("\n");
            #endif

            CRITICAL("Shader compile error");
        }  	
        return shaderId;
    }

    void CheckProgram( int programId )
    {
        int successful;
        glGetProgramiv( programId, GL_LINK_STATUS, &successful);
        if( successful ) return;

        char message[512];
        glGetProgramInfoLog( programId, sizeof(message), 0, message );
        DEBUG_CRITICAL( message );
    }

    void AttachShaderToProgram(int programId, int shaderId)
    {
        glAttachShader(programId, shaderId);
    }

    void LinkProgram(int programId)
    {
	    glLinkProgram(programId);
        CheckProgram(programId);
    }

    int BuildShaders(const char* source, int* shaders)
    {
        DEBUG_PATH;

        const char* shaderTypes[] = { "VERTEX", "CONTROL", "EVALUATION", "GEOMETRY", "FRAGMENT", "COMPUTE" };
        const int   shaderTypeCodes[] = { GL_VERTEX_SHADER, GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER, GL_COMPUTE_SHADER };
        int shadersCount = sizeof(shaderTypeCodes)/sizeof(int);

        int shaderIndex = 0;
        for( int i=0; i<shadersCount; i++ )
        {
            char* shaderPosition = strstr2( (char*)source, shaderTypes[i]);
            if(shaderPosition == NULL) continue;
            
            shaderPosition += strlen(shaderTypes[i]);
            DEBUG_ASSERT(*shaderPosition == ':');
            shaderPosition++;

            // Find next shader (if any) to compute the current shader length in characters.
	        int shaderLenght = 0;
	        for( int j=i+1; j<shadersCount; j++ )
	        {
                char* nextShaderPosition = strstr2( (char*)source, shaderTypes[j]);
                if(nextShaderPosition != NULL)
                {
		            shaderLenght = nextShaderPosition - shaderPosition;
		            break;
                }
	        }

	        if( shaderLenght==0 )
	        {
		        shaderLenght = strlen(shaderPosition);
	        }
		
            color(CMD_YELLOW, 0); DEBUG_PRINT("%s ", shaderTypes[i]);

	        int shaderId = GL::CompileShader(shaderPosition, shaderLenght, shaderTypeCodes[i]);
            DEBUG_ASSERT(shaderId);

            shaders[shaderIndex++] = shaderId;
        }
        DEBUG_PRINT("\n");
        ASSERT(shaderIndex > 0);
        shaders[shaderIndex] = 0;
        
        return shaderIndex;
    }

    int GetUniformLocation( int programId, char* name )
    {
        int location = glGetUniformLocation( programId, name );

        if( location<0 )
        {	
	        char buffer[200];
	        sprintf( buffer, "Uniform '%s' not found.\n", name );
	        DEBUG_CRITICAL( buffer );
        }
        return location;
    }
    void SetUniform( int programId, char* name, int value )
    {
        glUniform1i( GetUniformLocation(programId, name), value );
    }
    void SetUniform( int programId, char* name, float value )
    {
        glUniform1f( GetUniformLocation(programId, name), value );
    }
    void SetUniform( int programId, char* name, vec2 vector )
    {
        glUniform2fv( GetUniformLocation(programId, name), 1, vector.array );
    }
    void SetUniform( int programId, char* name, vec3& vector )
    {
        glUniform3fv( GetUniformLocation(programId, name), 1, vector.array );
    }
    void SetUniform( int programId, char* name, vec4& vector )
    {
        glUniform4fv( GetUniformLocation(programId, name), 1, vector.array );
    }
    void SetUniform( int programId, char* name, mat4& matrix )
    {
        glUniformMatrix4fv( GetUniformLocation(programId, name), 1, 0, matrix.array );
    }
    void SetUniform( int programId, char* name, mat3& matrix )
    {
        glUniformMatrix3fv( GetUniformLocation(programId, name), 1, 0, matrix.array );
    }

    int GetCurrentProgramId()
    {
        GLint _currentProgram;
        glGetIntegerv(GL_CURRENT_PROGRAM, &_currentProgram);
        return _currentProgram;
    }

    void UseProgram(int programId)
    {
        glUseProgram(programId);
    }
}// GL




namespace GL
{
    namespace GLSL
    {
        void Check(int requiredVersion = 4.3)
	    {
            DEBUG_TRACE(GL::renderer());
            DEBUG_TRACE(GL::version());
            DEBUG_ASSERT( GL::version() >= requiredVersion );

	        GL::Texturing::available();
	        //GL::Texturing::textureNonPowerOfTwoAvailable();
	        //GL::secondaryColorAvailable();
	        GL::swapControlAvailable();
	        GL::VBO::available();
            GL::GLSL::available();
	        GL::GLSL::tessellatorAvailable();
	        GL::GLSL::computeAvailable();
	        GL::MRT::available();

            if( GL::swapControlAvailable() )
            {
                GL::swapControl(0); // TODO check
            }

            if(!Verbose) 
            {
                return;
            }
    
		    if(true) 
            {
                int maxTextureUnits;
			    glGetIntegerv(GL_MAX_TEXTURE_UNITS, &maxTextureUnits);
                color(CMD_DARKGRAY, 0);
                DEBUG_TRACE(maxTextureUnits);
		    }

		    if( GL::GLSL::available() ) 
            {
                int maxTextureImageUnits;
			    glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &maxTextureImageUnits);
                color(CMD_DARKGRAY, 0);
                DEBUG_TRACE(maxTextureImageUnits);
		    }

		    if(glEnableVertexAttribArray) 
            {
                int maxVertexAttributes;
			    glGetIntegerv( GL_MAX_VERTEX_ATTRIBS, &maxVertexAttributes );
                color(CMD_DARKGRAY, 0);
                DEBUG_TRACE(maxVertexAttributes);
		    }

	        if( GL::GLSL::tessellatorAvailable() ) 
            {
                int maxPatchVertices;
		        glGetIntegerv(GL_MAX_PATCH_VERTICES, &maxPatchVertices);
                color(CMD_DARKGRAY, 0);
                DEBUG_TRACE(maxPatchVertices);
	        }

		    if(true)
            {
			    int maxDrawBuffers;
			    glGetIntegerv(GL_MAX_DRAW_BUFFERS, &maxDrawBuffers);
                color(CMD_DARKGRAY, 0);
                DEBUG_TRACE(maxDrawBuffers);
		    }
        } 
    } // GLSL

    void Initialize(int minimumVersion = 4.3)
    {
        GLSL::Check(minimumVersion);

		glEnable(GL_TEXTURE_2D);
		glEnable(GL_LIGHTING);
		
		glEnable(GL_DEPTH_TEST);
		glClearDepth(1.0f); 
		glDepthFunc(GL_LEQUAL);
		
		glShadeModel( GL_SMOOTH );	
		glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

		glDrawBuffer( GL_BACK );

		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );

		bool drawBackfaces = false;
		(drawBackfaces ? glDisable : glEnable)( GL_CULL_FACE );

		//GL::Texturing::unbind();
            GL::Texturing::unbind(0);
            GL::Texturing::unbind(1);
            GL::Texturing::unbind(2);
            glDisable(GL_BLEND);

		glMatrixMode(GL_TEXTURE);
		glLoadIdentity();
		
		// Reset VA state
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_SECONDARY_COLOR_ARRAY);
		glDisableClientState(GL_INDEX_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_EDGE_FLAG_ARRAY);

		// Reset VBO state
        if( GL::VBO::available() ) 
        {
			GL::VBO::unbind(GL_ARRAY_BUFFER);
			GL::VBO::unbind(GL_ELEMENT_ARRAY_BUFFER);
		}

        GL::UseProgram(0);
	}

}