#pragma once

//
// Part of the VBO and Texture logic is here.
// But it should be split into separate objects: VertexBuffer, Texture, etc.
// C++ Objects for OpenGL (Gluon, Glam, Glamourus)
// 

#include <gl/link.h>
#include <gl/buffer.h>
#include <gl/Texture.h>
#include <gl/VertexBuffer.h>
#include <gl/IndexBuffer.h>

namespace GL
{ 
    namespace VBO
    {
        void deallocate( uint& id )
        {
	        if( id ) glDeleteBuffers(1, &id);
	        id = 0;
        }

        void create( Buffer& buffer )
        {
	        glGenBuffers( 1, &buffer.id );
	        buffer.deallocator = deallocate;
        }

        void bind( Buffer& buffer, int type )//=GL_ARRAY_BUFFER )
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

        void allocate( Buffer& buffer, int size, void* data = NULL, int type = GL_ARRAY_BUFFER, int usage = GL_STATIC_DRAW )
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

        void update( Buffer& buffer, int size, void* data = NULL, int type = GL_ARRAY_BUFFER, int usage = GL_STATIC_DRAW )
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

        void upload( Buffer& buffer, int offset, int size, void* data, int type = GL_ARRAY_BUFFER )
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





