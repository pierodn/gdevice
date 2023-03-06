#pragma once

#define GL_GLEXT_PROTOTYPES // this makes ARB definitions to be included
#ifdef _WIN32
	#include "GL/gl.h"
	#include "GL/glext.h"
#else
	#include <OpenGL/glext.h>
#endif

//#include "engines/graphic/material.h"
#include "__temp/Light.h"
#include "__temp/Texture.h"
#include "gpu/Program.h"
#include "__temp/VertexBuffer.h"
#include "gpu/IndexBuffer.h"

namespace 
{
    PROC getProc( char* name, char* alternative1, char* alternative2 )
    {
	    PROC address = NULL;
	    address = wglGetProcAddress( name );
	    if( !address ) address = wglGetProcAddress( alternative1 );
	    if( !address ) address = wglGetProcAddress( alternative2 );
	    if( !address ) DEBUG_PRINT(TAB32 ": MISSING\n", name );
	    return address;
    }

    PROC getProc( char* name, char* alternative )
    {
	    PROC address = NULL;
	    address = wglGetProcAddress( name );
	    if( !address ) address = wglGetProcAddress( alternative );
	    if( !address ) DEBUG_PRINT(TAB32 ": MISSING\n", name );
	    return address;
    }

    PROC getProc( char* name )
    {
	    char alternative[80];
	    sprintf( alternative, "%sARB", name );
	    return getProc( name, alternative );
    }
}

/// ============================
//  OpenGL
/// ============================

#if defined(_WIN32)
	extern PFNGLMULTITEXCOORD1FPROC		pglMultiTexCoord1f;
	extern PFNGLMULTITEXCOORD2FPROC		pglMultiTexCoord2f;
	extern PFNGLMULTITEXCOORD3FPROC		pglMultiTexCoord3f;
	extern PFNGLMULTITEXCOORD4FPROC		pglMultiTexCoord4f;
	extern PFNGLACTIVETEXTUREPROC		pglActiveTexture;
	extern PFNGLCLIENTACTIVETEXTUREPROC	pglClientActiveTexture;
	#define glMultiTexCoord1f			pglMultiTexCoord1f
	#define glMultiTexCoord2f			pglMultiTexCoord2f
	#define glMultiTexCoord3f			pglMultiTexCoord3f
	#define glMultiTexCoord4f			pglMultiTexCoord4f
	#define glActiveTexture				pglActiveTexture
	#define glClientActiveTexture		pglClientActiveTexture
	PFNGLMULTITEXCOORD1FPROC			pglMultiTexCoord1f		= 0;
	PFNGLMULTITEXCOORD2FPROC			pglMultiTexCoord2f		= 0;
	PFNGLMULTITEXCOORD3FPROC			pglMultiTexCoord3f		= 0;
	PFNGLMULTITEXCOORD4FPROC			pglMultiTexCoord4f		= 0;
	PFNGLACTIVETEXTUREPROC				pglActiveTexture		= 0;
	PFNGLCLIENTACTIVETEXTUREPROC		pglClientActiveTexture	= 0;	
#endif

namespace GL
{ 
	inline const char* vendor()
	{
		return (const char *)glGetString(GL_VENDOR);
	}

	double version()
	{
		return atof((const char *)glGetString(GL_VERSION));
	}

	bool isExtensionSupported( const char *extension )
	{
		char* s = (char*) glGetString(GL_EXTENSIONS);
		while( (s=strstr2(s,extension)) ) {
			s += strlen(extension);
			if( *s==' ' || *s == '\0' ) break;
		}

		//DEBUG_PRINT( TAB32 ": %s\n", extension, s ? "OK" : "MISSING" );
		return (bool)s;
	}

	inline char* error()
	{
		GLenum error = glGetError();
		return 
            error == GL_NO_ERROR                        ? "No error" :
			error == GL_INVALID_ENUM		            ? "Invalid enum" :
			error == GL_INVALID_VALUE		            ? "Invalid value" :
			error == GL_INVALID_OPERATION               ? "Invalid operation" :
			error == GL_STACK_OVERFLOW                  ? "Stack overflow" :
			error == GL_STACK_UNDERFLOW                 ? "Stack underflow" :
			error == GL_OUT_OF_MEMORY                   ? "Out of memory" : 
            error == GL_INVALID_FRAMEBUFFER_OPERATION   ? "Invalid framebuffer operation" :
                                                          "Unknown error";
	}

    /// ============================
    //  Materials
    /// ============================
/*
	namespace Materials
	{
		void bind( Material& material )
		{
			glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   material.ambient.array );
			glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   material.diffuse.array );
			glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  material.emission.array );
			glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  material.specular.array );
			glMaterialf(  GL_FRONT_AND_BACK, GL_SHININESS, material.shininess );
		}
	};*/

    /// ============================
    //  Lighting
    /// ============================

	namespace Lighting
	{
		void bind( Light& light )
		{
			int id = GL_LIGHT0 + light.id;

			glEnable( id );
		
			if( light.positional )
			{
				glLightfv( id, GL_SPOT_DIRECTION, light.spot_direction.array );
				glLighti(  id, GL_SPOT_EXPONENT,  light.spot_exponent );
				glLighti(  id, GL_SPOT_CUTOFF,    light.spot_cutoff );
				glLightf(  id, GL_CONSTANT_ATTENUATION,  light.constant_attenuation );
				glLightf(  id, GL_LINEAR_ATTENUATION,    light.linear_attenuation );
				glLightf(  id, GL_QUADRATIC_ATTENUATION, light.quadratic_attenuation );
			}

			glLightfv( id, GL_POSITION, light.position.array );
			glLightfv( id, GL_AMBIENT,  light.ambient.array  );
			glLightfv( id, GL_DIFFUSE,  light.diffuse.array  );
			glLightfv( id, GL_SPECULAR, light.specular.array );
		}

		void unbind( Light& light )
		{
			glDisable( GL_LIGHT0 + light.id );
		}
	};

    /// ============================
    //  Texturing
    /// ============================

	namespace Texturing
	{
		bool available()
		{
			static int supported = -1;

			if( supported==-1 )
			{
				supported = isExtensionSupported("GL_ARB_multitexture");
						 
				#if defined(_WIN32)
				supported = supported
					&& (glMultiTexCoord1f		= (PFNGLMULTITEXCOORD1FPROC)		getProc("glMultiTexCoord1f"))
					&& (glMultiTexCoord2f		= (PFNGLMULTITEXCOORD2FPROC)		getProc("glMultiTexCoord2f"))
					&& (glMultiTexCoord3f		= (PFNGLMULTITEXCOORD3FPROC)		getProc("glMultiTexCoord3f"))
					&& (glMultiTexCoord4f		= (PFNGLMULTITEXCOORD4FPROC)		getProc("glMultiTexCoord4f"))
					&& (glActiveTexture			= (PFNGLACTIVETEXTUREPROC)			getProc("glActiveTexture"))
					&& (glClientActiveTexture	= (PFNGLCLIENTACTIVETEXTUREPROC)	getProc("glClientActiveTexture"))
				;
				#endif

				supported = supported 
					&& isExtensionSupported("GL_ARB_texture_env_combine")
					&& isExtensionSupported("GL_ARB_texture_env_dot3");

				//printf( "%-8s: %s\n", "FIXED", supported ? "OK" : "NOT AVAILABLE" );
			}

			return supported;
		}

		bool textureNonPowerOfTwoAvailable()
		{
			static int supported = -1;

			if( supported==-1 )
			{
				supported = isExtensionSupported("GL_ARB_texture_non_power_of_two" );
			}

			return supported;
		}

		template<typename texel>
		inline void getFormat( int& internalFormat, int& format, int& type )
		{
			switch( sizeof(texel) )
			{
			case sizeof(byte)	:	internalFormat = GL_LUMINANCE;	format = GL_LUMINANCE;	type = GL_UNSIGNED_BYTE;	break;
			case sizeof(byte)*3	:	internalFormat = GL_RGB;		format = GL_RGB;		type = GL_UNSIGNED_BYTE;	break;
			case sizeof(rgba)	:	internalFormat = GL_RGBA;		format = GL_RGBA;		type = GL_UNSIGNED_BYTE;	break;
			case sizeof(vec3)	:	internalFormat = GL_RGB32F;		format = GL_RGB;		type = GL_FLOAT;			break; 
			case sizeof(vec4)	:	internalFormat = GL_RGBA32F;	format = GL_RGBA;		type = GL_FLOAT;			break; 
			default				:	DEBUG_CRITICAL("Unknown texture format");
			}
		}

		inline int textureBound()
		{
			int current;
			glGetIntegerv(GL_TEXTURE_BINDING_2D, &current);
            return current;
		}

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
				glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, texture.size.width, texture.size.height, 
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
                int lods = log2(float(texture[0].size.width));
                for(int i = 0; i<lods; i++) {
                    int s = 1<<(lods-i);
				    glTexImage2D( GL_TEXTURE_2D, i, internalFormat, s, s, 0, format, type, texture[i].array );
                }
                texture[0].videomem_invalidated = false;
			}
		}


		void deallocate( uint& id )
		{
			glDeleteTextures(1, &id );
		}

		void bind()
		{
			vec4 constant;
			constant.rgb = vec3(0);	// color biasing 
			constant.a = 1.0;		// 4th detail amplitude, almost all GL implementations doesnt allow > 1.0
			glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, constant.array );
			
			glActiveTexture(GL_TEXTURE0);
	//		glTexEnvf(GL_TEXTURE_FILTER_CONTROL, GL_TEXTURE_LOD_BIAS, (distance*distance-1.0)*0.002 );
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);

			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB,  GL_DOT3_RGB);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_RGB,  GL_TEXTURE0);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE1_RGB,  GL_TEXTURE1);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_RGB, GL_SRC_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND1_RGB, GL_SRC_COLOR);

			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA,  GL_MODULATE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA,  GL_TEXTURE0);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE1_ALPHA,  GL_TEXTURE1);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND1_ALPHA, GL_SRC_ALPHA); 

			glActiveTexture(GL_TEXTURE1);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);

			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB,  GL_INTERPOLATE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_RGB,  GL_TEXTURE2);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE1_RGB,  GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE2_RGB,  GL_PRIMARY_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_RGB, GL_SRC_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND1_RGB, GL_SRC_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND2_RGB, GL_SRC_ALPHA);

			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA, GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_ALPHA, GL_SRC_ALPHA);

			glActiveTexture(GL_TEXTURE2);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);

			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB,  GL_INTERPOLATE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_RGB,  GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE1_RGB,  GL_CONSTANT);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE2_RGB,  GL_PRIMARY_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_RGB, GL_SRC_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND1_RGB, GL_SRC_COLOR);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND2_RGB, GL_SRC_COLOR);

			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_MODULATE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA, GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE1_ALPHA, GL_CONSTANT);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_ALPHA, GL_SRC_ALPHA);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND1_ALPHA, GL_SRC_ALPHA);

			glEnable(GL_BLEND);
			glBlendFunc( GL_ONE_MINUS_SRC_ALPHA, GL_ZERO );
		}

		void unbind( int slot )
		{
			glActiveTexture( GL_TEXTURE0 + slot );
			glDisable(GL_TEXTURE_2D);
		}

		void unbind()
		{
			unbind(0);
			unbind(1);
			unbind(2);
			glDisable(GL_BLEND);
		}
	};
};

/// ============================
//  GL_EXT_secondary_color
/// ============================

#if defined(_WIN32)
	extern	PFNGLSECONDARYCOLORPOINTERPROC	pglSecondaryColorPointer;
	#define	glSecondaryColorPointer			pglSecondaryColorPointer
	PFNGLSECONDARYCOLORPOINTERPROC			pglSecondaryColorPointer = 0;
#endif

namespace GL
{ 
	bool secondaryColorAvailable()
	{
		static int supported = -1;

		if( supported==-1 )
		{
			supported = isExtensionSupported("GL_EXT_secondary_color" );
		
			#if defined(_WIN32)
			supported = supported
				&& (glSecondaryColorPointer	= (PFNGLSECONDARYCOLORPOINTERPROC) getProc("glSecondaryColorPointer", "glSecondaryColorPointerEXT") )
			;
			#endif

			//printf( "%-8s: %s\n", "SECONDARY COLOR", supported ? "OK" : "NOT AVAILABLE" );
		}

		return supported;
	}
};


/// ============================
//  WGL_EXT_swap_control
/// ============================

#if defined(_WIN32)
	typedef void (APIENTRY *PFNWGLEXTSWAPCONTROLPROC)(int);
	typedef int (*PFNWGLEXTGETSWAPINTERVALPROC)(void);
	PFNWGLEXTSWAPCONTROLPROC		wglSwapIntervalEXT = 0;
	PFNWGLEXTGETSWAPINTERVALPROC	wglGetSwapIntervalEXT = 0;
#endif

namespace GL
{ 
	bool swapControlAvailable()
	{
		static int supported = -1;
		if( supported==-1 ) {
			supported = isExtensionSupported( "WGL_EXT_swap_control" );
			#if defined(_WIN32)
			supported = supported
				&& (wglSwapIntervalEXT		= (PFNWGLEXTSWAPCONTROLPROC)	 getProc("wglSwapInterval",	"wglSwapIntervalEXT") )
				&& (wglGetSwapIntervalEXT	= (PFNWGLEXTGETSWAPINTERVALPROC) getProc("wglGetSwapInterval", "wglGetSwapIntervalEXT") )
			;
			#endif
		}
		return supported;
	}

    void swapControl(int swapInterval)
    {
        wglSwapIntervalEXT(swapInterval);
    }
};

/// ============================================ 
//  Multiple Render Targets
/// ============================================

#if defined(_WIN32)
	extern PFNGLBINDFRAMEBUFFERPROC			pglBindFramebuffer;
	extern PFNGLBINDRENDERBUFFERPROC		pglBindRenderbuffer;
	//extern PFNGLBLITFRAMEBUFFERPROC			pglBlitFramebuffer;
	extern PFNGLCHECKFRAMEBUFFERSTATUSPROC	pglCheckFramebufferStatus;
	extern PFNGLDELETEFRAMEBUFFERSPROC		pglDeleteFramebuffers;
	extern PFNGLDELETERENDERBUFFERSPROC		pglDeleteRenderbuffers;
	extern PFNGLFRAMEBUFFERRENDERBUFFERPROC	pglFramebufferRenderbuffer;
	extern PFNGLFRAMEBUFFERTEXTUREPROC		pglFramebufferTexture;
	extern PFNGLFRAMEBUFFERTEXTURE1DPROC	pglFramebufferTexture1D;
	extern PFNGLFRAMEBUFFERTEXTURE2DPROC	pglFramebufferTexture2D;
	extern PFNGLFRAMEBUFFERTEXTURE3DPROC	pglFramebufferTexture3D;
	//extern PFNGLFRAMEBUFFERTEXTURELAYERPROC	pglFramebufferTextureLayer;
	extern PFNGLGENFRAMEBUFFERSPROC			pglGenFramebuffers;
	extern PFNGLGENRENDERBUFFERSPROC		pglGenRenderbuffers;
	extern PFNGLGENERATEMIPMAPPROC			pglGenerateMipmap;
	extern PFNGLGETFRAMEBUFFERATTACHMENTPARAMETERIVPROC	pglGetFramebufferAttachmentParameteriv;
	extern PFNGLGETRENDERBUFFERPARAMETERIVPROC	pglGetRenderbufferParameteriv;
	extern PFNGLISFRAMEBUFFERPROC			pglIsFramebuffer;
	extern PFNGLISRENDERBUFFERPROC			pglIsRenderbuffer;
	extern PFNGLRENDERBUFFERSTORAGEPROC		pglRenderbufferStorage;
	//extern PFNGLRENDERBUFFERSTORAGEMULTISAMPLEPROC			pglRenderbufferStorageMultisample;

	#define glBindFramebuffer				pglBindFramebuffer
	#define glBindRenderbuffer				pglBindRenderbuffer
	//#define glBlitFramebuffer				pglBlitFramebuffer
	#define glCheckFramebufferStatus		pglCheckFramebufferStatus
	#define glDeleteFramebuffers			pglDeleteFramebuffers
	#define glDeleteRenderbuffers			pglDeleteRenderbuffers
	#define glFramebufferRenderbuffer		pglFramebufferRenderbuffer
	#define glFramebufferTexture			pglFramebufferTexture
	#define glFramebufferTexture1D			pglFramebufferTexture1D
	#define glFramebufferTexture2D			pglFramebufferTexture2D
	#define glFramebufferTexture3D			pglFramebufferTexture3D
	//#define glFramebufferTextureLayer			pglFramebufferTextureLayer
	#define glGenFramebuffers				pglGenFramebuffers
	#define glGenRenderbuffers				pglGenRenderbuffers
	#define glGenerateMipmap				pglGenerateMipmap
	#define glGetFramebufferAttachmentParameteriv	pglGetFramebufferAttachmentParameteriv
	#define glGetRenderbufferParameteriv	pglGetRenderbufferParameteriv
	#define glIsFramebuffer					pglIsFramebuffer
	#define glIsRenderbuffer				pglIsRenderbuffer 
	#define glRenderbufferStorage			pglRenderbufferStorage 
	//#define glRenderbufferStorageMultisample	pglRenderbufferStorageMultisample

	PFNGLBINDFRAMEBUFFERPROC				glBindFramebuffer = 0;
	PFNGLBINDRENDERBUFFERPROC				glBindRenderbuffer = 0;
	//PFNGLBLITFRAMEBUFFERPROC				glBlitFramebuffer = 0;
	PFNGLCHECKFRAMEBUFFERSTATUSPROC			glCheckFramebufferStatus = 0;
	PFNGLDELETEFRAMEBUFFERSPROC				glDeleteFramebuffers = 0;
	PFNGLDELETERENDERBUFFERSPROC			glDeleteRenderbuffers = 0;
	PFNGLFRAMEBUFFERRENDERBUFFERPROC		glFramebufferRenderbuffer = 0;
	PFNGLFRAMEBUFFERTEXTUREPROC				glFramebufferTexture = 0;
	PFNGLFRAMEBUFFERTEXTURE1DPROC			glFramebufferTexture1D = 0;
	PFNGLFRAMEBUFFERTEXTURE2DPROC			glFramebufferTexture2D = 0;
	PFNGLFRAMEBUFFERTEXTURE3DPROC			glFramebufferTexture3D = 0;
	//PFNGLFRAMEBUFFERTEXTURELAYERPROC			glFramebufferTextureLayer = 0;
	PFNGLGENFRAMEBUFFERSPROC				glGenFramebuffers = 0;
	PFNGLGENRENDERBUFFERSPROC				glGenRenderbuffers = 0;
	PFNGLGENERATEMIPMAPPROC					glGenerateMipmap = 0;
	PFNGLGETFRAMEBUFFERATTACHMENTPARAMETERIVPROC	glGetFramebufferAttachmentParameteriv = 0;
	PFNGLGETRENDERBUFFERPARAMETERIVPROC		glGetRenderbufferParameteriv = 0;
	PFNGLISFRAMEBUFFERPROC					glIsFramebuffer = 0;
	PFNGLISRENDERBUFFERPROC					glIsRenderbuffer = 0;
	PFNGLRENDERBUFFERSTORAGEPROC			glRenderbufferStorage = 0;
	//PFNGLRENDERBUFFERSTORAGEMULTISAMPLEPROC		glRenderbufferStorageMultisample = 0;
#endif

namespace GL
{
    namespace MRT
    {
        bool available()
        {
	        static int supported = -1;

	        if( supported==-1 )
	        {
		        supported = GL::isExtensionSupported("GL_EXT_framebuffer_object") 
			             || GL::isExtensionSupported("GL_ARB_framebuffer_object");

		        #if defined(_WIN32)
		        supported = supported
			        && (glBindFramebuffer			= (PFNGLBINDFRAMEBUFFERPROC)			getProc("glBindFramebuffer",		"glBindFramebufferEXT" ))
			        && (glBindRenderbuffer			= (PFNGLBINDRENDERBUFFERPROC)			getProc("glBindRenderbuffer",		"glBindRenderbufferEXT" ))
			        //&& (glBlitFramebuffer			= (PFNGLBLITFRAMEBUFFERPROC)			getProc("glBlitFramebuffer",		"glBlitFramebufferEXT" ))
			        && (glCheckFramebufferStatus	= (PFNGLCHECKFRAMEBUFFERSTATUSPROC)		getProc("glCheckFramebufferStatus",	"glCheckFramebufferStatusEXT" ))
			        && (glDeleteFramebuffers		= (PFNGLDELETEFRAMEBUFFERSPROC)			getProc("glDeleteFramebuffers",		"glDeleteFramebuffersEXT" ))
			        && (glDeleteRenderbuffers		= (PFNGLDELETERENDERBUFFERSPROC)		getProc("glDeleteRenderbuffers",	"glDeleteRenderbuffersEXT" ))
			        && (glFramebufferRenderbuffer	= (PFNGLFRAMEBUFFERRENDERBUFFERPROC)	getProc("glFramebufferRenderbuffer","glFramebufferRenderbufferEXT" ))
			        && (glFramebufferTexture		= (PFNGLFRAMEBUFFERTEXTUREPROC)			getProc("glFramebufferTexture" )) // 3.2 only
			        && (glFramebufferTexture1D		= (PFNGLFRAMEBUFFERTEXTURE1DPROC)		getProc("glFramebufferTexture1D",	"glFramebufferTexture1DEXT" ))
			        && (glFramebufferTexture2D		= (PFNGLFRAMEBUFFERTEXTURE2DPROC)		getProc("glFramebufferTexture2D",	"glFramebufferTexture2DEXT" ))
			        && (glFramebufferTexture3D		= (PFNGLFRAMEBUFFERTEXTURE3DPROC)		getProc("glFramebufferTexture3D",	"glFramebufferTexture3DEXT" ))
			        //&& (glFramebufferTextureLayer	= (PFNGLFRAMEBUFFERTEXTURELAYERPROC)	getProc("glFramebufferTextureLayer","glFramebufferTextureLayerEXT" ))
			        && (glGenFramebuffers			= (PFNGLGENFRAMEBUFFERSPROC)			getProc("glGenFramebuffers",		"glGenFramebuffersEXT" ))
			        && (glGenRenderbuffers			= (PFNGLGENRENDERBUFFERSPROC)			getProc("glGenRenderbuffers",		"glGenRenderbuffersEXT" ))
			        && (glGenerateMipmap			= (PFNGLGENERATEMIPMAPPROC)				getProc("glGenerateMipmap",			"glGenerateMipmapEXT" ))
			        && (glGetFramebufferAttachmentParameteriv = (PFNGLGETFRAMEBUFFERATTACHMENTPARAMETERIVPROC)	getProc("glGetFramebufferAttachmentParameteriv",		"glGetFramebufferAttachmentParameterivEXT" ))
			        && (glGetRenderbufferParameteriv = (PFNGLGETRENDERBUFFERPARAMETERIVPROC) getProc("glGetRenderbufferParameteriv", "glGetRenderbufferParameterivEXT" ))
			        && (glIsFramebuffer				= (PFNGLISFRAMEBUFFERPROC)				getProc("glIsFramebuffer",			"glIsFramebufferEXT" ))
			        && (glIsRenderbuffer			= (PFNGLISRENDERBUFFERPROC)				getProc("glIsRenderbuffer",			"glIsRenderbufferEXT" ))
			        && (glRenderbufferStorage		= (PFNGLRENDERBUFFERSTORAGEPROC)		getProc("glRenderbufferStorage",	"glRenderbufferStorageEXT" ))
			        //&& (glRenderbufferStorageMultisample	= (PFNGLRENDERBUFFERSTORAGEMULTISAMPLEPROC)	getProc("glRenderbufferStorageMultisample",		"glRenderbufferStorageMultisampleEXT" ))
		        ;
		        #endif

		        //printf( "%-8s: %s\n", "MRT", supported ? "OK" : "NOT AVAILABLE" );
	        }

	        return supported;
        }

        void deallocate( uint& id )
        {
	        glDeleteFramebuffers( 1, &id );
        }
/*
        void bindTarget( VertexBuffer& vbo )
        {
	        if( vbo.id == 0 ) // used as FBO !!
	        {
		        // force first upload
		        GL::Texturing::bind(0, vbo.quartets );
		        GL::Texturing::bind(1, vbo.normals_ff );
		        GL::Texturing::bind(2, vbo.colors );
		        GL::Texturing::bind(3, vbo.mixmaps );
		        GL::Texturing::unbind(0);
		        GL::Texturing::unbind(1);
		        GL::Texturing::unbind(2);
		        GL::Texturing::unbind(3);
		        glBindTexture( GL_TEXTURE_2D, 0 );

		        glGenFramebuffers(1, &vbo.id);
		        vbo.deallocator = deallocate;

		        glBindFramebuffer(GL_FRAMEBUFFER, vbo.id);
		        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, vbo.quartets.id, 0);
		        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, vbo.normals_ff.id,  0);
		        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, vbo.colors.id,	0);
		        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, vbo.mixmaps.id,	0);
#if 0
		        if( glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE )
		        {
			        exit_on_error( "Framebuffer not complete!" );
		        }
#endif

		        int status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		        if( status != GL_FRAMEBUFFER_COMPLETE )
		        {
			        printf("status = 0x%x\n", status);

			        if( status == 0 )
			        {
				        printf("%s\n", GL::error() );
			        }

			        EXIT( 
				        status==GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT ?			"Incomplete attachment" :  
				        status==GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT ?	"Missing attachment" : 
				        status==GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT ?		"Incomplete dimensions" : 
				        status==GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT ?			"Incomplete formats" : 
				        status==GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER ?			"Incomplete draw buffer" : 
				        status==GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER ?			"Incomplete read buffer" : 
				        status==GL_FRAMEBUFFER_UNSUPPORTED ?					"Unsupported" : 
																		        "Unknown"
			        );
        			
		        }

	        }

	        glBindFramebuffer(GL_FRAMEBUFFER, vbo.id);

	        GLenum buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };
	        glDrawBuffers(4, buffers);

	        glPushAttrib(GL_VIEWPORT_BIT);
	        glViewport( 0,0, vbo.quartets.size.width, vbo.quartets.size.height );
        }
*/
        void unbindTarget()
        {
	        glPopAttrib();
	        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }
    }
}

/// ============================================ 
//  Render To Texture
/// ============================================

#if defined(_WIN32)
	extern PFNGLBEGINTRANSFORMFEEDBACKPROC		pglBeginTransformFeedback;
	extern PFNGLENDTRANSFORMFEEDBACKPROC		pglEndTransformFeedback;
	//extern PFNGLBINDBUFFERBASEPROC				pglBindBufferBase;
	//extern  PFNGLBINDBUFFEROFFSETPROC			pglBindBufferOffset;
	//extern PFNGLBINDBUFFERRANGEPROC				pglBindBufferRange;
	extern PFNGLTRANSFORMFEEDBACKVARYINGSPROC	pglTransformFeedbackVaryings;
	extern PFNGLGETTRANSFORMFEEDBACKVARYINGPROC	pglGetTransformFeedbackVarying;
	#define glBeginTransformFeedback			pglBeginTransformFeedback
	#define glEndTransformFeedback				pglEndTransformFeedback
	//#define glBindBufferBase					pglBindBufferBase 
	//#define glBindBufferOffset				pglBindBufferOffset 
	//#define glBindBufferRange					pglBindBufferRange 
	#define glTransformFeedbackVaryings			pglTransformFeedbackVaryings
	#define glGetTransformFeedbackVarying		pglGetTransformFeedbackVarying 
	PFNGLBEGINTRANSFORMFEEDBACKPROC				pglBeginTransformFeedback = 0;
	PFNGLENDTRANSFORMFEEDBACKPROC				pglEndTransformFeedback = 0;
	//PFNGLBINDBUFFERBASEPROC					pglBindBufferBase = 0;
	//PFNGLBINDBUFFEROFFSETPROC					pglBindBufferOffset = 0;
	//PFNGLBINDBUFFERRANGEPROC					pglBindBufferRange = 0;
	PFNGLTRANSFORMFEEDBACKVARYINGSPROC			pglTransformFeedbackVaryings = 0;
	PFNGLGETTRANSFORMFEEDBACKVARYINGPROC		pglGetTransformFeedbackVarying = 0;
#endif

namespace GL 
{
    namespace RTT
    {
        bool available()
        {
	        static int supported = -1;

	        if( supported==-1 )
	        {
		        supported = GL::isExtensionSupported("GL_EXT_transform_feedback") 
			             || GL::isExtensionSupported("GL_NV_transform_feedback");	// uguale a EXT?	
				         //|| isExtensionSupported("GL_ARB_transform_feedback2"); // addictional functionalities

		        #if defined(_WIN32)
		        supported = supported
			        && (glBeginTransformFeedback		= (PFNGLBEGINTRANSFORMFEEDBACKPROC)		getProc("glBeginTransformFeedback",		"glBeginTransformFeedbackEXT",		"glBeginTransformFeedbackNV"))
			        && (glEndTransformFeedback			= (PFNGLENDTRANSFORMFEEDBACKPROC)		getProc("glEndTransformFeedback",		"glEndTransformFeedbackEXT",		"glEndTransformFeedbackNV"))
			        //&& (glBindBufferBase				= (PFNGLBINDBUFFERBASEPROC)				getProc("glBindBufferBase",				"glBindBufferBaseEXT",				"glBindBufferBaseNV"))
			        //&& (glBindBufferOffset			= (PFNGLBINDBUFFEROFFSETPROC)			getProc("glBindBufferOffsetEXT"))
			        //&& (glBindBufferRange				= (PFNGLBINDBUFFERRANGEPROC)			getProc("glBindBufferRange",			"glBindBufferRangeEXT",				"glBindBufferRangeNV"))
			        && (glTransformFeedbackVaryings		= (PFNGLTRANSFORMFEEDBACKVARYINGSPROC)	getProc("glTransformFeedbackVaryings",	"glTransformFeedbackVaryingsEXT",	"glTransformFeedbackVaryingsNV"))
			        && (glGetTransformFeedbackVarying	= (PFNGLGETTRANSFORMFEEDBACKVARYINGPROC)getProc("glGetTransformFeedbackVarying","glGetTransformFeedbackVaryingEXT",	"glGetTransformFeedbackVaryingNV"))
		        ;
		        #endif

		        // TODO or EXT_framebuffer_object + ARB_pixel_buffer_object
		        // TODO or EXT_framebuffer_object + glCopyPixel
		        // TODO or write_to_backbuffer    + glCopyPixel

		        //printf( "%-8s: %s\n", "RTT", supported ? "OK" : "NOT AVAILABLE" );
	        }

	        return supported;
        }
    }
}

/// ============================================ 
//  Vertex Buffer Objects
/// ============================================

#if defined(_WIN32)
	extern PFNGLGENBUFFERSPROC				pglGenBuffers; 
	extern PFNGLBINDBUFFERPROC				pglBindBuffer; 
	extern PFNGLBUFFERDATAPROC				pglBufferData; 
	extern PFNGLBUFFERSUBDATAPROC			pglBufferSubData; 
	extern PFNGLDELETEBUFFERSPROC			pglDeleteBuffers; 
	extern PFNGLGETBUFFERPARAMETERIVPROC	pglGetBufferParameteriv; 
	extern PFNGLMAPBUFFERPROC				pglMapBuffer; 
	extern PFNGLUNMAPBUFFERPROC				pglUnmapBuffer; 
	#define glGenBuffers					pglGenBuffers
	#define glBindBuffer					pglBindBuffer
	#define glBufferData					pglBufferData
	#define glBufferSubData					pglBufferSubData
	#define glDeleteBuffers					pglDeleteBuffers
	#define glGetBufferParameteriv			pglGetBufferParameteriv
	#define glMapBuffer						pglMapBuffer
	#define glUnmapBuffer					pglUnmapBuffer
	PFNGLGENBUFFERSPROC						pglGenBuffers			= 0;
	PFNGLBINDBUFFERPROC						pglBindBuffer			= 0;
	PFNGLBUFFERDATAPROC						pglBufferData			= 0;
	PFNGLBUFFERSUBDATAPROC					pglBufferSubData		= 0;
	PFNGLDELETEBUFFERSPROC					pglDeleteBuffers		= 0;
	PFNGLGETBUFFERPARAMETERIVPROC			pglGetBufferParameteriv	= 0;
	PFNGLMAPBUFFERPROC						pglMapBuffer			= 0;
	PFNGLUNMAPBUFFERPROC					pglUnmapBuffer			= 0;
#endif

namespace GL
{
    namespace VBO
    {
        bool available()
        {
	        static int supported = -1;

	        if( supported==-1 )
	        {
		        // ?
		        supported = GL::isExtensionSupported("GL_ARB_vertex_buffer_object");

		        #if defined(_WIN32)
		        supported = supported
			        && (glGenBuffers			= (PFNGLGENBUFFERSPROC)				getProc("glGenBuffers"))
			        && (glBindBuffer			= (PFNGLBINDBUFFERPROC)				getProc("glBindBuffer"))
			        && (glBufferData			= (PFNGLBUFFERDATAPROC)				getProc("glBufferData"))
			        && (glBufferSubData			= (PFNGLBUFFERSUBDATAPROC)			getProc("glBufferSubData"))
			        && (glDeleteBuffers			= (PFNGLDELETEBUFFERSPROC)			getProc("glDeleteBuffers"))
			        && (glGetBufferParameteriv	= (PFNGLGETBUFFERPARAMETERIVPROC)	getProc("glGetBufferParameteriv"))
			        && (glMapBuffer				= (PFNGLMAPBUFFERPROC)				getProc("glMapBuffer"))
			        && (glUnmapBuffer			= (PFNGLUNMAPBUFFERPROC)			getProc("glUnmapBuffer"));
		        #endif

		        //printf( "%-8s: %s\n", "VBO", supported ? "OK" : "NOT AVAILABLE" );
	        }

	        return supported;
        }

        void deallocate( uint& id )
        {
	        if( id ) glDeleteBuffers(1, &id);
	        id = 0;
        }

        void create( Cacheable& buffer )
        {
            //PROFILE_ONCE;
            DEBUG_CHECKPOINT_ONCE(create);
	        glGenBuffers( 1, &buffer.id );
            DEBUG_TRACE_ONCE(buffer.id);
	        buffer.deallocator = deallocate;
            DEBUG_TRACE_ONCE(deallocate);
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
            DEBUG_TRACE_ONCE(buffer.id);
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
            DEBUG_TRACE_ONCE(buffer.videomem_invalidated);
	        if( buffer.videomem_invalidated )
	        {
                DEBUG_TRACE_ONCE(buffer.id);
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
                if( !vbo.id ) {
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
}

/// ============================================ 
//  OpenGL Shading Language
/// ============================================

#define MAX_PROGRAM_SOURCE_LENGTH	262144
#define MAX_SHADER_SOURCE_LENGTH	65536

#if defined(_WIN32)
	extern PFNGLCREATESHADERPROC			pglCreateShader;
	extern PFNGLSHADERSOURCEPROC			pglShaderSource;
	extern PFNGLCOMPILESHADERPROC			pglCompileShader;
	extern PFNGLCREATEPROGRAMPROC			pglCreateProgram;
	extern PFNGLATTACHSHADERPROC			pglAttachShader;
	extern PFNGLLINKPROGRAMPROC				pglLinkProgram;
	extern PFNGLUSEPROGRAMPROC				pglUseProgram;
	extern PFNGLDELETESHADERPROC			pglDeleteShader;
	extern PFNGLDELETEPROGRAMPROC			pglDeleteProgram;
	extern PFNGLGETSHADERIVPROC				pglGetShaderiv;
	extern PFNGLGETSHADERINFOLOGPROC		pglGetShaderInfoLog;
	extern PFNGLGETPROGRAMIVPROC			pglGetProgramiv;
	extern PFNGLGETPROGRAMINFOLOGPROC		pglGetProgramInfoLog;
	extern PFNGLGETUNIFORMLOCATIONPROC		pglGetUniformLocation;
	extern PFNGLUNIFORM1FPROC				pglUniform1f;
	extern PFNGLUNIFORM2FPROC				pglUniform2f;
	extern PFNGLUNIFORM3FPROC				pglUniform3f;
	extern PFNGLUNIFORM4FPROC				pglUniform4f;
	extern PFNGLUNIFORM1FVPROC				pglUniform1fv;
	extern PFNGLUNIFORM2FVPROC				pglUniform2fv;
	extern PFNGLUNIFORM3FVPROC				pglUniform3fv;
	extern PFNGLUNIFORM4FVPROC				pglUniform4fv;
	extern PFNGLUNIFORM1IPROC				pglUniform1i;
	extern PFNGLUNIFORMMATRIX2FVPROC		pglUniformMatrix2fv;
	extern PFNGLUNIFORMMATRIX3FVPROC		pglUniformMatrix3fv;
	extern PFNGLUNIFORMMATRIX4FVPROC		pglUniformMatrix4fv;
	extern PFNGLVERTEXATTRIBPOINTERPROC		pglVertexAttribPointer;
	extern PFNGLENABLEVERTEXATTRIBARRAYPROC	pglEnableVertexAttribArray;
	extern PFNGLDISABLEVERTEXATTRIBARRAYPROC	pglDisableVertexAttribArray;
	extern PFNGLBINDATTRIBLOCATIONPROC		pglBindAttribLocation;
	extern PFNGLGETATTRIBLOCATIONPROC		pglGetAttribLocation;
	extern PFNGLPATCHPARAMETERIPROC			pglPatchParameteri;
	extern PFNGLDRAWBUFFERSPROC				pglDrawBuffers;
	extern PFNGLBINDIMAGETEXTUREPROC		pglBindImageTexture;
	extern PFNGLDISPATCHCOMPUTEPROC			pglDispatchCompute;
	extern PFNGLDISPATCHCOMPUTEINDIRECTPROC pglDispatchComputeIndirect;
	#define glCreateShader					pglCreateShader
	#define glShaderSource					pglShaderSource
	#define glCompileShader					pglCompileShader
	#define glCreateProgram					pglCreateProgram
	#define glAttachShader					pglAttachShader
	#define glLinkProgram					pglLinkProgram
	#define glUseProgram					pglUseProgram
	#define glDeleteShader					pglDeleteShader
	#define glDeleteProgram					pglDeleteProgram
	#define glGetShaderiv					pglGetShaderiv
	#define glGetShaderInfoLog				pglGetShaderInfoLog
	#define glGetProgramiv					pglGetProgramiv
	#define glGetProgramInfoLog				pglGetProgramInfoLog
	#define glGetUniformLocation			pglGetUniformLocation
	#define glUniform1f						pglUniform1f
	#define glUniform2f						pglUniform2f 
	#define glUniform3f						pglUniform3f 
	#define glUniform4f						pglUniform4f
	#define glUniform1fv					pglUniform1fv
	#define glUniform2fv					pglUniform2fv
	#define glUniform3fv					pglUniform3fv
	#define glUniform4fv					pglUniform4fv
	#define glUniform1i						pglUniform1i
	#define glUniformMatrix2fv				pglUniformMatrix2fv
	#define glUniformMatrix3fv				pglUniformMatrix3fv
	#define glUniformMatrix4fv				pglUniformMatrix4fv
	#define glVertexAttribPointer			pglVertexAttribPointer
	#define glEnableVertexAttribArray		pglEnableVertexAttribArray
	#define glDisableVertexAttribArray		pglDisableVertexAttribArray
	#define glBindAttribLocation			pglBindAttribLocation
	#define glGetAttribLocation				pglGetAttribLocation
	#define glPatchParameteri				pglPatchParameteri
	#define glDrawBuffers					pglDrawBuffers
	#define glBindImageTexture				pglBindImageTexture
	#define glDispatchCompute				pglDispatchCompute
	#define glDispatchComputeIndirect		pglDispatchComputeIndirect
	PFNGLCREATESHADERPROC					glCreateShader = 0;
	PFNGLSHADERSOURCEPROC					glShaderSource = 0;
	PFNGLCOMPILESHADERPROC					glCompileShader = 0;
	PFNGLCREATEPROGRAMPROC					glCreateProgram = 0;
	PFNGLATTACHSHADERPROC					glAttachShader = 0;
	PFNGLLINKPROGRAMPROC					glLinkProgram = 0;
	PFNGLUSEPROGRAMPROC						glUseProgram = 0;
	PFNGLDELETESHADERPROC					glDeleteShader = 0;
	PFNGLDELETEPROGRAMPROC					glDeleteProgram = 0;
	PFNGLGETSHADERIVPROC					glGetShaderiv = 0;
	PFNGLGETSHADERINFOLOGPROC				glGetShaderInfoLog = 0;
	PFNGLGETPROGRAMIVPROC					glGetProgramiv = 0;
	PFNGLGETPROGRAMINFOLOGPROC				glGetProgramInfoLog = 0;
	PFNGLGETUNIFORMLOCATIONPROC				glGetUniformLocation = 0;
	PFNGLUNIFORM1FPROC						glUniform1f = 0;
	PFNGLUNIFORM2FPROC						glUniform2f = 0;
	PFNGLUNIFORM3FPROC						glUniform3f = 0;
	PFNGLUNIFORM4FPROC						glUniform4f = 0;
	PFNGLUNIFORM1FVPROC						glUniform1fv = 0;
	PFNGLUNIFORM2FVPROC						glUniform2fv = 0;
	PFNGLUNIFORM3FVPROC						glUniform3fv = 0;
	PFNGLUNIFORM4FVPROC						glUniform4fv = 0;
	PFNGLUNIFORM1IPROC						glUniform1i = 0;
	PFNGLUNIFORMMATRIX2FVPROC				glUniformMatrix2fv = 0;
	PFNGLUNIFORMMATRIX3FVPROC				glUniformMatrix3fv = 0;
	PFNGLUNIFORMMATRIX4FVPROC				glUniformMatrix4fv = 0;
	PFNGLVERTEXATTRIBPOINTERPROC			glVertexAttribPointer = 0;
	PFNGLENABLEVERTEXATTRIBARRAYPROC		glEnableVertexAttribArray = 0;
	PFNGLDISABLEVERTEXATTRIBARRAYPROC		glDisableVertexAttribArray = 0;
	PFNGLBINDATTRIBLOCATIONPROC				glBindAttribLocation = 0;
	PFNGLGETATTRIBLOCATIONPROC				pglGetAttribLocation = 0;
	PFNGLPATCHPARAMETERIPROC				pglPatchParameteri = 0;
	PFNGLDRAWBUFFERSPROC					pglDrawBuffers = 0;
	PFNGLBINDIMAGETEXTUREPROC				pglBindImageTexture = 0;
	PFNGLDISPATCHCOMPUTEPROC				pglDispatchCompute = 0;
	PFNGLDISPATCHCOMPUTEINDIRECTPROC		pglDispatchComputeIndirect = 0;
#endif

namespace GL
{
    namespace GLSL
    {
        double version()
        {
	        return atof((const char *)glGetString(GL_SHADING_LANGUAGE_VERSION));
        }

        bool available()
        {
	        static int supported = -1;

	        if( supported==-1 )
	        {
		        supported = GL::isExtensionSupported("GL_ARB_shader_objects");

		        supported = supported						 
			        && GL::isExtensionSupported("GL_ARB_vertex_shader")
			        && GL::isExtensionSupported("GL_ARB_fragment_shader")
			        && GL::isExtensionSupported("GL_ARB_shading_language_100");

		        #if defined(_WIN32)
		        supported = supported
			        && (glCreateShader =			(PFNGLCREATESHADERPROC)				getProc("glCreateShader",			"glCreateShaderObjectARB"))
			        && (glShaderSource =			(PFNGLSHADERSOURCEPROC)				getProc("glShaderSource",			"glShaderSourceARB"))
			        && (glCompileShader =			(PFNGLCOMPILESHADERPROC)			getProc("glCompileShader",			"glCompileShaderARB"))
			        && (glCreateProgram =			(PFNGLCREATEPROGRAMPROC)			getProc("glCreateProgram",			"glCreateProgramObjectARB"))
			        && (glAttachShader =			(PFNGLATTACHSHADERPROC)				getProc("glAttachShader",			"glAttachObjectARB"))
			        && (glLinkProgram =				(PFNGLLINKPROGRAMPROC)				getProc("glLinkProgram",			"glLinkProgramARB"))
			        && (glUseProgram =				(PFNGLUSEPROGRAMPROC)				getProc("glUseProgram",				"glUseProgramObjectARB"))
			        && (glDeleteShader =			(PFNGLDELETESHADERPROC)				getProc("glDeleteShader",			"glDeleteObjectARB"))
			        && (glDeleteProgram =			(PFNGLDELETEPROGRAMPROC)			getProc("glDeleteProgram",			"glDeleteObjectARB"))
			        && (glGetShaderiv =				(PFNGLGETSHADERIVPROC)				getProc("glGetShaderiv",			"glGetObjectParameterivARB"))
			        && (glGetShaderInfoLog =		(PFNGLGETSHADERINFOLOGPROC)			getProc("glGetShaderInfoLog",		"glGetInfoLogARB"))
			        && (glGetProgramiv =			(PFNGLGETPROGRAMIVPROC)				getProc("glGetProgramiv",			"glGetObjectParameterivARB"))
			        && (glGetProgramInfoLog =		(PFNGLGETPROGRAMINFOLOGPROC)		getProc("glGetProgramInfoLog",		"glGetInfoLogARB"))
			        && (glGetUniformLocation =		(PFNGLGETUNIFORMLOCATIONPROC)		getProc("glGetUniformLocation",		"glGetUniformLocationARB"))
			        && (glUniform1f =				(PFNGLUNIFORM1FPROC)				getProc("glUniform1f",				"glUniform1fARB"))
			        && (glUniform2f =				(PFNGLUNIFORM2FPROC)				getProc("glUniform2f",				"glUniform2fARB"))
			        && (glUniform3f =				(PFNGLUNIFORM3FPROC)				getProc("glUniform3f",				"glUniform3fARB"))
			        && (glUniform4f =				(PFNGLUNIFORM4FPROC)				getProc("glUniform4f",				"glUniform4fARB"))
			        && (glUniform1fv =				(PFNGLUNIFORM1FVPROC)				getProc("glUniform1fv",				"glUniform1fvARB"))
			        && (glUniform2fv =				(PFNGLUNIFORM2FVPROC)				getProc("glUniform2fv",				"glUniform2fvARB"))
			        && (glUniform3fv =				(PFNGLUNIFORM3FVPROC)				getProc("glUniform3fv",				"glUniform3fvARB"))
			        && (glUniform4fv =				(PFNGLUNIFORM4FVPROC)				getProc("glUniform4fv",				"glUniform4fvARB"))		
			        && (glUniform1i =				(PFNGLUNIFORM1IPROC)				getProc("glUniform1i",				"glUniform1iARB"))
			        && (glUniformMatrix2fv =		(PFNGLUNIFORMMATRIX2FVPROC)			getProc("glUniformMatrix2fv",		"glUniformMatrix2fvARB"))
			        && (glUniformMatrix3fv =		(PFNGLUNIFORMMATRIX3FVPROC)			getProc("glUniformMatrix3fv",		"glUniformMatrix3fvARB"))
			        && (glUniformMatrix4fv =		(PFNGLUNIFORMMATRIX4FVPROC)			getProc("glUniformMatrix4fv",		"glUniformMatrix4fvARB"))
			        && (glVertexAttribPointer =		(PFNGLVERTEXATTRIBPOINTERPROC)		getProc("glVertexAttribPointer",	"glVertexAttribPointerARB"))
			        && (glEnableVertexAttribArray =	(PFNGLENABLEVERTEXATTRIBARRAYPROC)	getProc("glEnableVertexAttribArray","glEnableVertexAttribArrayARB"))
			        && (glDisableVertexAttribArray =(PFNGLDISABLEVERTEXATTRIBARRAYPROC)	getProc("glDisableVertexAttribArray","glDisableVertexAttribArrayARB"))
			        && (glBindAttribLocation =		(PFNGLBINDATTRIBLOCATIONPROC)		getProc("glBindAttribLocation",		"glBindAttribLocationARB"))
			        && (glGetAttribLocation =		(PFNGLGETATTRIBLOCATIONPROC)		getProc("glGetAttribLocation",		"glGetAttribLocationARB"))
			        && (glDrawBuffers =				(PFNGLDRAWBUFFERSPROC)				getProc("glDrawBuffers",			"glDrawBuffersARB",					"glDrawBuffersATI"))
		        ;
		        #endif

		        //printf( "%-8s: %s\n", "GLSL", supported ? "OK" : "NOT AVAILABLE" );
	        }

	        return supported;
        }

        bool tessellatorAvailable()
        {
	        static int supported = -1;
	        if(supported==-1) {
		        supported = 
			            GL::isExtensionSupported("GL_ARB_texture_float") 
			        &&  GL::isExtensionSupported("GL_EXT_gpu_shader4")  // this should enable gl_vertexID
			        &&  GL::version() >= 4.1
			        && (glPatchParameteri =	(PFNGLPATCHPARAMETERIPROC) getProc("glPatchParameteri",	"glPatchParameteriARB"))
		        ;
	        }
	        return supported;
        }

        bool computeAvailable()
        {
	        static int supported = -1;

	        if( supported==-1 )
	        {
		        supported = 
			            GL::isExtensionSupported("GL_ARB_compute_shader") 
			        && (GL::isExtensionSupported("GL_ARB_shader_image_load_store") || GL::isExtensionSupported("GL_EXT_shader_image_load_store"))
			        &&  GL::version() >= 4.3
			        && (glBindImageTexture			= (PFNGLBINDIMAGETEXTUREPROC)			getProc("glBindImageTexture",	"glBindImageTextureEXT"))
			        && (glDispatchCompute			= (PFNGLDISPATCHCOMPUTEPROC)			getProc("glDispatchCompute"))
			        && (glDispatchComputeIndirect	= (PFNGLDISPATCHCOMPUTEINDIRECTPROC)	getProc("glDispatchComputeIndirect"))
		        ;
	        }

	        return supported;
        }

        void deallocateShader( uint& id )
        {
	        if( id ) glDeleteShader(id);
	        id = 0;
        }

        void deallocateProgram( uint& id )
        {
	        if( id ) glDeleteProgram(id);
	        id = 0;
        }

        bool attach( Program& program, char* source, int len, int type )
        {
	        for(; *source<=' '; source++) len--;
	        if(!source || len<=0) return false;

	        if(len>MAX_SHADER_SOURCE_LENGTH) {
		        DEBUG_CRITICAL("Length of shader source exceeds buffer");
	        }

	        char buffer[MAX_SHADER_SOURCE_LENGTH];
	        strncpy( buffer, source, len );
	        buffer[len] = 0;

	        Cacheable& shader = program.push().tail();
	        shader.id = glCreateShader( type );
	        shader.deallocator = deallocateShader;

	        glShaderSource( shader.id, 1, (const char**)&source, &len );
	        glCompileShader( shader.id );

	        int successful = 0;
	        glGetShaderiv( program.tail().id, GL_COMPILE_STATUS, &successful );

	        if( successful )
	        {
		        glAttachShader( program.id, program.tail().id );
	        }
        	
	        return successful;
        }

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
		        if(nIndex==0) break;

                if((focusLine < 0) || (focusLine >= 0) && (focusLine-6 <= line) && (line <= focusLine+3)) { 
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

        void checkProgram( Program& program )
        {
	        int successful;
	        glGetProgramiv( program.id, GL_LINK_STATUS, &successful);
	        if( successful ) return;

	        char message[512];
	        glGetProgramInfoLog( program.id, sizeof(message), 0, message );
	        DEBUG_CRITICAL( message );
        }

        void build( Program& program, char* source )
        {
            //PROFILE;
            DEBUG_PATH;

	        if(!source) {
		        DEBUG_CRITICAL("Source null");
	        }

	        if(strlen(source) > MAX_PROGRAM_SOURCE_LENGTH) {
		        DEBUG_CRITICAL("Length of program source exceeds buffer");
	        }

	        if( program.id==0 && source )
	        {
		        program.id = glCreateProgram();
		        program.deallocator = deallocateProgram;

		        char* shader_names[] = { "VERTEX", "CONTROL", "EVALUATION", "GEOMETRY", "FRAGMENT", "COMPUTE" };
		        int   shader_types[] = { GL_VERTEX_SHADER, GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER, GL_COMPUTE_SHADER };
		        int num_of_shaders = sizeof(shader_types)/sizeof(int);

		        char* shader_sources[20];
		        for(int i=0; i<num_of_shaders; i++) {
			        char temp[20];
			        strcpy(temp, shader_names[i]);
			        strcat(temp,":\r");
			        shader_sources[i] = strstr2( source, temp );
		        }

		        for( int i=0; i<num_of_shaders; i++ )
		        {
			        if( shader_sources[i] != NULL )
			        {
				        shader_sources[i] = strstr2( shader_sources[i], "\n" ) + 1;

				        int shader_lenght = 0;
				        for( int j=i+1; j<num_of_shaders; j++ )
				        {
					        if( shader_sources[j] != NULL ) 
					        {
						        shader_lenght = (int)(shader_sources[j]-shader_sources[i]);
						        break;
					        }
				        }

				        if( shader_lenght==0 )
				        {
					        shader_lenght = strlen( shader_sources[i] );
				        }
        			
				        bool hasShaderCompiled = attach(program, shader_sources[i], shader_lenght, shader_types[i]);
				        color(CMD_WHITE,CMD_CYAN); DEBUG_PRINT("%s", shader_names[i]);
						color(CMD_WHITE, 0); DEBUG_PRINT(" ");
                        
				        if( !hasShaderCompiled )
				        {
							color(CMD_LIGHTRED,0); DEBUG_PRINT("%s", "==> ERROR\n\n");

                            char message[512];
	                        glGetShaderInfoLog( program.tail().id, sizeof(message), NULL, message );
							printf("%s\n", message);

							// Parse to get the error line.
							// NOTE: Different drivers give different error messages
							// https://gamedev.stackexchange.com/questions/38685/is-this-a-reliable-method-of-parsing-glgetshaderinfolog
							int line = -1;
							if(line <= 0) {
								// NVIDIA error message
								char *p1 = 0, *p2 = 0;
								p1 = strchr(message, '(');
								if( p1>0 ) p2 = strchr(p1+1, ')');
								if( p1>0 && p2>0 ) line = strtol(p1+1, &p2, 10);
							}
							if(line <= 0) {
								// ATI/Intel error message
								char *p0 = 0, *p1 = 0, *p2 = 0;
								p0 = strchr(message, ':');
								if( p0>0 ) p1 = strchr(p0+1, ':');
								if( p1>0 ) p2 = strchr(p1+1, ':');
								if(p0>0 && p1>0 && p2>0) line = strtol(p1+1, &p2, 10);
							}
							ASSERT(line > 0);

                            #if defined(DEBUG_SHOW_GLSL_SOURCE)
					            printSource( shader_sources[i], shader_lenght, line ); 
                                printf("\n");
                            #endif

                            //CRITICAL("Shader compile error");
				        }
			        }
		        }

		        glLinkProgram( program.id );
		        checkProgram( program );
	        }
            DEBUG_PRINT("\n");
        }

        void bind( Program& program )
        {
            DEBUG_ASSERT( program.id != 0 );
            if( program.id != 0 ) {
		        glUseProgram( program.id );
            }
        }

        void unbind()
        {
            if( available() ) {
		        glUseProgram( 0 );
            }
        }

        int getLocation( Program& program, char* name )
        {
	        int location = glGetUniformLocation( program.id, name );

	        if( location<0 )
	        {	
		        char buffer[200];
		        sprintf( buffer, "Uniform '%s' not found.\n", name );
		        DEBUG_CRITICAL( buffer );
	        }
	        return location;
        }

        void set( Program& program, char* name, int value )
        {
	        glUniform1i( getLocation(program, name), value );
        }

        void set(  Program& program, char* name, float value )
        {
	        glUniform1f( getLocation(program, name), value );
        }

        void set(  Program& program, char* name, vec2 vector )
        {
	        glUniform2fv( getLocation(program, name), 1, vector.array );
        }

        void set(  Program& program, char* name, vec3& vector )
        {
	        glUniform3fv( getLocation(program, name), 1, vector.array );
        }

        void set(  Program& program, char* name, vec4& vector )
        {
	        glUniform4fv( getLocation(program, name), 1, vector.array );
        }

        void set(  Program& program, char* name, mat4& matrix )
        {
	        glUniformMatrix4fv( getLocation(program, name), 1, 0, matrix.array );
        }

        void set(  Program& program, char* name, mat3& matrix )
        {
	        glUniformMatrix3fv( getLocation(program, name), 1, 0, matrix.array );
        }
    }
}

