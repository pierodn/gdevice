#pragma once

//
// Glink - OpenGL Dynamic linker
// 

#define GL_GLEXT_PROTOTYPES
#ifdef _WIN32
	#include <gl/gl.h>
	#include <gl/glext.h>
#else
	#include <OpenGL/glext.h>
#endif



const bool Verbose = false;

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
	inline char* GetVendor()
	{
		return (char*)glGetString(GL_VENDOR);
	}

    inline char* GetRenderer()
	{
        return (char*)glGetString(GL_RENDERER);
    }

	float GetVersion()
	{
		return atof((char*)glGetString(GL_VERSION));
	}

	bool IsExtensionSupported( const char *extension )
	{
        std::string extensions = (char*) glGetString(GL_EXTENSIONS);

        bool found = (extensions.find(extension) != std::string::npos);

        if(Verbose)
        {
            color(CMD_DARKGRAY, 0);
		    DEBUG_PRINT(TAB32 ": %s\n", extension, found ? "OK" : "MISSING" );
        }
		return found;
	}
/*
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
    */

    const int FirstTime = -1;

    /// ============================
    //  Texturing
    /// ============================

	namespace Texturing
	{
		bool available()
		{
			static int supported = FirstTime;

			if( supported == FirstTime )
			{
				supported = IsExtensionSupported("GL_ARB_multitexture");
						 
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
					&& IsExtensionSupported("GL_ARB_texture_env_combine")
					&& IsExtensionSupported("GL_ARB_texture_env_dot3");

                if(Verbose)
                {
                    color(CMD_DARKGRAY, 0);
				    DEBUG_PRINT(TAB32 ": %s\n", "FIXED", supported ? "OK" : "NOT AVAILABLE" );
                }
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


		void deallocate( uint& id )
		{
			glDeleteTextures(1, &id );
		}

		void unbind( int slot )
		{
			glActiveTexture( GL_TEXTURE0 + slot );
			glDisable(GL_TEXTURE_2D);
		}

	};
};

/*
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
		static int supported = FirstTime;

		if( supported == FirstTime )
		{
			supported = IsExtensionSupported("GL_EXT_secondary_color" );
		
			#if defined(_WIN32)
			supported = supported
				&& (glSecondaryColorPointer	= (PFNGLSECONDARYCOLORPOINTERPROC) getProc("glSecondaryColorPointer", "glSecondaryColorPointerEXT") )
			;
			#endif
            
            if(Verbose)
            {
                color(CMD_DARKGRAY, 0);
			    DEBUG_PRINT(TAB32 ": %s\n", "SECONDARY COLOR", supported ? "OK" : "NOT AVAILABLE" );
		    }
        }

		return supported;
	}
};
*/


/// ============================
//  WGL_EXT_swap_control
/// ============================

#if defined(_WIN32)
	typedef void (APIENTRY *PFNWGLEXTSWAPCONTROLPROC)(int);
	typedef int (*PFNWGLEXTGETSWAPINTERVALPROC)(void);
	PFNWGLEXTSWAPCONTROLPROC		wglSwapIntervalEXT = 0;
	PFNWGLEXTGETSWAPINTERVALPROC	wglGetSwapIntervalEXT = 0;
#endif

// TODO rename as SetSwapBuffersCount, GetSwapBuffersCount
namespace GL
{ 
	bool swapControlAvailable()
	{
		static int supported = FirstTime;
		if( supported == FirstTime )
        {
			supported = IsExtensionSupported( "WGL_EXT_swap_control" );
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
	        static int supported = FirstTime;

	        if( supported == FirstTime )
	        {
		        supported = GL::IsExtensionSupported("GL_EXT_framebuffer_object") 
			             || GL::IsExtensionSupported("GL_ARB_framebuffer_object");

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

                if(Verbose)
                {
                    color(CMD_DARKGRAY, 0);
		            DEBUG_PRINT(TAB32 ": %s\n", "MRT", supported ? "OK" : "NOT AVAILABLE" );
                }
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

        void unbindTarget()
        {
	        glPopAttrib();
	        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }
        */
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
	        static int supported = FirstTime;

	        if( supported == FirstTime )
	        {
		        supported = GL::IsExtensionSupported("GL_EXT_transform_feedback") 
			             || GL::IsExtensionSupported("GL_NV_transform_feedback");	// uguale a EXT?	
				         //|| IsExtensionSupported("GL_ARB_transform_feedback2"); // addictional functionalities

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

                if(Verbose)
                {
                    color(CMD_DARKGRAY, 0);
		            DEBUG_PRINT(TAB32 ": %s\n", "RTT", supported ? "OK" : "NOT AVAILABLE" );
                }
	        }

	        return supported;
        }
    }
}

// TEMP
namespace GL 
{
    vec2 GetViewport()
    {
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);

        return vec2(viewport[2], viewport[3]);
    }

    void SetViewport(vec2& viewport)
    {
        glViewport( 0, 0, viewport.x, viewport.y );
    }

    void ClearBuffer() // TODO buffer param
	{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
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
	        static int supported = FirstTime;

	        if( supported == FirstTime )
	        {
		        // ?
		        supported = GL::IsExtensionSupported("GL_ARB_vertex_buffer_object");

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

                if(Verbose)
                {
                    color(CMD_DARKGRAY, 0);
		            DEBUG_PRINT(TAB32 ": %s\n", "VBO", supported ? "OK" : "NOT AVAILABLE" );
                }
	        }

	        return supported;
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
    extern PFNGLMEMORYBARRIERPROC           pglMemoryBarrier;
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
    #define glMemoryBarrier         		pglMemoryBarrier
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
    PFNGLMEMORYBARRIERPROC             		pglMemoryBarrier = 0;
#endif



namespace GL
{
    namespace GLSL
    {
        bool available()
        {
            static int supported = FirstTime;

            if( supported == FirstTime )
            {
	            supported = GL::IsExtensionSupported("GL_ARB_shader_objects");

	            supported = supported						 
		            && GL::IsExtensionSupported("GL_ARB_vertex_shader")
		            && GL::IsExtensionSupported("GL_ARB_fragment_shader")
		            && GL::IsExtensionSupported("GL_ARB_shading_language_100");

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

                if(Verbose)
                {
                    color(CMD_DARKGRAY, 0);
	                DEBUG_PRINT(TAB32 ": %s\n", "GLSL", supported ? "OK" : "NOT AVAILABLE" );
                }
            }

            return supported;
        }

        bool tessellatorAvailable()
        {
            static int supported = FirstTime;
            if( supported == FirstTime )
            {
	            supported = 
		                GL::IsExtensionSupported("GL_ARB_texture_float") 
		            &&  GL::IsExtensionSupported("GL_EXT_gpu_shader4")  // this should enable gl_vertexID
		            &&  GL::GetVersion() >= 4.1
		            && (glPatchParameteri =	(PFNGLPATCHPARAMETERIPROC) getProc("glPatchParameteri",	"glPatchParameteriARB"))
	            ;
            }
            return supported;
        }

        bool computeAvailable()
        {
            static int supported = FirstTime;

            if( supported == FirstTime )
            {
	            supported = 
		                GL::IsExtensionSupported("GL_ARB_compute_shader") 
		            && (GL::IsExtensionSupported("GL_ARB_shader_image_load_store") || GL::IsExtensionSupported("GL_EXT_shader_image_load_store"))
		            &&  GL::GetVersion() >= 4.3
		            && (glBindImageTexture			= (PFNGLBINDIMAGETEXTUREPROC)			getProc("glBindImageTexture",	"glBindImageTextureEXT"))
		            && (glDispatchCompute			= (PFNGLDISPATCHCOMPUTEPROC)			getProc("glDispatchCompute"))
		            && (glDispatchComputeIndirect	= (PFNGLDISPATCHCOMPUTEINDIRECTPROC)	getProc("glDispatchComputeIndirect"))
                    && (glMemoryBarrier			    = (PFNGLMEMORYBARRIERPROC)			    getProc("glMemoryBarrier"))
	            ;
            }

            return supported;
        }
    } // GLSL

}





namespace GL
{
    namespace GLSL
    {
        void Check(int requiredVersion = 4.3)
	    {
            DEBUG_TRACE(GL::GetRenderer());
            DEBUG_TRACE(GL::GetVersion());
            DEBUG_ASSERT(GL::GetVersion() >= requiredVersion);

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

    // TODO void Link(int minimumVersion = 4.3)
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
			//GL::VBO::unbind(GL_ARRAY_BUFFER);
			//GL::VBO::unbind(GL_ELEMENT_ARRAY_BUFFER);
            glBindBuffer( GL_ARRAY_BUFFER, 0 );
            glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );
		}

        //GL::UseProgram(0);
        glUseProgram(0);
	}

}