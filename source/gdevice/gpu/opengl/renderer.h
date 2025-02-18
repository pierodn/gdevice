#pragma once

#include "type/glsl.h"
#include "type/array.h"
#include "type/scene/node.h"

#include "os/profiler.h"
#include "os/timer.h"

#include "__temp/light.h"
#include "gpu/program.h"
#include "gpu/opengl/gl.h"
#include "gpu/controls.h"

#include "application/assets/worlds/planet1/parameters.h" // CLIPMAP_WINDOW, TEXTURE_RANGE
#include "application/assets/worlds/planet1/textures.h"
#include "application/assets/worlds/planet1/generate_terrain.glsl.h"
#include "application/assets/worlds/planet1/render_terrain.glsl.h"
#include "application/assets/worlds/planet1/render_sky.glsl.h"

class Renderer
{
public:
	Controls controls;

	Renderer()
	{
	}

	void initialize()
	{
		//
		// Hardware detection  
        // 
        DEBUG_TRACE(GL::renderer());
        DEBUG_TRACE(GL::version());
        DEBUG_PRINT("\n");

        DEBUG_ASSERT( GL::version() >= 4.3 );

	    GL::Texturing::available();
	    GL::Texturing::textureNonPowerOfTwoAvailable();
	    GL::secondaryColorAvailable();
	    GL::swapControlAvailable();
	    GL::VBO::available();
	    GL::GLSL::available();
	    GL::GLSL::tessellatorAvailable();
	    GL::GLSL::computeAvailable();
	    GL::MRT::available(); 
/*
		if(true) 
        {
            int maxTextureUnits;
			glGetIntegerv(GL_MAX_TEXTURE_UNITS, &maxTextureUnits);
			//DEBUG_PRINT( TAB32 ": %i\n", "Max texture units", maxTextureUnits );
		}

		if( GL::GLSL::available() ) 
        {
            int maxTextureImageUnits;
			glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &maxTextureImageUnits);
			//DEBUG_PRINT( TAB32 ": %i\n", "Max texture image units", maxTextureImageUnits );
		}

		if(glEnableVertexAttribArray) 
        {
			glGetIntegerv( GL_MAX_VERTEX_ATTRIBS, &maxVertexAttributes );
			DEBUG_PRINT( TAB32 ": %i\n", "Max vertex attributes", maxVertexAttributes );
		}

	    if( GL::GLSL::tessellatorAvailable() ) 
        {
            int maxPatchVertices;
		    glGetIntegerv(GL_MAX_PATCH_VERTICES, &maxPatchVertices);
		    //DEBUG_PRINT( TAB32 ": %i\n", "Max patch vertices", maxPatchVertices );
	    }

		if(true)
        {
			int maxDrawBuffers;
			glGetIntegerv(GL_MAX_DRAW_BUFFERS, &maxDrawBuffers);
			//DEBUG_PRINT( TAB32 ": %i\n", "Max draw buffers", maxDrawBuffers );
		}
*/
		//
		// Initialize GL state 
		//
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

		GL::Texturing::unbind();

		glMatrixMode(GL_TEXTURE);
		glLoadIdentity();

        if( GL::swapControlAvailable() ) {
            GL::swapControl(0); // TODO check
        }
		
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

        //
		// Initialize assets
		//
		initializeTerrainDetails();

		//
		// Initialize function
		//
        GL::GLSL::unbind();

        DEBUG_PRINT("\n");

		//
		// Flags
		// 
		controls.initialize();	// TEMP	
		controls.showLegenda();
	}

    vec2 GetViewport()
    {
        GLint viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);

        return vec2(viewport[2], viewport[3]);
    }
    
	void SetTarget( vec2 viewport )
	{
        glViewport( 0, 0, viewport.x, viewport.y );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	}

};
