#pragma once

// TODO remove dependency on these:
//#include "platforms/platform.h"
#include "parameters.h" // CLIPMAP_WINDOW, TEXTURE_RANGE
///////////////////////////////////////


#include "platforms/profiler.h"
#include "libraries/math/algebra.h"

#include "platforms/timer.h"

#include "engines/graphic/array.h"
#include "engines/graphic/material.h"
#include "engines/graphic/light.h"
#include "engines/graphic/program.h"
#include "engines/graphic/node.h"

//#include "engines/graphic/opengl/shaders/generate_gradientmap.glsl.h"
#include "engines/graphic/opengl/shaders/generate_terrain.glsl.h"
#include "engines/graphic/opengl/shaders/render_terrain.glsl.h"
#include "engines/graphic/opengl/shaders/render_sky.glsl.h"



#include "engines/graphic/opengl/gl.h"
//#include "engines/graphic/opengl/glsl.h"
//#include "engines/graphic/opengl/vbo.h"
//#include "engines/graphic/opengl/mrt.h"

#include "engines/graphic/asset.h"
#include "engines/graphic/controls.h"

#define MAX_LIGHTS 4 // 32


//  _____           _                 
// | __  |___ ___ _| |___ ___ ___ ___ 
// |    -| -_|   | . | -_|  _| -_|  _|
// |__|__|___|_|_|___|___|_| |___|_|  
//
class Renderer
{
public:
	Controls controls;

	mat4 ProjectionMatrix;
	mat4 ModelViewMatrix;
	mat3 NormalMatrix;
	mat4 ModelViewProjectionMatrix;
	mat4 inverseRotationMatrix;
	Array<mat4> ModelViewMatrixStack;

	Material* material;
	vec4 ambientLight;
	vec4 indirectLight;
	Light* lights[MAX_LIGHTS];

	float fog_density;

    Program programGenerateTerrain;
    Program programGenerateGradientMap;
	Program programRenderTerrain;
	Program programRenderSky;


	int maxVertexAttributes;
	int maxPatchVertices;
	int maxTextureUnits;
	int maxTextureImageUnits;


	Renderer()
	{
		ModelViewMatrixStack.allocate(2000,100);
	}

	void initialize()
	{
		//
		// Hardware detection  
        // 
#if defined(_DEBUG) // TEMP
		printh( "GRAPHICS ENGINE" );
		DEBUG_PRINT( TAB32 ": %s\n",   "OpenGL vendor",  GL::vendor() );
		DEBUG_PRINT( TAB32 ": %.1f\n", "OpenGL version", GL::version() );
		DEBUG_PRINT( TAB32 ": %.2f\n", "GLSL version",   GL::GLSL::version() );
#else 
        color(15,0); 
        DEBUG_PRINT( "OpenGL %.1f\n\n", GL::version() );
#endif

	    GL::Texturing::available();
	    GL::Texturing::textureNonPowerOfTwoAvailable();
	    GL::secondaryColorAvailable();
	    GL::swapControlAvailable();
	    GL::VBO::available();
	    GL::GLSL::available();
	    GL::GLSL::tessellatorAvailable();
	    GL::GLSL::computeAvailable();
	    GL::MRT::available(); 

		if(true) {
			glGetIntegerv(GL_MAX_TEXTURE_UNITS, &maxTextureUnits);
			//DEBUG_PRINT( TAB32 ": %i\n", "Max texture units", maxTextureUnits );
		}

		if( GL::GLSL::available() ) {
			glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, &maxTextureImageUnits);
			//DEBUG_PRINT( TAB32 ": %i\n", "Max texture image units", maxTextureImageUnits );
		}
/*
		if(glEnableVertexAttribArray) {
			glGetIntegerv( GL_MAX_VERTEX_ATTRIBS, &maxVertexAttributes );
			DEBUG_PRINT( TAB32 ": %i\n", "Max vertex attributes", maxVertexAttributes );
		}
*/
	    if( GL::GLSL::tessellatorAvailable() ) {
		    glGetIntegerv(GL_MAX_PATCH_VERTICES, &maxPatchVertices);
		    //DEBUG_PRINT( TAB32 ": %i\n", "Max patch vertices", maxPatchVertices );
	    }

		if(true) {
			int maxDrawBuffers;
			glGetIntegerv(GL_MAX_DRAW_BUFFERS, &maxDrawBuffers);
			//DEBUG_PRINT( TAB32 ": %i\n", "Max draw buffers", maxDrawBuffers );
		}

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
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		glDrawBuffer( GL_BACK );

		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );

		// BACKFACES ?
		(false ? glDisable : glEnable)( GL_CULL_FACE );

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
        if( GL::VBO::available() ) {
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
		DEBUG_ASSERT( GL::version() >= 4.3 );

        GL::GLSL::unbind();

        DEBUG_PRINT("\n");
		GL::GLSL::build( programRenderTerrain, render_terrain_glsl );
		GL::GLSL::build( programRenderSky, render_sky_glsl );
		GL::GLSL::build(programGenerateTerrain, generate_terrain_glsl);

		//
		// Flags
		// 
		controls.initialize();	// TEMP	
		controls.showLegenda();
	}
	

	// 
	// Hardware capabilities
	//

	bool isNonPowerOfTwoTexturesAvailable()
	{
		return GL::Texturing::textureNonPowerOfTwoAvailable();
	}


	//
	// Internal state
	//
	void setViewport( vec2& viewport )
	{
		glViewport( 0, 0, viewport.width, viewport.height );
	}

	vec2 getViewport()
	{
		GLint temp[4];
		glGetIntegerv( GL_VIEWPORT, temp );
		return vec2(temp[2]-temp[0], temp[3]-temp[1]);
	}

	void setProjection( double fov, double near_plane, double far_plane )
	{
		ProjectionMatrix = projection(getViewport(), fov, near_plane, far_plane);
		glMatrixMode( GL_PROJECTION );
		glLoadMatrixf( ProjectionMatrix.array );
		// TODO frustum culling
	}

	void setTransform( mat4& transform )
	{
		glMatrixMode( GL_MODELVIEW );
		glLoadMatrixf( transform.array );
	}

	void setFogDensity( float density )
	{
		fog_density = density;
	}

	void setAmbientLight( vec4& ambientLight )
	{
		this->ambientLight = ambientLight;
		glLightModelfv( GL_LIGHT_MODEL_AMBIENT, ambientLight.array );
	}

	int vertices;


	void clear( vec4& color = vec4(0) )
	{
		vertices = 0;

		// TEMP 
		// glPolygonMode( GL_FRONT_AND_BACK, _getFlag(WIREFRAME) ? GL_LINE : GL_FILL );
		glPolygonMode( GL_FRONT_AND_BACK, false ? GL_LINE : GL_FILL );

		glClearColor( color.r, color.g, color.b , color.a );	//
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		for( int i=0; i<MAX_LIGHTS; i++ ) {
			lights[i] = NULL;
		}

		ModelViewMatrix = mat4(1);
	}









public:

	void draw( Node& node )
	{
        // ****  NOTE  **** ////////////
        // Do not send in pipeline a node the very first time is generated as it might not be updated (and valid) yet.
        // From 2nd generation on, the node might be updated or not, yet valid. In case it is not updated, it will be the previous version.
        // This is workaround to avoid waiting for glDispatchCompute() to finish (syncronization GPU/CPU).
        bool skipRenderingThisNodeThisTime = false;
        ////////////////////////////////

		ModelViewMatrixStack.push( ModelViewMatrix );
		ModelViewMatrix *= TransformationMatrix( node.transform );
		ModelViewProjectionMatrix = ProjectionMatrix * ModelViewMatrix;
		NormalMatrix = inverseTranspose( mat3(ModelViewMatrix) );
        /*
            // TEST Removing scale component from NormalMatrix to fix the gradient discountinuity problem.
            mat4 ModelViewMatrixWithoutScale = ModelViewMatrixStack.tail();
                // TransformationMatrix(node.transform);
            	mat4 t = mat4(1.0); //mat4(transform.scale);
                t.translation.xyz = //transform.scale * 
                                    node.transform.position;
                t = RotationMatrix(node.transform.rotation) * t;
                ModelViewMatrixWithoutScale *= t;
            NormalMatrix = inverseTranspose( mat3(ModelViewMatrixWithoutScale) );
        */


		setTransform( ModelViewMatrix );

		if( node.light ) {
			lights[ node.light->id ] = node.light;
		}
/*
		if( node.material ) {
			material = node.material;

			indirectLight	= material->emission + material->ambient * ambientLight;

			for( int i=0; i<MAX_LIGHTS; i++ )
			{
				if( lights[i] != NULL )
				{
					lights[i]->FrontLightProduct_ambient  = material->ambient  * lights[i]->ambient;
					lights[i]->FrontLightProduct_specular = material->specular * lights[i]->specular;
				}
			}
		}
        */
/* // is this causing occasional crashes at startup?
		if( !node.impostor.empty() )
		{
			// TODO draw impostor quad
		} 
		else */ if( node.vbo ) {
			if( node.videomem_invalidated )
			{
                DEBUG_WARNING_ONCE("Found invalid tile - need to generate first!");

				generate2( &node );

				node.videomem_invalidated = false;
				node.hostmem_invalidated = true;

                DEBUG_RUN_ONCE(
                    skipRenderingThisNodeThisTime = true;
                );
			}

            if( !skipRenderingThisNodeThisTime ) { // TEMP
                drawMesh( node );
            } else {
                DEBUG_TRACE("Skipping drawMesh (might still be generating the node)");
            }
		}

		for( int i=0; i<node.children.size(); i++ )
		{
			draw( *node.children[i] );
		}

		ModelViewMatrix = ModelViewMatrixStack.pop();
	}


	void drawMesh( Node& node )
	{		
		DEBUG_ASSERT( node.vbo );
        DEBUG_ASSERT( node.ibo );

        // NOTE: In order to avoid syncronization CPU/GPU, invalidation is not managed properly.
        // So possible invalid nodes can't be skipped here (as invalidation is not reliable).
/*
        // TEMP Seeking for the reason of the crash
        if( node.videomem_invalidated ) {
            DEBUG("node.videomem_invalidated\n");
            return;
        } */
        if( node.vbo == NULL ) {
            DEBUG_WARNING("node.vbo == NULL");
            return;
        }
        if( node.ibo == NULL ) {
            DEBUG_WARNING("node.ibo == NULL");
            return;
        } /*
        if( node.vbo->videomem_invalidated ) {
            printf("node.vbo->videomem_invalidated\n");
            return;
        }
        if( node.ibo->videomem_invalidated ) {
            printf("node.ibo->videomem_invalidated\n");
            return;
        } */

		Node& tile = node;
		Node& clipmap = *tile.parent;
		Node& heightmap = *clipmap.parent;

		vec2 tileOffset = heightmap.transform.position.xy + clipmap.transform.position.xy + tile.transform.position.xy;

		//float tileSize = clipmap.tileSize;
		//vec2 location = vec2(heightmap.location.x, heightmap.location.y); // TODO handle the loss of data
		//vec2 world_tile_offset = -location + (clipmap.transform.position.xy + tile.transform.position.xy) * tileSize;

		VertexBuffer& vbo = *node.vbo;
		IndexBuffer&  ibo = *node.ibo;


		//
		// select effective lod
		//

		int lod = 0;

		//
		// Binding program
		//
        GL::GLSL::bind(programRenderTerrain);

		//
		// Binding uniforms
		// 
		#define NUM_CAUSTICS  32
 		float scale = node.vbo->mixmaps.scale.s;

		for( int i=0; i<Controls::CONTROLCOUNT; i++ ) 
		{
			int size = controls.literals[i].size() <=1 ? 2 : controls.literals[i].size();
			GL::GLSL::set( programRenderTerrain, controls.literals[i][0], controls.values[Controls::Bindings[i]] % size ); 
		}
		
		//
		// Binding attributes
		//
		GL::GLSL::set( programRenderTerrain, "quartetsTU",  0); GL::Texturing::bind( 0, vbo.quartets );
		GL::GLSL::set( programRenderTerrain, "gradientsTU", 1); GL::Texturing::bind( 1, vbo.gradients );
		GL::GLSL::set( programRenderTerrain, "colorsTU",    2); GL::Texturing::bind( 2, vbo.colors ); 
		GL::GLSL::set( programRenderTerrain, "mixmapsTU",   3); GL::Texturing::bind( 3, vbo.mixmaps ); 
		GL::GLSL::set( programRenderTerrain, "detailsTU",   4); GL::Texturing::bind( 4, details );
		GL::GLSL::set( programRenderTerrain, "detailsDxTU", 5); GL::Texturing::bind( 5, detailsDx );
		GL::GLSL::set( programRenderTerrain, "detailsDyTU", 6); GL::Texturing::bind( 6, detailsDy );

//#define NAME(var) #var
      //GL::GLSL::set( programRenderTerrain, NAME(detailsDy) "TU", 6);
//        programRenderTerrain->set(NAME(detailsDy) "TU", detailsDy, 6);
        // SET( program, detailsDy, 6)

		// LOD transitions blending
		GL::GLSL::set( programRenderTerrain, "tileOffset", tileOffset );
		GL::GLSL::set( programRenderTerrain, "kernelSize", float(CLIPMAP_WINDOW/2) );
		GL::GLSL::set( programRenderTerrain, "scale", scale );

		// transformation
		GL::GLSL::set( programRenderTerrain, "ModelViewProjectionMatrix",	ModelViewProjectionMatrix );
		GL::GLSL::set( programRenderTerrain, "ModelViewMatrix",				ModelViewMatrix );
		GL::GLSL::set( programRenderTerrain, "NormalMatrix",				NormalMatrix );

		// lighting
		GL::GLSL::set( programRenderTerrain, "Light0_position",				lights[0]->position );

		// light scattering 
		GL::GLSL::set( programRenderTerrain, "viewport", getViewport() );
		GL::GLSL::set( programRenderTerrain, "InverseRotationProjection", inverseRotationMatrix * inverseProjection(ProjectionMatrix) );
        
        //long visibileDistance = ((Heightmap)heightmap).visibility();
        //float visibileDistance = 4*(1<<(heightmap.children.size()));
        float visibileDistance = 2*(1<<(heightmap.children.size()));
        //DEBUG_TRACE_ONCE(visibileDistance);
        //DEBUG_TRACE_ONCE(heightmap.children.size());
        GL::GLSL::set( programRenderTerrain, "visibileDistance", visibileDistance );
		//GL::GLSL::set( programRenderTerrain, "AbsoluteTime",	float(Timer::absoluteTime()) );

        const int HeightBlendView = 2;
        if( controls.values[Controls::Bindings[Controls::DEBUGMODE]] == HeightBlendView ) {
            GL::GLSL::set(programRenderTerrain, "defaultColorR", vec4(1.0, 0.0, 0.0, 0.0));
            GL::GLSL::set(programRenderTerrain, "defaultColorG", vec4(0.0, 1.0, 0.0, 0.0));
            GL::GLSL::set(programRenderTerrain, "defaultColorB", vec4(0.0, 0.0, 1.0, 0.0));
            GL::GLSL::set(programRenderTerrain, "defaultColorA", vec4(0.5, 0.5, 0.5, 0.0));
        } else {
            GL::GLSL::set(programRenderTerrain, "defaultColorR", vec4(0.36, 0.30, 0.26, 0.0));  // Light stone
            GL::GLSL::set(programRenderTerrain, "defaultColorG", vec4(0.28, 0.24, 0.20, 0.0));  // Pebbles
            GL::GLSL::set(programRenderTerrain, "defaultColorB", vec4(0.34, 0.26, 0.22, 0.0));  // Brownish stone
            GL::GLSL::set(programRenderTerrain, "defaultColorA", vec4(0.50, 0.38, 0.30, 0.0));  // Dirty sand
        }





		//
		// Drawing
		//
			
		GL::VBO::bind( vbo );
	    GL::VBO::bind( ibo ); // Needs an index buffer bound

		glPatchParameteri(GL_PATCH_VERTICES, 4);		
		glDrawElements(
			GL_PATCHES, 
			ibo.lods[lod].count, // It is 1024 for tiles 17x17
			GL_UNSIGNED_SHORT, 
			0
		); 

		//
		// unbind 
		// 
		GL::Texturing::unbind(3);
		GL::Texturing::unbind(4);
		GL::Texturing::unbind(5);
		GL::GLSL::unbind();

		GL::VBO::unbind( GL_ARRAY_BUFFER );
		GL::VBO::unbind( GL_ELEMENT_ARRAY_BUFFER );

		//
		// vertex counter
		//
		this->vertices += ibo.lods[lod].count;
	}


	void generate2( Node* tile )
	{
		//glFlush();
		//glFinish();
		VertexBuffer& vbo = *tile->vbo;

        GL::GLSL::bind( programGenerateTerrain );

        GL::GLSL::set( programGenerateTerrain, "offset",	tile->object_id.xy );
		GL::GLSL::set( programGenerateTerrain, "size",		tile->object_id.z );

        // creation of texture, uploading and binding to texture unit
		GL::Texturing::bind( 0, vbo.quartets );
		GL::Texturing::bind( 1, vbo.gradients );
		GL::Texturing::bind( 2, vbo.colors ); 
		GL::Texturing::bind( 3, vbo.mixmaps ); 
		glBindImageTexture(0, vbo.quartets.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);  
		glBindImageTexture(1, vbo.gradients.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
		glBindImageTexture(2, vbo.colors.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
		glBindImageTexture(3, vbo.mixmaps.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

        // binding texture unit
	    //GL::GLSL::set(programGenerateTerrain, "detailsTU",   4); glBindImageTexture(4, details.id,	    0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);  
		//GL::GLSL::set(programGenerateTerrain, "detailsDxTU", 5); glBindImageTexture(5, detailsDx[0].id, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);
		//GL::GLSL::set(programGenerateTerrain, "detailsDyTU", 6); glBindImageTexture(6, detailsDy[0].id, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);
        //GL::Texturing::bind( 4, details );
        //glBindImageTexture(4, details.id,	0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32F);  

		// execute
		glDispatchCompute( vbo.quartets.size.x/2+1, vbo.quartets.size.y/2+1, 1 );
		//glFlush();
		//glFinish();
	}

	void drawSky()
	{
		GL::GLSL::bind( programRenderSky );
        for(int i = Controls::GAMMA; i <= Controls::VIGNETTING; i++) {
			int size = controls.literals[i].size() <= 1 ? 2 : controls.literals[i].size();
			GL::GLSL::set( programRenderSky, controls.literals[i][0], controls.values[Controls::Bindings[i]] % size ); 
		}
		GL::GLSL::set( programRenderSky, "viewport",	                getViewport() );
		GL::GLSL::set( programRenderSky, "InverseRotationProjection",   inverseRotationMatrix * inverseProjection(ProjectionMatrix) );
		GL::GLSL::set( programRenderSky, "Light0_position",		        lights[0]->position );
		GL::GLSL::set( programRenderSky, "AbsoluteTime",	            float(Timer::absoluteTime()) );
		glDrawArrays(GL_POINTS, 0, 1);
	}
/*
    void generateDetailGradient() 
    {
        GL::GLSL::bind(programGenerateGradientMap);

        // Creation of texture, uploading and binding to texture unit
		GL::Texturing::bind(0, details);
		GL::Texturing::bind(1, detailsDx[0]);
		GL::Texturing::bind(2, detailsDy[0]); 

		// Binding to image unit
		GL::GLSL::set(programGenerateGradientMap, "detailsTU",   0); glBindImageTexture(0, details.id,	 0, GL_FALSE, 0, GL_READ_ONLY,  GL_RGBA32F);  
		GL::GLSL::set(programGenerateGradientMap, "detailsDxIU", 1); glBindImageTexture(1, detailsDx[0].id, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
		GL::GLSL::set(programGenerateGradientMap, "detailsDyIU", 2); glBindImageTexture(2, detailsDy[0].id, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

        // Execute
	    glDispatchCompute( details.size.x, details.size.y, 1 );
    }
*/


};
