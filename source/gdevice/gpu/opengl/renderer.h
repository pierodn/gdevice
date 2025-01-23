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


#define MAX_LIGHTS 4 // 32


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

	Light light;

    Program programGenerateTerrain;
    Program programGenerateGradientMap;
	Program programRenderTerrain;
	Program programRenderSky;

	Renderer()
	{
		ModelViewMatrixStack.allocate(2000,100);
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
/*
		if(glEnableVertexAttribArray) 
        {
			glGetIntegerv( GL_MAX_VERTEX_ATTRIBS, &maxVertexAttributes );
			DEBUG_PRINT( TAB32 ": %i\n", "Max vertex attributes", maxVertexAttributes );
		}
*/
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

	bool isNonPowerOfTwoTexturesAvailable()
	{
		return GL::Texturing::textureNonPowerOfTwoAvailable();
	}



	void setViewport( vec2& viewport )
	{
		glViewport( 0, 0, viewport.x, viewport.y );
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
	}

    void setModelViewMatrix( mat4& mat)
    {
        ModelViewMatrix = mat4(1);
    }

    // buffer logic
	void clearBuffer( vec4& color = vec4(0) )
	{
		glClearColor( color.r, color.g, color.b , color.a );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	}

    // uniform logic
    void setLight(Light& light)
    {
        this->light = light;
    }




int verticesCount;

int traverse( Node& node, bool skipRendering = false )
{
    verticesCount = 0;
    

    traverse( node, verticesCount, skipRendering );

    return verticesCount;
}

// traverse(node, Transform::Update, Imposter::Update, Geometry::Update, Drawable::Draw ... )
void traverse( Node& node, int& verticesCount, bool skipRendering = false )
{
    Transform* pTransform = dynamic_cast<Transform*>(&node);
    if (pTransform) 
    {
	    ModelViewMatrixStack.push( ModelViewMatrix );
	    ModelViewMatrix *= TransformationMatrix( pTransform->transform );
	    ModelViewProjectionMatrix = ProjectionMatrix * ModelViewMatrix;
	    NormalMatrix = inverseTranspose( mat3(ModelViewMatrix) );
    }

    Impostor* pImpostor = dynamic_cast<Impostor*>(&node);
    Geometry* pGeometry = dynamic_cast<Geometry*>(&node);
    

    if( pImpostor && !pImpostor->impostorTexture.empty() )
	{
        // TODO AND is the impostor still valid?
		// TODO draw impostor quad: pImpostor->Render();
	}  
    else if(pGeometry && pGeometry->ibo)
    {
        Updatable* pUpdatable = dynamic_cast<Updatable*>(&node);
        if( pUpdatable && pUpdatable->videomem_invalidated )
		{
            pUpdatable->Update(*this);
			//Update( &node );

			pUpdatable->videomem_invalidated = false;
			pUpdatable->hostmem_invalidated = true;

            if( pImpostor && !pImpostor->impostorTexture.empty() )
            {
                // TODO pImpostor->Update()
            }
		}

        if( !skipRendering)
        {
            Renderable* pRenderable = dynamic_cast<Renderable*>(&node);
            if(pRenderable) 
            {
                RenderTarget target;
                pRenderable->Render(*this, target);
            }
        }
	}

    Parent* pParent = dynamic_cast<Parent*>(&node);
    if( pParent )
    {
	    for( int i=0; i<pParent->children.size(); i++ )
	    {
            Node* childAsNode = dynamic_cast<Node*>(pParent->children[i]);
            if(childAsNode) 
            {
		        traverse( *childAsNode );
            }
	    }
    }

    if (pTransform) 
    {
	    ModelViewMatrix = ModelViewMatrixStack.pop();
    }
}


// TODO 
void RenderTerrainTile( VertexBuffer& vbo, IndexBuffer& ibo, vec2 tileOffset, float visibileDistance  )
{		
    GL::GLSL::bind(programRenderTerrain);

	for( int i=0; i<Controls::CONTROLCOUNT; i++ ) 
	{
		int size = controls.literals[i].size() <=1 ? 2 : controls.literals[i].size();
		GL::GLSL::set( programRenderTerrain, controls.literals[i][0], controls.values[Controls::Bindings[i]] % size ); 
	}
	
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
	GL::GLSL::set( programRenderTerrain, "scale", vbo.mixmaps.scale.s );

	// transformation
	GL::GLSL::set( programRenderTerrain, "ModelViewProjectionMatrix",	ModelViewProjectionMatrix );
	GL::GLSL::set( programRenderTerrain, "ModelViewMatrix",				ModelViewMatrix );
	GL::GLSL::set( programRenderTerrain, "NormalMatrix",				NormalMatrix );

    // tessellation
    GL::GLSL::set( programRenderTerrain, "tessellationRange", 0.20f );
    GL::GLSL::set( programRenderTerrain, "tessellationFactor", 0.25f );
    GL::GLSL::set( programRenderTerrain, "tessellationDisplacement", 0.012f );

	// lighting
	GL::GLSL::set( programRenderTerrain, "Light0_position",				light.position );

	// light scattering 
	GL::GLSL::set( programRenderTerrain, "viewport", getViewport() );
	GL::GLSL::set( programRenderTerrain, "InverseRotationProjection", inverseRotationMatrix * inverseProjection(ProjectionMatrix) );
    
    GL::GLSL::set( programRenderTerrain, "visibileDistance", visibileDistance );
	GL::GLSL::set( programRenderTerrain, "AbsoluteTime",	float(Timer::absoluteTime()) );

    const int HeightBlendView = 2;
    if( controls.values[Controls::Bindings[Controls::DEBUGMODE]] == HeightBlendView )
    {
        GL::GLSL::set(programRenderTerrain, "defaultColorR", vec4(1.0, 0.0, 0.0, 0.0));
        GL::GLSL::set(programRenderTerrain, "defaultColorG", vec4(0.0, 1.0, 0.0, 0.0));
        GL::GLSL::set(programRenderTerrain, "defaultColorB", vec4(0.0, 0.0, 1.0, 0.0));
        GL::GLSL::set(programRenderTerrain, "defaultColorA", vec4(1.0, 0.0, 0.0, 0.0));
    } else {
        GL::GLSL::set(programRenderTerrain, "defaultColorR", vec4(0.36, 0.30, 0.26, 0.0));  // Light stone
        GL::GLSL::set(programRenderTerrain, "defaultColorG", vec4(0.28, 0.24, 0.20, 0.0));  // Pebbles
        GL::GLSL::set(programRenderTerrain, "defaultColorB", vec4(0.34, 0.26, 0.22, 0.0));  // Brownish stone
        GL::GLSL::set(programRenderTerrain, "defaultColorA", vec4(0.50, 0.38, 0.30, 0.0));  // Dirty sand
    }
		
	GL::VBO::bind( vbo );
    GL::VBO::bind( ibo ); // Needs an index buffer bound (for now)

    int lod = 0; // TODO

	glPatchParameteri(GL_PATCH_VERTICES, 4);		
	glDrawElements(
		GL_PATCHES, 
		ibo.lods[lod].count, // It is 1024 for tiles 17x17
		GL_UNSIGNED_SHORT, 
		0
	); 

	GL::Texturing::unbind(3);
	GL::Texturing::unbind(4);
	GL::Texturing::unbind(5);
	GL::GLSL::unbind();

	GL::VBO::unbind( GL_ARRAY_BUFFER );
	GL::VBO::unbind( GL_ELEMENT_ARRAY_BUFFER );

	verticesCount += ibo.lods[lod].count;
}
 


void GenerateTerrainTile(const vec4& tileID, VertexBuffer& vbo)
{
    GL::GLSL::bind( programGenerateTerrain );

    GL::GLSL::set( programGenerateTerrain, "offset",	tileID.xy );
    GL::GLSL::set( programGenerateTerrain, "size",		tileID.z );

    GL::Texturing::bind( 0, vbo.quartets );
    GL::Texturing::bind( 1, vbo.gradients );
    GL::Texturing::bind( 2, vbo.colors ); 
    GL::Texturing::bind( 3, vbo.mixmaps ); 

    glBindImageTexture(0, vbo.quartets.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);  
    glBindImageTexture(1, vbo.gradients.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
    glBindImageTexture(2, vbo.colors.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
    glBindImageTexture(3, vbo.mixmaps.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

    glDispatchCompute( vbo.quartets.size.x/2+1, vbo.quartets.size.y/2+1, 1 );
}

void drawSky()
{
	GL::GLSL::bind( programRenderSky );

    for(int i = Controls::GAMMA; i <= Controls::VIGNETTING; i++)
    {
		int size = controls.literals[i].size() <= 1 ? 2 : controls.literals[i].size();
		GL::GLSL::set( programRenderSky, controls.literals[i][0], controls.values[Controls::Bindings[i]] % size ); 
	}
	GL::GLSL::set( programRenderSky, "viewport",	                getViewport() );
	GL::GLSL::set( programRenderSky, "InverseRotationProjection",   inverseRotationMatrix * inverseProjection(ProjectionMatrix) );
	GL::GLSL::set( programRenderSky, "Light0_position",		        light.position );
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
