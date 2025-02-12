#pragma once

#include "os/thread.h"
#include "os/mutex.h"
#include "type/scene/terrain/clipmap.h"
#include "gpu/renderer.h"
#include "gpu/program.h"

/*
    #include "application/assets/worlds/planet1/generate_terrain.glsl.h"
    #include "application/assets/worlds/planet1/render_terrain.glsl.h"
    #include "application/assets/worlds/planet1/render_sky.glsl.h"
*/



struct Heightmap : public Node, Parent, Child, Transform 
{
private:
	int   tileRes;
    dvec2 location;
    
    IndexBuffer ibo;

    // TODO
    Program generateProgram;
    Program generateGradientMapProgram; // TODO
	Program renderProgram;
	Program renderSkyProgram;



public:

	Heightmap()
	{
		tileRes = 17;

		transform.position.xy += CLIPMAP_ODDITY ? 0.0f : -1.0f;
		transform.scale    = vec3(1);
	}

	~Heightmap()
	{
		setLODs(0);		
	}

    bool Initialize()
    {
        GL::GLSL::build( generateProgram, generate_terrain_glsl);
        // TODO GL::GLSL::build( generateGradientMapProgram, generate_grandienmap_glsl );
        GL::GLSL::build( renderProgram, render_terrain_glsl );
		GL::GLSL::build( renderSkyProgram, render_sky_glsl );
        return true;
    }

	void setTileResolution(int resolution)
	{
        float ex = log2(float(resolution-1));
        DEBUG_ASSERT(ex - int(ex) < 1E-9); // Must be 2^N+1
		this->tileRes = (1 << int(ex)) + 1;	
		ibo.init(this->tileRes);
	}

    int getTileResolution()
	{
		return tileRes;
	}

	void setLODs( int lods )
	{
		int d = lods - children.size();
		for(; d>0; d-- ) children.push( new Clipmap( this, children.size(), tileRes, &ibo/*tstrips*/ ) );
		for(; d<0; d++ ) { delete children[children.size()-1]; children.pop(); }
	}

    long GetVisibilityDistance()
	{
		return 4*(1<<(children.size()));
	}

	// TODO: Returns z from surface
	float moveAt(dmat3 camera) 
    {
        int coarsestLODToBeUpdated = -1; // TODO: Use.
        for( int lod = 0; lod<children.size(); lod++ )
        {
            Clipmap* clipmap = (Clipmap*)children[lod];
			if( clipmap->move(-camera.position.xy) ) 
            {
				coarsestLODToBeUpdated = lod;
			}
		}

        transform.rotation = vec3( -camera.rotation.x, -camera.rotation.y, -camera.rotation.z );
	    transform.position.z = -camera.position.z;

        //bindTileWithItsCoarserTile();

        // TODO: Return height detected (?)
		return 0.0; 
    }



    void GenerateTile(const vec4& tileID, VertexBuffer& vbo)
    {
        GL::GLSL::bind( generateProgram );

        GL::GLSL::set( generateProgram, "offset",	tileID.xy );
        GL::GLSL::set( generateProgram, "size",		tileID.z );

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
/*
*/
/*
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
        mat4 ModelViewProjectionMatrix = ProjectionMatrix * ModelViewMatrix;
        mat3 NormalMatrix = inverseTranspose( mat3(ModelViewMatrix) );
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
	    GL::GLSL::set( programRenderTerrain, "viewport", viewport  );
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

	    // TODO verticesCount += ibo.lods[lod].count;
    }
*/




};