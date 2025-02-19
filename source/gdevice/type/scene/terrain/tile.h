#pragma once

#include "type/scene/node.h"

struct NodeTransform;

struct Tile : public Node, Transform, Child, Geometry, Updatable, Renderable
{
    vec4 tileID;

    Program* generator;
    Program* renderer;

    void Render(Renderer& renderer, NodeTransform& nodeTransform, SceneState& sceneState)
    {
        // Ensure it's overriding the virtual function of the base class.
        //static_cast<void(Renderable::*)(Renderer& renderer, RenderTarget& target)>(&Tile::Render);
    nodeTransform.inverseRotationMatrix;

        DEBUG_ASSERT( vbo );
        DEBUG_ASSERT( ibo );
        DEBUG_ASSERT( parent );

        Transform* pClipmap = dynamic_cast<Transform*>(parent);
        DEBUG_ASSERT( pClipmap );

        Child* pClipmapAsChild = dynamic_cast<Child*>(pClipmap);
        DEBUG_ASSERT( pClipmapAsChild );
        DEBUG_ASSERT( pClipmapAsChild->parent );

        // TODO Transform* pHeightmap = Cast(pClipmapAsChild->parent);
        Transform* pHeightmap = dynamic_cast<Transform*>(pClipmapAsChild->parent);
        Parent* pHeightmapAsParent = dynamic_cast<Parent*>(pClipmapAsChild->parent);

        vec2 tileOffset = pHeightmap->transform.position.xy + pClipmap->transform.position.xy + transform.position.xy;

        float visibileDistance = 2*(1<<(pHeightmapAsParent->children.size()));

        // TODO use the local renderer program instead of the Renderer 
        //renderer.RenderTerrainTile(*this->renderer, nodeTransform, *vbo, *ibo, tileOffset, visibileDistance);
        __RenderTerrainTile(renderer, *this->renderer, nodeTransform, sceneState, *vbo, *ibo, tileOffset, visibileDistance);
    }

void __RenderTerrainTile(Renderer& renderer, Program& programRenderTerrain, NodeTransform& nodeTransform, SceneState& sceneState, VertexBuffer& vbo, IndexBuffer& ibo, vec2 tileOffset, float visibileDistance  )
{		
nodeTransform.ModelViewMatrix;

    GL::GLSL::bind(programRenderTerrain);

	for( int i=0; i<Controls::CONTROLCOUNT; i++ ) 
	{
		int size = Controls::GetInstance().literals[i].size() <=1 ? 2 : Controls::GetInstance().literals[i].size();
		GL::GLSL::set( programRenderTerrain, Controls::GetInstance().literals[i][0], Controls::GetInstance().values[Controls::Bindings[i]] % size ); 
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
    mat4 ModelViewProjectionMatrix = nodeTransform.ProjectionMatrix * nodeTransform.ModelViewMatrix;
    mat3 NormalMatrix = inverseTranspose( mat3(nodeTransform.ModelViewMatrix) );
	GL::GLSL::set( programRenderTerrain, "ModelViewProjectionMatrix",	ModelViewProjectionMatrix );
	GL::GLSL::set( programRenderTerrain, "ModelViewMatrix",				nodeTransform.ModelViewMatrix );
	GL::GLSL::set( programRenderTerrain, "NormalMatrix",				NormalMatrix );

    // tessellation
    GL::GLSL::set( programRenderTerrain, "tessellationRange", 0.20f );
    GL::GLSL::set( programRenderTerrain, "tessellationFactor", 0.25f );
    GL::GLSL::set( programRenderTerrain, "tessellationDisplacement", 0.012f );

	// lighting
	GL::GLSL::set( programRenderTerrain, "Light0_position",	sceneState.light.position );

	// light scattering 
    vec2 viewport = renderer.GetViewport();
	GL::GLSL::set( programRenderTerrain, "viewport", viewport  );
	GL::GLSL::set( programRenderTerrain, "InverseRotationProjection", nodeTransform.inverseRotationMatrix * inverseProjection(nodeTransform.ProjectionMatrix) );
    
    GL::GLSL::set( programRenderTerrain, "visibileDistance", visibileDistance );
	GL::GLSL::set( programRenderTerrain, "AbsoluteTime",	float(Timer::absoluteTime()) );

    const int HeightBlendView = 2;
    if( Controls::GetInstance().values[Controls::Bindings[Controls::DEBUGMODE]] == HeightBlendView )
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
 


    void Update() 
    {
        DEBUG_ASSERT(vbo);
        DEBUG_ASSERT(generator);

        GL::GLSL::bind( *generator );
        GL::GLSL::set( *generator, "offset",	tileID.xy );
        GL::GLSL::set( *generator, "size", tileID.z );

        GL::Texturing::bind( 0, vbo->quartets );
        GL::Texturing::bind( 1, vbo->gradients );
        GL::Texturing::bind( 2, vbo->colors ); 
        GL::Texturing::bind( 3, vbo->mixmaps ); 

        glBindImageTexture(0, vbo->quartets.id,	    0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);  
        glBindImageTexture(1, vbo->gradients.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
        glBindImageTexture(2, vbo->colors.id,	    0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
        glBindImageTexture(3, vbo->mixmaps.id,	    0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

        glDispatchCompute( vbo->quartets.size.x/2 + 1, vbo->quartets.size.y/2 + 1, 1 );
    }

};

