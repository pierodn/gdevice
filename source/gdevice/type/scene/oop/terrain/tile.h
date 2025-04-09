#pragma once

#include "type/scene/oop/node.h"

#include "application/assets/worlds/planet1/textures.h"

struct NodeTransform;

struct Tile : public Node, Transform, Child, Geometry, Updatable, Renderable
{
    vec4 tileID;

    Program* generator;
    Program* renderer;

    void Render(NodeState& nodeState, SceneState& sceneState)
    {
        // Ensure it's overriding the virtual function of the base class.
        static_cast<void(Renderable::*)(NodeState&, SceneState&)>(&Tile::Render);

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

        Program& programRenderTerrain = *this->renderer;
        VertexBuffer& vbo = *this->vbo;
        IndexBuffer& ibo = *this->ibo;

        GL::GLSL::bind(programRenderTerrain);

	    for( int i=0; i<Controls::CONTROLCOUNT; i++ ) 
	    {
		    int size = Controls::GetInstance().literals[i].size() <=1 ? 2 : Controls::GetInstance().literals[i].size();
		    renderer->Set(Controls::GetInstance().literals[i][0], Controls::GetInstance().values[Controls::Bindings[i]] % size); 
	    }
    	
	    renderer->Set("quartetsTU",  0); GL::Texturing::bind( 0, vbo.quartets );
	    renderer->Set("gradientsTU", 1); GL::Texturing::bind( 1, vbo.gradients );
	    renderer->Set("colorsTU",    2); GL::Texturing::bind( 2, vbo.colors ); 
	    renderer->Set("mixmapsTU",   3); GL::Texturing::bind( 3, vbo.mixmaps ); 
	    renderer->Set("detailsTU",   4); GL::Texturing::bind( 4, details );
	    renderer->Set("detailsDxTU", 5); GL::Texturing::bind( 5, detailsDx );
	    renderer->Set("detailsDyTU", 6); GL::Texturing::bind( 6, detailsDy );

    //#define NAME(var) #var
      //renderer->Set(NAME(detailsDy) "TU", 6);
    //        programRenderTerrain->set(NAME(detailsDy) "TU", detailsDy, 6);
        // SET( program, detailsDy, 6)

	    // LOD transitions blending
        renderer->Set("tileOffset", tileOffset);
        renderer->Set("kernelSize", float(CLIPMAP_WINDOW/2));
	    renderer->Set("scale", vbo.mixmaps.scale.s);

	    // transformation
        mat4 ModelViewProjectionMatrix = nodeState.ProjectionMatrix * nodeState.ModelViewMatrix;
        mat3 NormalMatrix = inverseTranspose(mat3(nodeState.ModelViewMatrix));
	    renderer->Set("ModelViewProjectionMatrix", ModelViewProjectionMatrix);
	    renderer->Set("ModelViewMatrix", nodeState.ModelViewMatrix);
	    renderer->Set("NormalMatrix", NormalMatrix);

        // tessellation
        renderer->Set("tessellationKernelStart", 0.04f ); // 0.03
        renderer->Set("tessellationKernelRange", 0.16f ); // 0.10
        renderer->Set("tessellationVanishPower", 0.50f ); // less is more
        renderer->Set("tessellationMaxLevel", 0.25f );
        renderer->Set("tessellationDisplacement", 0.010f );
        renderer->Set("povZ", pHeightmap->transform.position.z );

	    // lighting
	    renderer->Set("Light0_position",	vec4(sceneState.sun, 0.0) );

	    // light scattering 
        vec2 viewport = GL::GetViewport();
	    renderer->Set("viewport", viewport  );
	    renderer->Set("InverseRotationProjection", nodeState.inverseRotationMatrix * inverseProjection(nodeState.ProjectionMatrix) );
        
        renderer->Set("visibileDistance", visibileDistance );
	    renderer->Set("AbsoluteTime",	float(Timer::absoluteTime()) );

        const int HeightBlendView = 2;
        if( Controls::GetInstance().values[Controls::Bindings[Controls::DEBUGMODE]] == HeightBlendView )
        {
            renderer->Set("defaultColorR", vec4(1.0, 0.0, 0.0, 0.0));
            renderer->Set("defaultColorG", vec4(0.0, 1.0, 0.0, 0.0));
            renderer->Set("defaultColorB", vec4(0.0, 0.0, 1.0, 0.0));
            renderer->Set("defaultColorA", vec4(1.0, 0.0, 0.0, 0.0));
        } else {
            renderer->Set("defaultColorR", vec4(0.36, 0.30, 0.26, 0.0));  // Light stone
            renderer->Set("defaultColorG", vec4(0.28, 0.24, 0.20, 0.0));  // Pebbles
            renderer->Set("defaultColorB", vec4(0.34, 0.26, 0.22, 0.0));  // Brownish stone
            renderer->Set("defaultColorA", vec4(0.50, 0.38, 0.30, 0.0));  // Dirty sand
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

