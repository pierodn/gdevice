#pragma once

#include <glpp/program.h>
#include <glpp/glpp.h>
#include <type/scene/oop/node.h>
#include "textures.h"


struct NodeTransform;

struct Tile : public Node, Transform, Child, Geometry, Updatable, Renderable
{
    vec4 tileID;

    Program* generator;
    Program* renderer;

    void Update()
    {
        DEBUG_ASSERT(vbo);
        DEBUG_ASSERT(generator);

        generator->Use();
        generator->SetUniform("offset", tileID.xy);
        generator->SetUniform("size", tileID.z);

        // TODO generator->SetOutputBuffer()
        GL::Texturing::bind( 0, vbo->quartets );
        GL::Texturing::bind( 1, vbo->gradients );
        GL::Texturing::bind( 2, vbo->colors ); 
        GL::Texturing::bind( 3, vbo->mixmaps ); 
        glBindImageTexture(0, vbo->quartets.id,	    0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);  
        glBindImageTexture(1, vbo->gradients.id,	0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
        glBindImageTexture(2, vbo->colors.id,	    0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
        glBindImageTexture(3, vbo->mixmaps.id,	    0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

        // TODO generator->Run(vbo->quartets.size/2.0 + 1.0);
        glDispatchCompute( vbo->quartets.size.x/2 + 1, vbo->quartets.size.y/2 + 1, 1 );
    }

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

        VertexBuffer& vbo = *this->vbo;
        IndexBuffer& ibo = *this->ibo;

        renderer->Use();

	    for( int i=0; i<Controls::CONTROLCOUNT; i++ ) 
	    {
		    int size = Controls::GetInstance().literals[i].size() <=1 ? 2 : Controls::GetInstance().literals[i].size();
		    renderer->SetUniform(Controls::GetInstance().literals[i][0], Controls::GetInstance().values[Controls::Bindings[i]] % size); 
	    }
    	
	    renderer->SetUniform("quartetsTU",  0); GL::Texturing::bind( 0, vbo.quartets );
	    renderer->SetUniform("gradientsTU", 1); GL::Texturing::bind( 1, vbo.gradients );
	    renderer->SetUniform("colorsTU",    2); GL::Texturing::bind( 2, vbo.colors ); 
	    renderer->SetUniform("mixmapsTU",   3); GL::Texturing::bind( 3, vbo.mixmaps ); 
	    renderer->SetUniform("detailsTU",   4); GL::Texturing::bind( 4, details );
	    renderer->SetUniform("detailsDxTU", 5); GL::Texturing::bind( 5, detailsDx );
	    renderer->SetUniform("detailsDyTU", 6); GL::Texturing::bind( 6, detailsDy );

    //#define NAME(var) #var
      //renderer->Set(NAME(detailsDy) "TU", 6);
    //        programRenderTerrain->set(NAME(detailsDy) "TU", detailsDy, 6);
        // SET( program, detailsDy, 6)

	    // LOD transitions blending
        renderer->SetUniform("tileOffset", tileOffset);
        renderer->SetUniform("kernelSize", float(CLIPMAP_WINDOW/2));
	    renderer->SetUniform("scale", vbo.mixmaps.scale.s);

	    // transformation
        mat4 ModelViewProjectionMatrix = nodeState.ProjectionMatrix * nodeState.ModelViewMatrix;
        mat3 NormalMatrix = inverseTranspose(mat3(nodeState.ModelViewMatrix));
	    renderer->SetUniform("ModelViewProjectionMatrix", ModelViewProjectionMatrix);
	    renderer->SetUniform("ModelViewMatrix", nodeState.ModelViewMatrix);
	    renderer->SetUniform("NormalMatrix", NormalMatrix);

        // tessellation
        renderer->SetUniform("tessellationKernelStart", 0.04f ); // 0.03
        renderer->SetUniform("tessellationKernelRange", 0.16f ); // 0.10
        renderer->SetUniform("tessellationVanishPower", 0.50f ); // less is more
        renderer->SetUniform("tessellationMaxLevel", 0.25f );
        renderer->SetUniform("tessellationDisplacement", 0.010f );
        renderer->SetUniform("povZ", pHeightmap->transform.position.z );

	    // lighting
	    renderer->SetUniform("Light0_position",	vec4(sceneState.sun, 0.0) );

	    // light scattering 
        vec2 viewport = GL::GetViewport();
	    renderer->SetUniform("viewport", viewport  );
	    renderer->SetUniform("InverseRotationProjection", nodeState.inverseRotationMatrix * inverseProjection(nodeState.ProjectionMatrix) );
        
        renderer->SetUniform("visibileDistance", visibileDistance );
	    renderer->SetUniform("AbsoluteTime",	float(Timer::absoluteTime()) );

        const int HeightBlendView = 2;
        if( Controls::GetInstance().values[Controls::Bindings[Controls::DEBUGMODE]] == HeightBlendView )
        {
            renderer->SetUniform("defaultColorR", vec4(1.0, 0.0, 0.0, 0.0));
            renderer->SetUniform("defaultColorG", vec4(0.0, 1.0, 0.0, 0.0));
            renderer->SetUniform("defaultColorB", vec4(0.0, 0.0, 1.0, 0.0));
            renderer->SetUniform("defaultColorA", vec4(1.0, 0.0, 0.0, 0.0));
        } else {
            renderer->SetUniform("defaultColorR", vec4(0.36, 0.30, 0.26, 0.0));  // Light stone
            renderer->SetUniform("defaultColorG", vec4(0.28, 0.24, 0.20, 0.0));  // Pebbles
            renderer->SetUniform("defaultColorB", vec4(0.34, 0.26, 0.22, 0.0));  // Brownish stone
            renderer->SetUniform("defaultColorA", vec4(0.50, 0.38, 0.30, 0.0));  // Dirty sand
        }
    		
	    GL::VBO::bind( vbo );
        GL::VBO::bind( ibo ); // Needs an index buffer bound (for now)

        int lod = 0; // TODO

        // TODO renderer->Run()
	    glPatchParameteri(GL_PATCH_VERTICES, 4);		
	    glDrawElements(
		    GL_PATCHES, 
		    64*64, //ibo.lods[lod].count, // It is 1024 for tiles 17x17
		    GL_UNSIGNED_SHORT, 
		    0
	    ); 

//        glDrawArrays(GL_PATCHES, 0, 64*64);

	    GL::Texturing::unbind(3);
	    GL::Texturing::unbind(4);
	    GL::Texturing::unbind(5);
	    //GL::UseProgram(0); // unbind
            glUseProgram(0);

	    GL::VBO::unbind( GL_ARRAY_BUFFER );
	    GL::VBO::unbind( GL_ELEMENT_ARRAY_BUFFER );

	    // TODO verticesCount += ibo.lods[lod].count;
    }
 




};

