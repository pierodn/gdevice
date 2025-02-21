#pragma once

#include "gpu/controls.h"

#include "type/scene/node.h"

//#include "gpu/opengl/renderer.h"
#include "gpu/opengl/gl.h"
#include "__temp/light.h"

#include "application/assets/worlds/planet1/render_sky.glsl.h"

struct NodeTransform
{
	mat4 ProjectionMatrix;
    mat4 inverseRotationMatrix; // TODO inverseCameraRotationMatrix;
	mat4 ModelViewMatrix;
	Array<mat4> ModelViewMatrixStack;   

    NodeTransform()
    {
		ModelViewMatrixStack.allocate(2000,100);
    }
};

struct SceneState
{
    Light light;
};

struct Scene : Node, Transform, Parent
{
    NodeTransform nodeTransform;    
    SceneState state;

    Program programRenderSky;

    Renderer* renderer; // TODO remove

    Scene()
	{
	}

    void Initialize()
    {
        GL::GLSL::build( programRenderSky, render_sky_glsl );
    }

    


    int Traverse( bool skipRendering = false )
    {
        int vertexCount = 0;
        return Traverse( *this, vertexCount, skipRendering );
    }

    int Traverse( Node& node, int& vertexCount, bool skipRendering )
    {
        Transform* pTransform = dynamic_cast<Transform*>(&node);
        if (pTransform) 
        {
            nodeTransform.ModelViewMatrixStack.push( nodeTransform.ModelViewMatrix );
	        nodeTransform.ModelViewMatrix *= TransformationMatrix( pTransform->transform );
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
                pUpdatable->Update();
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
                    pRenderable->Render(*renderer, nodeTransform, state);
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
		            Traverse( *childAsNode, vertexCount, skipRendering );
                }
	        }
        }

        if (pTransform) 
        {
            nodeTransform.ModelViewMatrix = nodeTransform.ModelViewMatrixStack.pop();
        }

        return vertexCount;
    }



    // TODO rename as drawGI
    void drawSky()
    {
	    GL::GLSL::bind( programRenderSky );

        for(int i = Controls::GAMMA; i <= Controls::VIGNETTING; i++)
        {
		    int size = Controls::GetInstance().literals[i].size() <= 1 ? 2 : Controls::GetInstance().literals[i].size();
		    GL::GLSL::set( programRenderSky, Controls::GetInstance().literals[i][0], Controls::GetInstance().values[Controls::Bindings[i]] % size ); 
	    }

        vec2 viewport = GL::GetViewport();
	    GL::GLSL::set( programRenderSky, "viewport",	                viewport  );
	    GL::GLSL::set( programRenderSky, "InverseRotationProjection",   nodeTransform.inverseRotationMatrix * inverseProjection(nodeTransform.ProjectionMatrix) );
	    GL::GLSL::set( programRenderSky, "Light0_position",		        state.light.position );
	    GL::GLSL::set( programRenderSky, "AbsoluteTime",	            float(Timer::absoluteTime()) );
    	
        glDrawArrays(GL_POINTS, 0, 1);
    }


};