#pragma once

#include "type/scene/node.h"
#include "gpu/opengl/renderer.h"

struct Scene : Node, Transform, Parent
{
/*
    Controls controls;

    vec2 viewport;
	mat4 ProjectionMatrix;

    mat4 inverseRotationMatrix;

	mat4 ModelViewMatrix;
	Array<mat4> ModelViewMatrixStack;   

    Light light;

    Program programGenerateTerrain;
    Program programGenerateGradientMap;
	Program programRenderTerrain;
	Program programRenderSky;
*/
    Renderer* renderer;

    Scene()
	{
		//ModelViewMatrixStack.allocate(2000,100);
	}

    void Initialize(Renderer* renderer)
    {
        this->renderer = renderer;
        renderer->initialize();
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
	        renderer->ModelViewMatrixStack.push( renderer->ModelViewMatrix );
	        renderer->ModelViewMatrix *= TransformationMatrix( pTransform->transform );
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
                pUpdatable->Update(*renderer);
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
                    pRenderable->Render(*renderer, target);
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
	        renderer->ModelViewMatrix = renderer->ModelViewMatrixStack.pop();
        }

        return vertexCount;
    }


};