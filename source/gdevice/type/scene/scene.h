#pragma once

#include "type/scene/node.h"
#include "gpu/opengl/renderer.h"

struct Scene : Node, Transform, Parent
{
/*
	mat4 ProjectionMatrix;
	mat4 ModelViewMatrix;
	mat3 NormalMatrix;s
	mat4 ModelViewProjectionMatrix;
	mat4 inverseRotationMatrix;
	Array<mat4> ModelViewMatrixStack;    

    Scene()
	{
		ModelViewMatrixStack.allocate(2000,100);
	}
*/
    
/*
    // TODO int traverse(node, Transform::Update, Imposter::Update, Geometry::Update, Drawable::Draw ... )
    int vertices;

    void traverse( Node& node, bool skipRendering = false )
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
*/
};