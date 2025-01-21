#pragma once

#include "type/scene/node.h"


struct Tile : public Node, Transform, Child, Geometry, Updatable, Renderable
{
    vec4 tileID;

    void Render(Renderer& renderer, RenderTarget& target)
    {
        // Ensure it's overriding 
        static_cast<void(Renderable::*)(Renderer& renderer, RenderTarget& target)>(&Tile::Render);

        DEBUG_ASSERT( vbo );
        DEBUG_ASSERT( ibo );
        DEBUG_ASSERT( parent );

        Transform* pClipmap = dynamic_cast<Transform*>(parent);
        DEBUG_ASSERT( pClipmap );

        Child* pClipmapAsChild = dynamic_cast<Child*>(pClipmap);
        DEBUG_ASSERT( pClipmapAsChild );
        DEBUG_ASSERT( pClipmapAsChild->parent );

        Transform* pHeightmap = dynamic_cast<Transform*>(pClipmapAsChild->parent);
        Parent* pHeightmapAsParent = dynamic_cast<Parent*>(pClipmapAsChild->parent);

        vec2 tileOffset = pHeightmap->transform.position.xy + pClipmap->transform.position.xy + transform.position.xy;

        float visibileDistance = 2*(1<<(pHeightmapAsParent->children.size()));


        renderer.RenderTerrainTile(*vbo, *ibo, tileOffset, visibileDistance);
    }

    void Update(Renderer& renderer) 
    {
        DEBUG_ASSERT(vbo);
        renderer.GenerateTerrainTile(tileID, *vbo);
    }
};
