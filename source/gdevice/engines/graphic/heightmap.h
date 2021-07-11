#pragma once

#include "platforms/thread.h"
#include "platforms/mutex.h"

#include "engines/graphic/clipmap.h"
#include "engines/graphic/skybox.h"
#include "engines/graphic/renderer.h"



// TODO: Rename as Terrain (?)
struct Heightmap : public Node
{
private:
	int   tileRes;
	IndexBuffer ibo;

public:
	int getTileResolution()
	{
		return tileRes;
	}

	void setTileResolution(int resolution)
	{
        float ex = log2(float(resolution-1));
        DEBUG_ASSERT(ex - int(ex) < 1E-9); // Must be 2^N+1
		this->tileRes = (1 << int(ex)) + 1;	
		ibo.init(this->tileRes);
	}

	Heightmap()
	{
		vbo = NULL;
		tileRes = 17;

		transform.position.xy += CLIPMAP_ODDITY ? 0.0f : -1.0f;
		transform.scale    = vec3(1);
	}

	~Heightmap()
	{
		setLODs(0);		
	}

	void setLODs( int lods )
	{
		int d = lods - children.size();
		for(; d>0; d-- ) children.push( new Clipmap( this, children.size(), tileRes, &ibo/*tstrips*/ ) );
		for(; d<0; d++ ) { delete children[children.size()-1]; children.pop(); }
	}

    long visibility()
	{
		return 4*(1<<(children.size()));
	}

	// TODO: Returns z from surface
	float moveAt(dmat3 camera) 
    {
        int coarsestLODToBeUpdated = -1; // TODO: Use.
        for( int lod = 0; lod<children.size(); lod++ ) {
            Clipmap* clipmap = (Clipmap*)children[lod];
			if( clipmap->move(-camera.position.xy) ) {
				coarsestLODToBeUpdated = lod;
			}
		}

        transform.rotation = vec3( -camera.rotation.x, -camera.rotation.y, -camera.rotation.z );
	    transform.position.z = -camera.position.z;

        //bindTileWithItsCoarserTile();

        // TODO: Return height detected (?)
		return 0.0; 
    }

    void generateInvalidatedTiles(Renderer& renderer) 
    {
        DEBUG_PROFILE_ONCE(5.0);
        for(int lod=children.size()-1; lod>=0; lod--) {
            Clipmap* clipmap = (Clipmap*)children[lod];
	        for(int i=0; i<CLIPMAP_SIZE; i++)
	        for(int j=0; j<CLIPMAP_SIZE; j++) {
                Node& node = clipmap->getTile(i,j); 

			    if( node.videomem_invalidated ) {
				    renderer.generate2( &node );
				    node.videomem_invalidated = false;
				    node.hostmem_invalidated = true;
                }
	        }
        }
	}
/*
    // NOTE: Might be invariant
    void bindTileWithItsCoarserTile() 
    {
        for(int lod=0; lod<children.size()-2; lod++) {

            Clipmap* clipmap = (Clipmap*)children[lod];
            Clipmap* coarserClipmap = (Clipmap*)children[lod+1];

            for(int i=0; i<CLIPMAP_SIZE; i++)
            for(int j=0; j<CLIPMAP_SIZE; j++) {

                ivec2 ij = ivec2(i,j)/2 + CLIPMAP_SIZE/2 - 1;
                ivec2 qu = ivec2(i % 2, j % 2);
                Node& coarserNode = coarserClipmap->getTile(ij.x, ij.y);

                Node& node = clipmap->getTile(i, j);
                node.coarser = &coarserNode;
                node.coarserQuadrant = vec2(qu.x, qu.y);
            }
        }  
    }
*/


};