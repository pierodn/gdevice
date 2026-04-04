#pragma once

#include <os/thread.h>
#include <os/mutex.h>
#include <gl/program.h>

#include <type/scene/oop/terrain/clipmap.h>

// TODO this should be a configuration to the generic Heightmap instance
#include "shaders/generate_terrain.glsl.h"
#include "shaders/render_terrain.glsl.h"


struct Heightmap : public Node, Parent, Child, Transform 
{
private:
	int   tileRes;
    dvec2 location;
    
    IndexBuffer ibo;

    Program generator;
    Program generateGradientMapProgram; // TODO
	Program renderer;

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
		InitializeTerrainDetails();

        //generator.Build(generateTerrainSource);
        generator.Build(generate_terrain_glsl);
        renderer.Build(render_terrain_glsl);

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
		for(; d>0; d-- ) children.push_back( new Clipmap(this, children.size(), tileRes, &ibo, &generator, &renderer) );
		for(; d<0; d++ ) delete children[children.size()-1];
	}

    long GetVisibilityDistance()
	{
		return 4*(1<<(children.size()));
	}

	// TODO: Returns z from surface
	float moveAt(dmat3 camera) 
    {
        int coarsestLODToBeUpdated = -1; // TODO: Use.
        for( unsigned int lod = 0; lod<children.size(); lod++ )
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


};