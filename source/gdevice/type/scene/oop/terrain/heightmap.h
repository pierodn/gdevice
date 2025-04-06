#pragma once

#include "os/thread.h"
#include "os/mutex.h"
#include "type/scene/oop/terrain/clipmap.h"
#include "gpu/program.h"

#include "application/assets/worlds/planet1/generate_terrain.glsl.h"
#include "application/assets/worlds/planet1/render_terrain.glsl.h"


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

        GL::GLSL::build( generator, generate_terrain_glsl);
        // TODO GL::GLSL::build( generateGradientMapProgram, generate_grandienmap_glsl );
        GL::GLSL::build( renderer, render_terrain_glsl );
		//GL::GLSL::build( renderSkyProgram, render_sky_glsl );

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
		for(; d>0; d-- ) children.push( new Clipmap(this, children.size(), tileRes, &ibo, &generator, &renderer) );
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


};