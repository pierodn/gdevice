#pragma once

#include "engines/graphic/node.h"


struct Skybox : public Node
{
	Node faces[6];

	VertexBuffer vbo;
	IndexBuffer	 ibo;


	Skybox()
	{
		int tileRes = 2;

		vbo.init( tileRes );
		ibo.init( tileRes );
		
		for( int i=0; i<6; i++ )
		{
			faces[i].vbo = &vbo;
			faces[i].ibo = &ibo;
			faces[i].lod_index = 0;
		}
	}

	void init( int lods )
	{
		// TODO

	}

	~Skybox()
	{

	}

};