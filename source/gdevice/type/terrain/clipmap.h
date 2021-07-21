#pragma once

#include "parameters.h"
#include "type/node.h"

struct Clipmap : public Node
{
	Node tiles[ CLIPMAP_SIZE * CLIPMAP_SIZE ];

	int lod;
    ivec2 coarserTile;
    int quadrantInCoarserTile;
    

	dvec2 previousLocation;

	inline Node& getTile( int i, int j ) 
	{
		return tiles[i+j*CLIPMAP_SIZE];
	}

	Clipmap( Node* parent, int lod, int tileRes, IndexBuffer* ibo )
	{
		this->parent = parent;
		this->lod = lod;
		this->tileSize = 1<<lod;
		//this->material = NULL; // inherits from father (heightmap)

		//transform.position = vec3(0,0,0); // modified by move()
		transform.scale = vec3( tileSize, tileSize, 1 );
		
		// Forces invalidation
		previousLocation = dvec2(1E10,1E10);

		for(int i=0; i<CLIPMAP_SIZE; i++)
		for(int j=0; j<CLIPMAP_SIZE; j++) {
			Node& tile = getTile(i,j);
			tile.transform.position = vec3( i-CLIPMAP_SIZE/2, j-CLIPMAP_SIZE/2, 0 );
			tile.transform.rotation = vec3( 0, 0, 0 );
			tile.videomem_invalidated = true;
			//tile.hostmem_invalidated = true; //

			tile.vbo = new VertexBuffer();
			VertexBuffer& vbo = *tile.vbo;

			vbo.init( tileRes );
			vbo.mixmaps.scale.st = tileSize;

			tile.ibo = ibo;
			tile.lod_index = 0;

			tile.parent = this;
		}
	}

	~Clipmap()
	{
		for(int i=0; i<CLIPMAP_SIZE; i++)
		for(int j=0; j<CLIPMAP_SIZE; j++) {
			delete getTile(i,j).vbo;
		}
	}

	bool move( dvec2 location )
	{
		transform.position.x = scrollValue( location.x, 2*tileSize );
		transform.position.y = scrollValue( location.y, 2*tileSize );

		bool bx = transform.position.x >= tileSize;
		bool by = transform.position.y >= tileSize;
        quadrantInCoarserTile = int(bx) + int(by)*2;
 
		transform.position.xy /= tileSize;
		transform.position.xy -= float(CLIPMAP_ODDITY);

		int di = -(tileValue(location.x, 2*tileSize) - tileValue(previousLocation.x, 2*tileSize))*2;
		int dj = -(tileValue(location.y, 2*tileSize) - tileValue(previousLocation.y, 2*tileSize))*2;

		bool invalidation = di!=0 || dj!=0;
		
		if( invalidation ) {
			scrollTiles( di, dj, location );
		}

		updateChildren( bx, by );

		previousLocation = location;
		return invalidation;
	}

	// TODO: Toroidal access.
	void scrollTiles( int di, int dj, dvec2 location )
	{
		int i0,i1,is, j0,j1,js; 
		
		if( di>=0 ) {
			i0 = 0;
			i1 = CLIPMAP_SIZE;
			is = +1;
		} else {
			i0 = CLIPMAP_SIZE-1;
			i1 = -1;
			is = -1;
		}
		
		if( dj>=0 ) {
			j0 = 0;
			j1 = CLIPMAP_SIZE;
			js = +1;
		} else {
			j0 = CLIPMAP_SIZE-1;
			j1 = -1;
			js = -1;
		}
		
		for( int j=j0; j!=j1; j+=js )
		for( int i=i0; i!=i1; i+=is ) {
			int p = i + di;
			int q = j + dj;

			if( p>=0 && p<CLIPMAP_SIZE && q>=0 && q<CLIPMAP_SIZE) {
				Node& a = getTile(i,j);
				Node& b = getTile(p,q);
				VertexBuffer* t = a.vbo; a.vbo = b.vbo; b.vbo = t;
				bool  n = a.videomem_invalidated; a.videomem_invalidated = b.videomem_invalidated; b.videomem_invalidated = n;
				      n = a.hostmem_invalidated;  a.hostmem_invalidated = b.hostmem_invalidated;   b.hostmem_invalidated = n;
			} else {
				getTile(i,j).videomem_invalidated = true;
				getTile(i,j).hostmem_invalidated = true;
			}

			vec2 point;
			point.x = -(tileValue(location.x, 2*tileSize) + CLIPMAP_ODDITY/2.0f)* 2*tileSize + (i-CLIPMAP_SIZE/2)*tileSize;
			point.y = -(tileValue(location.y, 2*tileSize) + CLIPMAP_ODDITY/2.0f)* 2*tileSize + (j-CLIPMAP_SIZE/2)*tileSize;
			getTile(i,j).object_id = vec4( point, tileSize, 0.0 );
		}
	}

	void updateChildren( bool bx, bool by )
	{
		children.reset();

		// Part of the nested clipmap which is always present in lod==0
		int n1 = (CLIPMAP_SIZE - CLIPMAP_KERNEL +1)/2;			// 4 // 6
		int n2 = (CLIPMAP_SIZE + CLIPMAP_KERNEL -3)/2;			// 6 // 11

		if( lod==0 ) {
			for( int y=n1; y<=n2; y++ )
			for( int x=n1; x<=n2; x++ )
				children.push( &tiles[x + y*CLIPMAP_SIZE] );
		}

		if( lod==0 ||  bx ) {
			for( int y=n1; y<=n2; y++ )
				children.push( &tiles[n2+1 + y*CLIPMAP_SIZE] );
		}
		if( lod==0 ||  bx ||  by ) {
			children.push( &tiles[n2+1 + (n2+1)*CLIPMAP_SIZE] );
		}
		if( lod==0 ||         by ) {
			for( int x=n1; x<=n2; x++ )
				children.push( &tiles[x + (n2+1)*CLIPMAP_SIZE] );
		}

		if( lod==0 || !bx ||  by ) {
			children.push( &tiles[n1-1 + (n2+1)*CLIPMAP_SIZE] );
		}

		if( lod==0 || !bx ) {
			for( int y=n1; y<=n2; y++ )
				children.push( &tiles[n1-1 + y*CLIPMAP_SIZE] );
		}

		if( lod==0 || !bx || !by ) {
			children.push( &tiles[n1-1 + (n1-1)*CLIPMAP_SIZE] );
		}

		if( lod==0 ||        !by ) {
			for( int x=n1; x<=n2; x++ )
				children.push( &tiles[x + (n1-1)*CLIPMAP_SIZE] );
		}

		if( lod==0 ||  bx || !by ) {
			children.push( &tiles[n2+1 + (n1-1)*CLIPMAP_SIZE] );
		}

		for( int y=2+CLIPMAP_BUFFERING; y<=n1-2; y++ )
		for( int x=2+CLIPMAP_BUFFERING; x<=CLIPMAP_SIZE-2-CLIPMAP_BUFFERING-1; x++ )
			children.push( &tiles[x + y*CLIPMAP_SIZE] );

		for( int y=n1-1; y<=n2+1; y++ )
		for( int x=2+CLIPMAP_BUFFERING; x<=n1-2; x++ )
			children.push( &tiles[x + y*CLIPMAP_SIZE] );

		for( int y=n1-1; y<=n2+1; y++ )
		for( int x=n2+2; x<=CLIPMAP_SIZE-2-CLIPMAP_BUFFERING-1; x++ )
			children.push( &tiles[x + y*CLIPMAP_SIZE] );

		for( int y=n2+2; y<=CLIPMAP_SIZE-2-CLIPMAP_BUFFERING-1; y++ )
		for( int x=2+CLIPMAP_BUFFERING; x<=CLIPMAP_SIZE-2-CLIPMAP_BUFFERING-1; x++ )
			children.push( &tiles[x + y*CLIPMAP_SIZE] );
	}

};