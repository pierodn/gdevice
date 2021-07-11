#pragma once 

// Window
#define WINDOW_WIDTH		(1280+16)	// 1280x720 =>	youtube HD
#define WINDOW_HEIGHT		(720+34)
#define FULLSCREEN_WIDTH	1280
#define FULLSCREEN_HEIGHT	600
#define	FOV					90

// World scale
#define POV_HEIGHT			(0.180*4)
#define METERS_PER_TILE		2.0

// Terrain, heightmap, clipmap
#define TILE_RESOLUTION		17 //33 //2^N+1. NOTE: Change also in generate_terrain.glsl
#define CLIPMAPS_COUNT		5		// 12
#define CLIPMAP_WINDOW		10
#define CLIPMAP_BUFFERING	0
#define CLIPMAP_SIZE		(CLIPMAP_BUFFERING + 2 + CLIPMAP_WINDOW + 2 + CLIPMAP_BUFFERING)
#define CLIPMAP_KERNEL		(CLIPMAP_WINDOW/2)
#define CLIPMAP_ODDITY		(CLIPMAP_KERNEL%2)

// World parameters
#define FOG_DENSITY			0.00010  // 0.0004

// Rendering
#define TARGET_FPS          60
#define TESSELLATION		0.15	//0.15
#define NEAR_CLIP_PLANE		(POV_HEIGHT/100)
#define TEXTURE_RANGE		0.002

// Control parameters
#define ROTATION_SPEED		100.0
#define TRANSLATION_SPEED	(POV_HEIGHT / 1000.0)
#define BOOST_FACTOR        10.0
#define WARP_FACTOR         300.0
#define TIME_START          16.00
#define TIME_SPEED			1.00
#define FPS_FREQUENCY       5