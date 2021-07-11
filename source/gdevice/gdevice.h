#pragma once

#define DEBUG_SHOW_GLSL_SOURCE	

#include "engines/graphic/heightmap.h" 
#include "platforms/platform.h"
#include "platforms/application.h"
#include "libraries/math/algebra.h"
#include "engines/graphic/node.h"

#include "parameters.h"
#include "platforms/keyboard.h"


class GDevice : public Application<GDevice>
{
	float time; 
	dmat3 camera; // TODO: mat3
	Node scene;
	Heightmap heightmap;
	Light sun;

    
	double wheelSpeed; // TODO: Get rid.
	
public:

	void onOpen(Window<GDevice>& window)
    {
        heightmap.setTileResolution(TILE_RESOLUTION);
	    heightmap.setLODs(CLIPMAPS_COUNT);	
	    heightmap.light = &sun;

        scene.transform.rotation = vec3( -90, 0, 0 );
	    scene.children.push( &heightmap );

        wheelSpeed = 4;
	    camera.position = 100.0 * dvec3(0,0,1);
	    camera.rotation = dvec3(0,0,0);
    	
        time = TIME_START;
        sun.id = 0;
	    sun.spot_direction = rotate((time-7)*360/24, 0.7f, 1.7f) * vec3(0.1,-0.8,0.1);   
    }

	void onSize(Window<GDevice>& window)
    {
	    window.renderer->setViewport( window.size ); 
	    window.renderer->setProjection( FOV, NEAR_CLIP_PLANE, heightmap.visibility() );
    }

	void onDraw(Window<GDevice>& window, double elapsed)
    {
        //DEBUG_PROFILE(0.020);
        Renderer& renderer = *window.renderer;
	    
        //////////////////////////
        // Update by user input.
        //
	    // Camera rotation
	    dvec3 mouse_rotation	= dvec3(window.mouseDeltaY(), 0, window.mouseDeltaX());
        dvec3 cursors_rotation	= dvec3(Key::cursorDeltaY(),  0, Key::cursorDeltaX());
	    dvec3 input_rotation	= double(!window.isPointerVisible) * mouse_rotation + 0.01 * cursors_rotation;
	    camera.rotation   -= ROTATION_SPEED * input_rotation;			// counter-clockwise
	    camera.rotation.z  = mod( camera.rotation.z, 360.0 );			// positive horizontal axis
	    camera.rotation.x  = clamp( camera.rotation.x, -90.0, +90.0 );	// clamped to [-90,+90] in order to avoid gimbal lock

        // Camera position
	    dmat3 T = transpose(rotate(camera.rotation));
	    dvec3 forward =  T.yAxis;
	    dvec3 strafe  = -T.xAxis;
	    float speedFactor = (Key(LSHIFT).isPressed() ? BOOST_FACTOR : 1) * (Key(LCONTROL).isPressed() ? WARP_FACTOR : 1) * (wheelSpeed*wheelSpeed);
	    double ts = TRANSLATION_SPEED * speedFactor;
	    camera.position += forward*(Key('W').isPressed()? +ts : Key('S').isPressed()? -ts : 0) + strafe*(Key('A').isPressed()? +ts : Key('D').isPressed()? -ts : 0);
	    camera.position.z += (Key(' ').isPressed()? +ts : Key('Z').isPressed()? -ts : 0); // TODO could be included into strafe

	    // Time
	    float speed = (Key('E').isPressed() ? +1 : Key('R').isPressed() ? -1 : 0) * float(TIME_SPEED) * (wheelSpeed);
	    time = fmod( time + speed*float(elapsed), 24.0f );
	    sun.spot_direction = rotate((time-7)*360/24, 40, 0)*vec3(0.2, -0.8, 0.1);
    	
	    // TODO: Resolve this dependency.
	    for(int i=0; i<Controls::CONTROLCOUNT; i++) {
		    Key(Controls::Bindings[i]).isPressed();
	    }

        ////////////////////////////
        // Update scene
        //
        double level = heightmap.moveAt( camera );
        heightmap.generateInvalidatedTiles(renderer);

        // TODO: Sleep while compute shader is working.
/*      double amountSleep = 1.0/TARGET_FPS - elapsed;
        DEBUG_TRACE(amountSleep);
        if( amountSleep>0 ) {
            Sleep( amountSleep ); 
        }*/

        

	    ////////////////////////////
	    // Rendering
	    //
        // TODO: renderer.barrier 
	    renderer.setFogDensity( FOG_DENSITY ); // TODO: Move to world properties.
	    renderer.clear();
	    renderer.inverseRotationMatrix = transpose( RotationMatrix(scene.transform.rotation) * RotationMatrix(heightmap.transform.rotation) );
        // TODO: renderer.setTarget<rgba>( NULL );
	    renderer.draw( scene );
	    renderer.drawSky();


	    // 
	    // Misc controls
	    //
	    wheelSpeed = clamp( wheelSpeed + window.mouseDeltaWheel()/200.0, 1.0, 31.6228 );

	    if( Key('L').isJustPressed() ) window.toogleFullscreen();
	    if( Key('P').isJustPressed() ) window.togglePointer();
	    if( Key(ESCAPE).isPressed() ) window.close();

	    char controls_string[255];
	    renderer.controls.getString(controls_string);

	    static Timer fpsTimer;
	    static dvec3 previous_pos; 
	    static int counter = 0;
	    static float interval;
	    if( (++counter) % FPS_FREQUENCY == 0 ) interval = fpsTimer.elapsed()/FPS_FREQUENCY;
	    int fps = (int)(1.0/interval);

	    window.setTitle( "FPS=%i%s Speed=%.2fkmh (x%i) Renderer=[%s] Time=%02i.%02i GPS=( %s)", 
		    fps, fps<100 ? "  " : "",
		    distance(camera.position, previous_pos) * METERS_PER_TILE/interval * 3600/1000,
		    int(speedFactor),
		    controls_string, 
		    (int)time, (int)(fract(time)*60),
		    str(METERS_PER_TILE * camera.position) 
	    );
	    previous_pos = camera.position;

        DEBUG_TRACE_ONCE( renderer.vertices );
    }

    int run() {
        return window.runPlainMessageLoop();
    }
};

