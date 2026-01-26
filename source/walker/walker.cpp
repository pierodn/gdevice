#pragma once

#define DEBUG_SHOW_GLSL_SOURCE	

#include "application/assets/worlds/planet1/parameters.h" // CLIPMAP_WINDOW, TEXTURE_RANGE


#include "os/platform.h"
#include "os/application.h"
#include "os/window.h"
#include "os/keyboard.h"

#include "type/glsl.h"
#include "type/scene/oop/scene.h"
#include "type/scene/oop/terrain/heightmap.h"






class Walker : gd::Application<Walker>
{
    dmat3 camera;
    double speedFactor;

    Scene scene;
    Heightmap heightmap;
    float time; 
    vec3 sun;

public:

	void onOpen(Window<Walker>& window)
    {
        GL::Initialize();

        heightmap.setTileResolution(TILE_RESOLUTION);
	    heightmap.setLODs(CLIPMAPS_COUNT);
        heightmap.Initialize();

        scene.children.push(&heightmap);
        scene.transform.rotation = vec3(-90, 0, 0);
        scene.Initialize();

        speedFactor = 40;
	    camera.position = dvec3(0,0,10.0);
	    camera.rotation = dvec3(0,0,0);
    	
        time = TIME_START;
	    sun = rotate((time-7)*360/24, 0.7f, 1.7f) * vec3(0.1,-0.8,0.1); 

        // TODO add InputControl

        DEBUG_PRINT("\n");
		Controls::GetInstance().ShowLegenda();
    }

	void onSize(Window<Walker>& window)
    {
        scene.nodeState.ProjectionMatrix = projection( window.size, FOV, NEAR_CLIP_PLANE, heightmap.GetVisibilityDistance() );
    }

	void onDraw(Window<Walker>& window, double elapsed)
    {
        ///////////////////////////////
        // Update scene by user input.
        //

        // TODO camera = transform(camera, window, ROTATION_SPEED);
	        // Camera rotation
	        dvec3 mouse_rotation     = dvec3(window.mouseDeltaY(), 0, window.mouseDeltaX());
            dvec3 cursors_rotation   = dvec3(Key::cursorDeltaY(), 0, Key::cursorDeltaX()); // TODO window.HasFocus()
	        dvec3 input_rotation     = double(!window.isPointerVisible) * mouse_rotation + 0.01 * cursors_rotation;
	        camera.rotation   -= ROTATION_SPEED * input_rotation;		// counter-clockwise
	        camera.rotation.z  = mod(camera.rotation.z, 360.0);			// positive horizontal axis
	        camera.rotation.x  = clamp(camera.rotation.x, -90.0, +90.0);// clamped to [-90,+90] in order to avoid gimbal lock

            // Camera position
	        dmat3 T = transpose(rotate(camera.rotation));
	        dvec3 forward =  T.yAxis;
	        dvec3 strafe  = -T.xAxis;
	        double factor = (Key(LSHIFT).isPressed() ? BOOST_FACTOR : 1) * (Key(LCONTROL).isPressed() ? WARP_FACTOR : 1) * pow(speedFactor, 2.0);
	        double speed = factor * TRANSLATION_SPEED;
	        camera.position += forward*(Key('W').isPressed()? +speed : Key('S').isPressed()? -speed : 0) + strafe*(Key('A').isPressed()? +speed : Key('D').isPressed()? -speed : 0);
	        camera.position.z += (Key(' ').isPressed()? +speed : Key('Z').isPressed()? -speed : 0); // TODO could be included into strafe



        /////////////////////////////
        // Update time by user input

        // TODO
        // time = 
	    float timeSpeed = (Key('E').isPressed() ? +1 : Key('R').isPressed() ? -1 : 0) * float(TIME_SPEED) * (speedFactor);
	    time = fmod( time + timeSpeed*float(elapsed), 24.0f );
	    sun = rotate((time-7)*360.0f/24.0f, 40.0f, 0.0f) * vec3(0.2, -0.8, 0.1);

	    // TODO: Resolve this dependency.
	    for(int i=0; i<Controls::CONTROLCOUNT; i++)
        {
		    Key(Controls::Bindings[i]).isPressed();
	    }



        ////////////////////////////
        // Update scene
        //
        double level = heightmap.moveAt(camera); // moves the world around the camera


	    ////////////////////////////
	    // Rendering
	    //
        GL::SetViewport(window.size); // TODO GL::SetTarget<rgba>(window)
        GL::ClearBuffer();

        scene.nodeState.ModelViewMatrix = mat4(1);
	    scene.nodeState.inverseRotationMatrix = transpose(RotationMatrix(scene.transform.rotation) * RotationMatrix(heightmap.transform.rotation));
        scene.sceneState.sun = sun;

        // ****  NOTE  **** ////////////
        // Do not send in pipeline a node the very first time is generated as it might not be updated (and valid) yet.
        // From 2nd generation on, the node might be updated or not, yet valid. In case it is not updated, it will be the previous version.
        // This is workaround to avoid waiting for glDispatchCompute() to finish (syncronization GPU/CPU).
        static bool skipRenderingOnce = true;
        int vertexCount = scene.Traverse( skipRenderingOnce );

        // TODO place this node in th scene graph
        scene.RenderSkydome();

        skipRenderingOnce = false;


	    // 
	    // Misc controls
	    //
	    speedFactor = clamp(speedFactor + window.mouseDeltaWheel()/200.0, 0.4, 31.6228);
	    if(Key(ESCAPE).isJustPressed()) window.togglePointer();
        if(Key('L').isJustPressed()) window.toogleFullscreen();
	    if(Key('X').isPressed()) window.close();

        //
        // Show rendering state
        //
	    char controls[255];
	    Controls::GetInstance().GetStatusString(controls);

	    static Timer fpsTimer;
	    static dvec3 previous_pos; 
	    static int counter = 0;
	    static float interval;
	    if((++counter) % FPS_FREQUENCY == 0) 
            interval = fpsTimer.elapsed()/FPS_FREQUENCY;
	    int fps = (int)(1.0/interval);

		char camera_position[256];
		STR(METERS_PER_TILE * camera.position, camera_position);

		char camera_rotation[256];
		STR(camera.rotation, camera_rotation);

        // TODO FPS=%i%s CPU=%i%s GPU=%i%s
	    window.setTitle("FPS=%i%s Debug=[%s] Speed=%.2fkmh (x%i) Time=%02i.%02i Location=(%s) Direction=(%s)", 
		    fps, fps<100 ? "  " : "",
            controls, 
		    distance(camera.position, previous_pos) * METERS_PER_TILE*100/interval * 3600/1000,
		    int(speedFactor),
		    int(time), int(fract(time)*60),
		    camera_position,
			camera_rotation
	    );
	    previous_pos = camera.position;

        // TODO DEBUG_TRACE_ONCE(renderer.verticesCount);
    }

    int run() {
        return window.runPlainMessageLoop();
    }

private:
// TODO make private static functions for
// - Updating the scene by using the user input
// - Rendering the scene
// - Display debug stuff
// TODO make these functions tasks that can run in CPU or GPU as system functions over components.


};



int main()
{
    return Walker().run();
}
