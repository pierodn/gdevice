// gdevice tech demo

#ifdef _MSC_VER
#pragma comment( lib, "opengl32.lib" )
#pragma warning(disable: 4244) // conversion with possible loss of data
#pragma warning(disable: 4305) // truncation from double to float
#pragma warning(disable: 4800) // forcing value to bool
#endif

#define TEXTURE_SIZE 256 // TEMP

#include "frameworks/application.h"
//#include "engines/graphic/renderer/camera.h"
//#include "renderer/rocktexture.h"
//#include "renderer/voxel.h"
#include "libraries/math/algebra.h"

// TEMP
#include <stdio.h>
//#include "renderer/opengl/teapot.h"
#include "engines/graphic/__mesh.h"





#include "libraries/math/functors.h"

float cup_x[] = { 0.0, 0.5, 1.0 };
float cup_y[] = {-0.5, 1.0, 0.7 };

float teapot_x[] = { 0.0, 0.5, 1.0 };
float teapot_y[] = {-0.5, 1.0, 0.7 };

float teapot_sil[] = { 0.00,0.6, 0.20,1.0, 0.80,0.7, 0.85,0.1, 0.90,0.2, 1.0,0.0 };
float teapot_spout_rotz[] = { 0.0,0, 0.5,45, 0.7,0, 1.0,0 };

vec3 teapot_spout[] = { vec3(0,0,0), vec3(1,0,0), vec3(2,0,0) };






class HelloworldApplication : public Application
{
private:
	float rotationSpeed;
	float translationSpeed;

    int range;

	Mesh wall;
	Mesh tube;
	Mesh cone;
	Mesh ball;
	Mesh ring;
	Mesh cube;
	Mesh polyhedra[7][7];
	Mesh cup;
	Mesh teapot;
	Mesh teapot_spout;



public:

	void onOpen( Renderer& renderer )
	{
        rotationSpeed = 500;
        translationSpeed = 0.04;
		renderer.camera.position = vec3(0);
		renderer.camera.orientation = vec3(0);

        range = 1;


		for( int i=2; i<7; i++ )
		for( int j=2; j<7; j++ )
		{
			polyhedra[i][j].genBall(i,j).setup();
		}

		wall.genWall( 16 ).mix(polyhedra[4][6],0.03).setupNormals().setupIndices( );
		tube.genTube(5,5, 0.1 ).setup(); // tube.vertices.xloop = false; tube.setup(); // FIX capping
		cone.genCone(19).setup();
		ball.genBall(19,13).setup();
		ring.genRing(23,11).setup();
		cube.genCube( 0.02, 0.02 ).setup(); // cube.resample<lerp>(20,20); // DEBUG feature or bug ?



		cup.genTube(46,22).tubemul( unit, unit,
			add<zero,interpolate<float, cup_x, cup_y, coserp>> 
			).setup();
		// cup.resample<hermite3>(22,11); //.resize(27,27); // DEBUG 

		// TODO
		// cup.mul<cylindrical>(unit*unit, unit*unit, unit*cubic<teapot> )

		static Timer timer;
		printf("delta=%f\n", timer.delta() );


//for(int zz=0;zz<10;zz++)
		teapot.genTube(46,22).scale(vec3(1,0.8,1)).tubemul( unit, unit, 
			add<zero,cubicinterpolate<teapot_sil>>
			).setup();
//		printf("delta=%f\n", timer.delta() ); // 2.29

/*
for(int zz=0;zz<10;zz++)
		teapot.genTube(46,22).scale(vec3(1,0.8,1)).test_shapeY( 0.00,0.6, 0.20,1.0, 0.80,0.7, 0.85,0.1, 0.90,0.2, 1.0,0.0 ).setup();
		printf("delta=%f\n", timer.delta() ); // 2.18
*/
/* //NOTE this is the right candidate for it allows runtime generation
for(int zz=0;zz<10;zz++)
		teapot.genTube(46,22).scale(vec3(1,0.8,1)).test2_shapeY(teapot_sil).setup();
		printf("delta=%f\n", timer.delta() ); // 2.19
*/

		teapot_spout.genTube(12).scale( vec3(0.3,1,0.3) )
			.twist_y(zero, zero, cubicinterpolate<teapot_spout_rotz> )
			.rotate(0,0,-45)
			.setup();

		//teapot_spout.yAxialShape( 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.99,1.0,1.01,1.02);
		// renderer.add( TEAPOT_SPOUT, Mesh.genTube(12).scale(vec3(0.4,1,0.4)) );

	}

	void onDraw( Renderer& renderer )
	{
		//
		// Control Camera
		//
		vec3 rotation;
		rotation.pitch = window.mouseDeltaY();
		rotation.yaw   = window.mouseDeltaX();
		rotation.roll  = 0;

		renderer.camera.rotate( rotation * rotationSpeed );

		const float boostSpeed = 20;
		float tSpeed = translationSpeed * (window.key(VK_LSHIFT) ? boostSpeed : 1.0 );
		if( window.key('W') ) renderer.camera.move( +tSpeed );
        if( window.key('S') ) renderer.camera.move( -tSpeed );
        if( window.key('A') ) renderer.camera.strafe( +tSpeed );
        if( window.key('D') ) renderer.camera.strafe( -tSpeed );
		if( window.key(' ') ) renderer.camera.position += vec3(0.0,+tSpeed,0.0);

		// no mouse system
		const float rotspeed = 0.002f;
		float rs = rotspeed * rotationSpeed * (window.key(VK_LSHIFT) ? boostSpeed : 1.0 );
		if( window.key(VK_LEFT) )  renderer.camera.rotate( vec3(   0, -rs, 0 ) );
        if( window.key(VK_RIGHT) ) renderer.camera.rotate( vec3(   0, +rs, 0 ) );
        if( window.key(VK_UP) )    renderer.camera.rotate( vec3( -rs,   0, 0 ) );
        if( window.key(VK_DOWN) )  renderer.camera.rotate( vec3( +rs,   0, 0 ) );


        //
		// Rendering
		//
        glClearColor( 0.0f, 0.0f, 0.5f ,1.0f );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

/*		http://www.gamedev.net/reference/articles/article1591.asp

		// First I like to draw the subject. Basically you are just drawing the subject directly in front of you
		glTranslatef(0.0f, 0.0f, -cameradist);
		DrawObject();

		//You reset the view and now any scenery or world that you want to draw is drawn.
		glLoadIdentity();
		glRotatef(360.0f - yaw, 0.0f, 1.0f, 0.0f);
		glTranslatef(-camerax, 0.0f, -cameraz);
		DrawScenery();
*/


		glLoadIdentity();

		// orientation

        glRotatef( renderer.camera.orientation.pitch, 1.0f, 0.0f, 0.0f );
		glRotatef( renderer.camera.orientation.yaw,   0.0f, 1.0f, 0.0f );
		glRotatef( renderer.camera.orientation.roll,  0.0f, 0.0f, 1.0f );


		// scrolling
//
		vec3 location;
		location.x = (int)renderer.camera.position.x;
		location.y = (int)renderer.camera.position.y;
		location.z = (int)renderer.camera.position.z;
		vec3 scrolling = renderer.camera.position - location;

/*
		glTranslatef( 1.0-scrolling.x, 1.0-scrolling.y, 1.0-scrolling.z );


		vec3 r;

		for( r.x=-range; r.x<=range; r.x++ )
		for( r.y=-range; r.y<=range; r.y++ )
		for( r.z=-range; r.z<=range; r.z++ )
		{
			vec3 voxel = r + location;

			if( voxel.x==0 && voxel.y==0 && voxel.z==0 )
			{
				glColor3f( 1.0f, 0.0f, 0.0f );
				glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
			}
			else
			{
				glColor3f( 0.0f, 0.0f, 1.0f );
				glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			}

			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			glBegin( GL_QUAD_STRIP );
				glVertex3d(r.x+0.0, r.z+0.0, r.y+0.0);
				glVertex3d(r.x+0.0, r.z+1.0, r.y+0.0);

				glVertex3d(r.x+1.0, r.z+0.0, r.y+0.0);
				glVertex3d(r.x+1.0, r.z+1.0, r.y+0.0);

				glVertex3d(r.x+1.0, r.z+0.0, r.y+1.0);
				glVertex3d(r.x+1.0, r.z+1.0, r.y+1.0);

				glVertex3d(r.x+0.0, r.z+0.0, r.y+1.0);
				glVertex3d(r.x+0.0, r.z+1.0, r.y+1.0);

				glVertex3d(r.x+0.0, r.z+0.0, r.y+0.0);
				glVertex3d(r.x+0.0, r.z+1.0, r.y+0.0);
			glEnd();


			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
			
			renderer.useShader( 1 );
			//renderer.setProgramVariable( 1, "color0", 0.0, 0.8, 0.0, 1.0 );
			//renderer.setProgramVariable( 1, "color1", 0.0, 0.0, 0.0, 1.0);
			//renderer.setProgramVariable( 1, "color2", 0.0, 0.8, 0.0, 1.0);

            //glPushMatrix();
            //glRotatef(theta, 1.0f, 0.0f, 1.0f);
            //drawSphere( 4, 5, 1.0 );
			prism4->draw();
            //glPopMatrix();

		}

		//glPushMatrix();
*/

		glTranslatef( -renderer.camera.position.x, -renderer.camera.position.y, -renderer.camera.position.z );
/*
		glPushMatrix();
        glTranslatef(2, -10, -3 );
        glRotatef(0.2, 1, 0, 0);
        glRotatef(0.1, 0, 1, 0);
		drawTeapot();
		glPopMatrix();
*/
		renderer.useShader(1);

			for( int i=2; i<7; i++ )
			for( int j=2; j<7; j++ )
			{
				Node node( &polyhedra[i][j] );
				node.place( i*3.0,j*3.0, i*2.0 );
				node.orientate( 0.0,0.0,0.0 );
				renderer.draw( node );
			}

			//Node node
			//wall.position = vec3( -6, 0, -8 );
			renderer.draw( Node(&wall).place(-6,0,-8) );

			//tube.position = vec3( -4, 0, -8 );
			renderer.draw( Node(&tube).place(-4,0,-8) );

			//cone.position = vec3( -2, 0, -8 );
			renderer.draw( Node(&cone).place(-2,0,-8) );

			//ball.position = vec3( 0, 0, -8 );
			renderer.draw( Node(&ball).place(0,0,-8) );

			//ring.position = vec3( 2, 0, -8 );
			renderer.draw( Node(&ring).place(2,0,-8) );

			//cube.position = vec3( 4, 0, -8 );
			renderer.draw( Node(&cube).place(4,0,-8) );
		
			//cup.position = vec3( 6, 0, -8 );
			renderer.draw( Node(&cup).place(6,0,-8) );

			//teapot.position = vec3( 6, 2, -8 );
			renderer.draw( Node(&teapot).place(6,2,-8) );

			//teapot_spout.position = vec3( 8, 2, -8 );
			renderer.draw( Node(&teapot_spout).place(8,2,-8) );

		renderer.useShader(0);

		
		
		


		// 
		// Misc controls
		//
		for(int i=0; i<renderer.getFlagCount(); i++ )
			if( window.key(VK_F1+i)==1 ) renderer.changeFlagState(i);

		//
		// GUI panel
		//
		glDisable(GL_LIGHTING);

		glLoadIdentity();
		glTranslatef( 0.0f, 0.0f, -0.5f );
		glColor3f( 1.0f, 1.0f, 1.0f );

		renderer.moveCursorAt(0,0);

		renderer.print( "camera.orientation = %s", str(renderer.camera.orientation) );
		renderer.print( "   camera.position = %s", str(renderer.camera.position) );
		renderer.print( "     tile location = %s", str(location) );
		renderer.print( "    tile scrolling = %s", str(scrolling) );
		renderer.print( "");
		renderer.print( "--- RENDERING CONTROL FLAGS ---" );
		for(int i=0; i<renderer.getFlagCount(); i++ )
			renderer.print( "[F%i] %10s = %s", i+1, renderer.getFlagName(i), renderer.getFlagStateDesc(i) );
		renderer.print( "" );

		renderer.print( "--- STATS ---" );
		renderer.print( "vertices = %i", renderer.vertices );
        renderer.print( "     fps = %d", renderer.fps );
		
		glEnable(GL_LIGHTING);
	}

};


int main ()
{
	//test_glsl(); {fgetc(stdin);}
	return HelloworldApplication().run();
}



