#pragma once

// TODO fix dependencies
#include <windows.h>
#define F1			VK_F1
#define F2			VK_F2
#define	F3			VK_F3
#define F4			VK_F4
#define F5			VK_F5
#define F6			VK_F6
#define F7			VK_F7
#define F8			VK_F8
#define F9			VK_F9
#define	F10			VK_F10
#define F11			VK_F11
#define F12			VK_F12
///////////////////////////////////////////

// TODO fix hard coded key bindings
static const unsigned char BINDINGS[] = 
{ 
	F1, F2, F3, F4,
	F5, F6, F7, F8, F9, F10, F11, F12, 
	'G', 'C', 'U', 'T', 'V' 
};
    
//
//   _____         _           _     
//  |     |___ ___| |_ ___ ___| |___ 
//  |   --| . |   |  _|  _| . | |_ -|
//  |_____|___|_|_|_| |_| |___|_|___|
//                                   
struct Controls 
{
	static const unsigned char* Bindings; 

	enum constants {	WIREFRAME, DEBUGMODE, BUMPS, TESSELLATOR, 
						DIFFUSE, SPECULAR, INDIRECT, SKY, FRESNEL, 
                        SHADOWS, VOLUMETRIC, SCATTERING,
						GAMMA, CONTRAST, UNSATURATE, TINT, VIGNETTING, 
						CONTROLCOUNT };

	int* values;
	Array<char*> literals[CONTROLCOUNT];

	void initialize() 
	{	
		// NOTE first value is the name of the uniform
		// NOTE "light0" stands for the main light source and is normally the sun
		literals[WIREFRAME]	    .push("Wireframe");
        literals[DEBUGMODE]	    .push("DebugMode").push("Color").push("HeightBlend").push("Normal").push("Light");
		literals[BUMPS]	        .push("Bumps");
		literals[TESSELLATOR]	.push("Tessellator");   // FIX
		literals[DIFFUSE]		.push("Diffuse");		// Lambertian
        literals[SPECULAR]		.push("Specular");	
		literals[INDIRECT]		.push("Indirect");		// Indirect lambertian contribute
		literals[SKY]			.push("Sky");			// Sky dome contribute	
		literals[FRESNEL]		.push("Fresnel");
        literals[SHADOWS]		.push("Shadows");		// Low-Frequency lambertian filter
        literals[VOLUMETRIC]	.push("Volumetric");
		literals[SCATTERING]	.push("Scattering");	// Fog + halo
		literals[GAMMA]			.push("Gamma");
		literals[CONTRAST]		.push("Contrast");
		literals[UNSATURATE]	.push("Unsaturate");
		literals[TINT]			.push("Tint");
		literals[VIGNETTING]	.push("Vignetting");

		// default values
		values[Bindings[WIREFRAME]]	    = 0;
        values[Bindings[DEBUGMODE]]	    = 0;
		values[Bindings[BUMPS]]         = 1;
		values[Bindings[TESSELLATOR]]	= 1;
		values[Bindings[DIFFUSE]]		= 1;
        values[Bindings[SPECULAR]]		= 1;
		values[Bindings[INDIRECT]]		= 1;
		values[Bindings[SKY]]			= 1;
		values[Bindings[FRESNEL]]		= 1;
        values[Bindings[SHADOWS]]		= 1;
		values[Bindings[VOLUMETRIC]]	= 0;
        values[Bindings[SCATTERING]]	= 1;
		values[Bindings[GAMMA]]			= 1;
		values[Bindings[CONTRAST]]		= 0;
		values[Bindings[UNSATURATE]]	= 0;
		values[Bindings[TINT]]			= 0;
		values[Bindings[VIGNETTING]]	= 1;
	}

	void getString(char* string)
	{
		for( int i=0; i<CONTROLCOUNT; i++ )
		{
			int size = literals[i].size();
			int value = values[Bindings[i]] % (size <= 1 ? 2 : size);

			string[i]  =	size <= 0 ?	'?' :
							size == 1 ?	(value == 0 ? '-' : literals[i][0][0] ) :
										(value == 0 ? '-' : literals[i][value][0]);
		}
		string[CONTROLCOUNT] = 0;
	}

	void showLegenda()
	{
		color(CMD_WHITE, 0); 
        printh("CONTROLS");
		for( int i=0; i<CONTROLCOUNT; i++ )
		{
			bool fx  = VK_F1<= Bindings[i] && Bindings[i] <= VK_F9;
			bool fxx = VK_F1<= Bindings[i] && Bindings[i] <= VK_F24;
			if( fxx ) printf("F%i%s", Bindings[i] - VK_F1 + 1, fx ? " " : "");
				 else printf( "%c  ", Bindings[i] );

			printf(" => %s", literals[i][0] );

			if( literals[i].size() > 1 )
			{
				printf("\t[Off, ");
				for( int j=1; j<literals[i].size(); j++ ) 
				{
					printf("%s", (literals[i])[j]);
					if( j<literals[i].size()-1 ) printf(", ");
				}
				printf("]");
			}
			printf("\n" );
		}
        DEBUG_PRINT("\n");
	}
};

const unsigned char* Controls::Bindings = BINDINGS;