#pragma once

#include <vector>
#include <string>

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

static const unsigned char BINDINGS[] = 
{ 
	F1, F2, F3, F4,
	F5, F6, F7, F8, 
	F9, F10, F11, F12, 
	'H', 
    'G', 'C', 'U', 'T', 'V' 
};
                       
struct Controls 
{
	static const unsigned char* Bindings; 

	enum constants {	WIREFRAME, DEBUGMODE, DIFFUSE, SPECULAR,
						FRESNEL, SKY, INDIRECT, SCATTERING, 
                        TESSELLATOR, BUMPS, SHADOWS, PBR, 
                        HEATMAP,
						GAMMA, CONTRAST, UNSATURATE, TINT, VIGNETTING, 
						CONTROLCOUNT };

	int* keyCounters; // Dependency on Key::getCounters()
    std::vector<std::string> literals[CONTROLCOUNT];

	void GetStatusString(char* string)
	{
		for( int i=0; i<CONTROLCOUNT; i++ )
		{
			int size = literals[i].size();
			int value = keyCounters[Bindings[i]] % (size <= 1 ? 2 : size);

			string[i]  =	size <= 0 ?	'?' :
							size == 1 ?	(value == 0 ? '-' : literals[i][0][0] ) :
										(value == 0 ? '-' : literals[i][value][0]);
		}
		string[CONTROLCOUNT] = 0;
	}

	void ShowLegenda()
	{
		color(CMD_WHITE, 0); 
        printh("CONTROLS");
		for( int i=0; i<CONTROLCOUNT; i++ )
		{
			bool fx  = VK_F1<= Bindings[i] && Bindings[i] <= VK_F9;
			bool fxx = VK_F1<= Bindings[i] && Bindings[i] <= VK_F24;
			if( fxx ) printf("F%i%s", Bindings[i] - VK_F1 + 1, fx ? " " : "");
				 else printf( "%c  ", Bindings[i] );

            printf(" => %s", literals[i][0].c_str() );

			if( literals[i].size() > 1 )
			{
				printf("\t[Off, ");
				for( unsigned int j=1; j<literals[i].size(); j++ ) 
				{
					printf("%s", (literals[i])[j].c_str());
					if( j<literals[i].size()-1 ) printf(", ");
				}
				printf("]");
			}
			printf("\n" );
		}
        DEBUG_PRINT("\n");
	}

    static Controls& GetInstance()
    {
        static Controls instance;
        return instance;
    }

private:
    Controls()
    {
        literals[WIREFRAME]	    .push_back("Wireframe");
        literals[DEBUGMODE]	    .push_back("DebugMode");
        literals[DEBUGMODE]	    .push_back("Color");
        literals[DEBUGMODE]	    .push_back("HeightBlend");
        literals[DEBUGMODE]	    .push_back("Normal");
        literals[DEBUGMODE]	    .push_back("Light");
		literals[DIFFUSE]		.push_back("Diffuse");		// Direct light: Lambertian (classic or some other)
        literals[SPECULAR]		.push_back("Specular");		// Direct light: Specular (classic or PBR)

        literals[FRESNEL]		.push_back("Fresnel");		// Ambient light: fresnel
		literals[SKY]			.push_back("Sky");			// Ambient light: sky
        literals[INDIRECT]		.push_back("Indirect");		// Ambient light: direct light bouncing back
		literals[SCATTERING]	.push_back("Scattering");	// Light scattering

		literals[TESSELLATOR]	.push_back("Tessellator");	// Micropolygons
		literals[BUMPS]	        .push_back("Bumps");
		literals[SHADOWS]		.push_back("Shadows");		// Low-Frequency lambertian filter
        literals[PBR]			.push_back("PBR");

        literals[HEATMAP]		.push_back("Heatmap");

		literals[GAMMA]			.push_back("Gamma");
		literals[CONTRAST]		.push_back("Contrast");
		literals[UNSATURATE]	.push_back("Unsaturate");
		literals[TINT]			.push_back("Tint");
		literals[VIGNETTING]	.push_back("Vignetting");

		//
		// default values
		//
        keyCounters = Key::getCounters(); // TODO fix

		keyCounters[Bindings[WIREFRAME]]	= 0;
        keyCounters[Bindings[DEBUGMODE]]	= 0;
		keyCounters[Bindings[DIFFUSE]]		= 1;
        keyCounters[Bindings[SPECULAR]]		= 1;
        
        keyCounters[Bindings[FRESNEL]]		= 1;
		keyCounters[Bindings[SKY]]			= 1;
        keyCounters[Bindings[INDIRECT]]		= 1;
		keyCounters[Bindings[SCATTERING]]	= 1;
        
		keyCounters[Bindings[TESSELLATOR]]	= 1;
		keyCounters[Bindings[BUMPS]]        = 1;
		keyCounters[Bindings[SHADOWS]]		= 1;
		keyCounters[Bindings[PBR]]			= 0;

        keyCounters[Bindings[HEATMAP]]		= 1;
		
		keyCounters[Bindings[GAMMA]]		= 1;
		keyCounters[Bindings[CONTRAST]]		= 0;
		keyCounters[Bindings[UNSATURATE]]	= 0;
		keyCounters[Bindings[TINT]]			= 0;
		keyCounters[Bindings[VIGNETTING]]	= 1;
	}
};

const unsigned char* Controls::Bindings = BINDINGS;