#pragma once

// TODO ////////////

// TODO output to file
#include "platforms/platform.h"
//#include <stdio.h> // printf
//#include <string.h> // strlen
//#include "libraries/math/algebra.h" // min(a,b)

#define _TAB32 "%-32s"

struct Log 
{
/*
    // http://forums.codeguru.com/showthread.php?539241-RESOLVED-can-anyone-tell-me-all-colors-const-for-SetConsoleTextAttribute()-function
    static const int BLACK          = 0;
    static const int BLUE			= 1;
    static const int GREEN			= 2;
    static const int CYAN			= 3;
    static const int RED			= 4;
    static const int MAGENTA		= 5;
    static const int BROWN			= 6;
    static const int LIGHTGRAY		= 7;
    static const int DARKGRAY		= 8;
    static const int LIGHTBLUE		= 9;
    static const int LIGHTGREEN	    = 10;
    static const int LIGHTCYAN		= 11;
    static const int LIGHTRED		= 12;
    static const int LIGHTMAGENTA	= 13;
    static const int YELLOW		    = 14;
    static const int WHITE			= 15;

    // Win32 only
    static void setColor(int foreground, int background)
    {
         WORD wColor = ((background & 0x0F) << 4) + (foreground & 0x0F);
         SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), wColor);
    }
*/
/*
    inline static void error(char* message) {
        setColor(0,4);
        printf( "%s\n", message );
        setColor(15,0);
    }


    inline void debug(char* message) {
#if defined(_DEBUG)
        printf( "%s\n", message );
#endif
    }

    inline void debug(char* variable, char* value) {
#if defined(_DEBUG)
        printf( _TAB32 ": %s\n",   variable,  value );
#endif
    }
    */
/*
	void println( const char* string="" )
	{
		printf( "%s\n", string );
	}

	void printh( const char* header="" )
	{
		static char line[] = "------------------------------------------------------";
		println();
		// TODO avoid printing spaces around the header string when it is empty ("")
		printf( "------  %s  %s\n", header, line + min(6+2+2+strlen(header), sizeof(line)) );
		println();
	}
*/

};