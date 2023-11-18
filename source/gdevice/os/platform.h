#pragma once

// TODO
#if defined(WIN32)
	//#include "platform/win32/win32window.h"
    #define WIN32_LEAN_AND_MEAN
    #define WIN32_EXTRA_LEAN
    #include <windows.h>
#elif defined(__APPLE__)
    //#include "platform/apple/applewindow.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    //#include "platform/linux/linuxwindow.h"
#else
    #error unknown platform!
#endif

//////////////////////////////////
// Configuration
#if defined(_DEBUG)
	#pragma message(" -----> DEBUG")
#endif

#if defined(_CONSOLE)
	#pragma message(" -----> CONSOLE")
	#pragma comment(linker, "/SUBSYSTEM:CONSOLE")
#else
	#pragma comment(linker, "/SUBSYSTEM:WINDOWS")
#endif

#if defined(_MSC_VER)
	#pragma comment( lib, "opengl32.lib" )
	#pragma warning(disable: 4996) // unsafe(?) crt functions 
	#pragma warning(disable: 4244) // conversion with possible loss of data
	#pragma warning(disable: 4305) // truncation from double to float
	#pragma warning(disable: 4309) // truncation of constant value
	#pragma warning(disable: 4800) // forcing value to bool
#elif
	#error unknown compiler!
#endif

//
// CRT Entry point
// 
int main();

///////////////////
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
//#include <stdarg.h> 

//#include "type/glsl.h"

/////////////////////////////////////////
// Console
/////////////////////////////////////////

#define TAB32 "%-32s"

void printh( const char* header="" )
{
	static char line[] = "------------------------------------------------------";
	printf("\n");
	// TODO avoid printing spaces around the header string when it is empty ("")
	printf( "------  %s  %s\n", header, line + min(6+2+2+strlen(header), sizeof(line)) );
	printf("\n");
}


// http://forums.codeguru.com/showthread.php?539241-RESOLVED-can-anyone-tell-me-all-colors-const-for-SetConsoleTextAttribute()-function
static const int CMD_BLACK          = 0;
static const int CMD_BLUE			= 1;
static const int CMD_GREEN			= 2;
static const int CMD_CYAN			= 3;
static const int CMD_RED			= 4;
static const int CMD_MAGENTA		= 5;
static const int CMD_BROWN			= 6;
static const int CMD_LIGHTGRAY		= 7;
static const int CMD_DARKGRAY		= 8;
static const int CMD_LIGHTBLUE		= 9;
static const int CMD_LIGHTGREEN	    = 10;
static const int CMD_LIGHTCYAN		= 11;
static const int CMD_LIGHTRED		= 12;
static const int CMD_LIGHTMAGENTA	= 13;
static const int CMD_YELLOW		    = 14;
static const int CMD_WHITE			= 15;


static void color(int foreground, int background)
{
#if defined(WIN32)
     WORD wColor = ((background & 0x0F) << 4) + (foreground & 0x0F);
     SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), wColor);
#else
    // TODO Implement
#endif
}

inline char* _REMOVE_PATH(char* filename) {
	char* temp = max( strrchr(filename,'/'), strrchr(filename,'\\') );
	return temp ? temp+1 : filename;
}

inline void _PRINT_PATH(char* function, char* filename, int line) {
    color(CMD_GREEN,0); printf("%X ", GetCurrentThreadId()); 
    color(CMD_LIGHTGRAY,0); printf("%s:%d ", _REMOVE_PATH(filename), line);
    color(CMD_CYAN,0); printf("%s() ", function); 
    color(15,0);
}


inline void printf(int d)			{ printf("%d", d); };
inline void printf(unsigned int d)  { printf("%d", d); };
inline void printf(float f)			{ printf("%.2f", f); };
inline void printf(double f)		{ printf("%.4f", f); };
inline void printf(void* x)			{ printf("%08X", x); };

/////////////////////////////////////////
// Debug macros
/////////////////////////////////////////
    #define ASSERT(statement)     (void)((statement) || (_ASSERT(#statement, __FUNCTION__, __FILE__, __LINE__),0))
    #define CRITICAL(message)     _CRITICAL((message), __FUNCTION__, __FILE__, __LINE__)

// NOTE: Enabling DEBUG macros for RELEASE too, makes it not crash. TODO: Fix bug.
#if defined(_DEBUG)
    #define DEBUG_PRINT                 printf
    #define DEBUG_ASSERT(statement)     (void)((statement) || (_ASSERT(#statement, __FUNCTION__, __FILE__, __LINE__),0))
    #define DEBUG_CHECKPOINT_AUTO       _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__); color(CMD_WHITE,CMD_CYAN); printf(__COUNTER__); color(15,0);
    #define DEBUG_CHECKPOINT(tag)       _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__); color(CMD_WHITE,CMD_CYAN); printf("%s\n", #tag); color(15,0);
    #define DEBUG_PATH                  _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__);
	#define DEBUG_TRACE(variable)       _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__); color(CMD_LIGHTGRAY,0); printf("%s ", #variable); color(CMD_DARKGRAY,0); printf("= "); color(CMD_YELLOW,0); printf(variable); printf("\n");
    #define DEBUG_WARNING(message)      _WARNING((message), __FUNCTION__, __FILE__, __LINE__)
    #define DEBUG_CRITICAL(message)     _CRITICAL((message), __FUNCTION__, __FILE__, __LINE__)
    #define DEBUG_PROFILE(expected)     Profiler __profiler(expected, __FUNCTION__, __FILE__, __LINE__)
    #define DEBUG_CHECKPOINT_AUTO_ONCE  {static bool t = false; if(!t) { _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__); color(CMD_WHITE,CMD_CYAN); printf(__COUNTER__); color(15,0); t = true;}}
    #define DEBUG_CHECKPOINT_ONCE(tag)  {static bool t = false; if(!t) { _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__); color(CMD_WHITE,CMD_CYAN); printf("%s\n", #tag); color(15,0); t = true;}}
    #define DEBUG_TRACE_ONCE(variable)  {static bool t = false; if(!t) { _PRINT_PATH(__FUNCTION__, __FILE__, __LINE__); color(CMD_LIGHTGRAY,0); printf("%s ", #variable); color(CMD_DARKGRAY,0); printf("= "); color(CMD_YELLOW,0); printf(variable); printf("\n"); t = true;}}
    #define DEBUG_WARNING_ONCE(message) {static bool t = false; if(!t) { _WARNING((message),  __FUNCTION__, __FILE__, __LINE__); t = true;}}
    #define DEBUG_PROFILE_ONCE(expected) static bool t__FUNCTION__=false; Profiler __profiler(expected, __FUNCTION__, __FILE__, __LINE__, t__FUNCTION__, true); t__FUNCTION__=true;
    #define DEBUG_RUN_ONCE(code)        {static bool t = false; if(!t) { code; t = true;}}
#else // TODO write on a log file
    #define DEBUG_PRINT  
    #define DEBUG_ASSERT(statement)
    #define DEBUG_CHECKPOINT_AUTO
    #define DEBUG_CHECKPOINT(tag)
    #define DEBUG_PATH
    #define DEBUG_TRACE(variable)
    #define DEBUG_WARNING(message)
    #define DEBUG_CRITICAL(message)     _CRITICAL((message), __FUNCTION__, __FILE__, __LINE__)
    #define DEBUG_PROFILE(expected)
    #define DEBUG_CHECKPOINT_AUTO_ONCE
    #define DEBUG_CHECKPOINT_ONCE(tag)
    #define DEBUG_TRACE_ONCE(variable)
    #define DEBUG_WARNING_ONCE(message)
    #define DEBUG_PROFILE_ONCE(expected)
    #define DEBUG_RUN_ONCE(code)          {static bool t = false; if(!t) { code; t = true;}}
#endif

inline void _ASSERT(char* statement, char* function, char* filename, int line)
{
    _PRINT_PATH(function, filename, line);
    color(12,0); printf("ASSERTION FAILED: %s\n", statement );
    getc(stdin);
    exit(-1);
}

// When something happened but termination is not required.
inline void _WARNING(char* message, char* function, char* filename, int line)
{
    _PRINT_PATH(function, filename, line);
    color(12,0); printf("WARNING: %s\n", message ); color(15,0);
}

// When termination is required.
inline void _CRITICAL(char* message, char* function, char* filename, int line)
{
    _PRINT_PATH(function, filename, line);
    color(12,0); printf("CRITICAL ERROR: %s\n", message ); color(15,0);
	getc(stdin);
    exit(-1);
}

///////////////////////////////////////////////////////////
// String utilities
///////////////////////////////////////////////////////////

char* GetLastErrorAsString()
{
    char* buffer = NULL; // TODO Who will deallocate this buffer?
#if defined(WIN32)
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                                 NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&buffer, 0, NULL);
#else
    // TODO Implement
#endif
    return buffer;
}

char* strstr2( char* src, const char* sub )
{
	for( int i; *src; src++ )
	{
		for( i=0; src[i]==sub[i] && src[i]; i++ ); 
		if( !sub[i] ) return src;
	}
	return 0;
}

#if !defined(CONSOLE)
	int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
	{
		return main();
	}
#endif







