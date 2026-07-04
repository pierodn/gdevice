#pragma once

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
    //#error unknown platform!
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>

//
// CRT Entry point
// 
int main();