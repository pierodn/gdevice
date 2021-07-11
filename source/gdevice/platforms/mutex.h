#pragma once

#if defined _WIN32
	#include "platforms/win32/mutex.h"
#elif defined __APPLE__
    #include "platforms/apple/applemutex.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "platforms/linux/linuxmutex.h"
#else
    #error unknown platform!
#endif