#pragma once

#if defined _WIN32
	#include "os/win32/mutex.h"
#elif defined __APPLE__
    #include "os/apple/applemutex.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "os/linux/linuxmutex.h"
#else
    #error unknown platform!
#endif