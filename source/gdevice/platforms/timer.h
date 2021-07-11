#pragma once

#if defined _WIN32
	#include "platforms/win32/timer.h"
#elif defined __APPLE__
    #include "platforms/apple/appletimer.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "platforms/linux/linuxtimer.h"
#else
    #error unknown platform!
#endif
