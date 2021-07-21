#pragma once

#if defined _WIN32
	#include "os/win32/thread.h"
#elif defined __APPLE__
    #include "os/apple/applethread.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "os/linux/linuxtimer.h"
#else
    #error unknown platform!
#endif