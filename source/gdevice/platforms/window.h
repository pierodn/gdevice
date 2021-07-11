#pragma once

#if defined(WIN32)
	#include "platforms/win32/window.h"
#elif defined(__APPLE__)
    #include "platforms/apple/applewindow.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "platforms/linux/linuxwindow.h"
#else
    #error unknown platform!
#endif
