#pragma once

#if defined(WIN32)
	#include "os/win32/window.h"
#elif defined(__APPLE__)
    #include "os/apple/applewindow.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "os/linux/linuxwindow.h"
#else
    #error unknown platform!
#endif
