#pragma once

#if defined _WIN32
	#include "platforms/win32/keyboard.h"
#elif defined __APPLE__
    #include "platforms/apple/applekeyboard.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "platforms/linux/linuxkeyboard.h"
#else
    #error unknown platform!
#endif