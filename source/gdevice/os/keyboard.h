#pragma once

#if defined _WIN32
	#include "os/win32/keyboard.h"
#elif defined __APPLE__
    #include "os/apple/applekeyboard.h"
#elif defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
    #include "os/linux/linuxkeyboard.h"
#else
    #error unknown platform!
#endif