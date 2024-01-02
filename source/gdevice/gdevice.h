#pragma once

#ifdef EXPORTING_GDEVICE_DLL
#define GDEVICE_API __declspec(dllexport)
#else
#define GDEVICE_API __declspec(dllimport)
#endif

#ifdef __cpluscplus
extern "C" {
#endif


GDEVICE_API int summit(int, int);



#ifdef __cpluscplus
}
#endif