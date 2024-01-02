#include <windows.h>

#define EXPORTING_GDEVICE_DLL
#include "gdevice.h"

#if 0
BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved )
{/*
    switch ( ul_reason_for_call )
    {
        case DLL_PROCESS_ATTACHED: // A process is loading the DLL.
        break;
        case DLL_THREAD_ATTACHED: // A process is creating a new thread.
        break;
        case DLL_THREAD_DETACH: // A thread exits normally.
        break;
        case DLL_PROCESS_DETACH: // A process unloads the DLL.
        break;
    }*/
    return TRUE;
}
#endif

//_declspec (dllexport) 
int summit(int a, int b)
{
    return(a+b);
}