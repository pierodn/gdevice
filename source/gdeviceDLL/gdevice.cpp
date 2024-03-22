
#define __declspec(dllexport)

extern "C"
{
    __declspec(dllexport) void* gdCreateWindow()
    {
        return 0;
    }
}