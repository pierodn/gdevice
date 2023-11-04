#if 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

#include "gdevice.h"

int main()
{
    return GDevice().run();
}

#else

#include "type/glsl_tests.h"

int main()
{
    test_all();
}

#endif