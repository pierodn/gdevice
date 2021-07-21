#pragma once

#include <windows.h>
#include "type/glsl.h"

class Timer // : public Timer
{
private:
    double itime;
    int64 timeCounter;
    int64 elapsedCounter;
    int64 frequency;

public:

    Timer()
    {
        QueryPerformanceFrequency( (LARGE_INTEGER*)&frequency );
        reset();
    }

    void reset()
    {
        QueryPerformanceCounter( (LARGE_INTEGER*)&timeCounter );
        elapsedCounter = timeCounter;
        itime = 0.0;
    }

	double time()
    {
        int64 counter;
        QueryPerformanceCounter( (LARGE_INTEGER*)&counter );
        itime += ( counter - timeCounter ) / (double) frequency;
        timeCounter = counter;
        return itime;
    }

    static double absoluteTime()
    {
        int64 counter;
		int64 frequency;
        QueryPerformanceCounter( (LARGE_INTEGER*)&counter );
		QueryPerformanceFrequency( (LARGE_INTEGER*)&frequency );
        return (double) counter / (double) frequency;
    }

    double elapsed()
    {
        int64 counter;
        QueryPerformanceCounter( (LARGE_INTEGER*)&counter );
        double elapsed = (counter - elapsedCounter) / (double)frequency;
        elapsedCounter = counter;
        return elapsed;
    }

    void wait( double seconds )
    {
        Sleep( int(seconds*1000) );
    }
};
