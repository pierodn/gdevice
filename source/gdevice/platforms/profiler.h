#pragma once

#include "platforms/timer.h"

// TODO FIX
// See http://preshing.com/20111203/a-c-profiling-module-for-multithreaded-apis/
                               
class Profiler
{
    Timer timer;
    float expected;
    char* function;
    char* filename;
    int   line;
    bool  skip;
    bool  verbose;

public:
    Profiler(float expected, char* function, char* filename, int line, bool skip = false, bool verbose = false)
    {
        this->expected = expected;
        this->function = function;
        this->filename = filename;
        this->line = line;
        this->skip = skip;
        this->verbose = verbose;
/*
        if( skip ) return;

        if( verbose ) {
            _PRINT_PATH(function, filename, line); 
            color(0, CMD_BROWN); printf("PROFILING"); 
            color(CMD_DARKGRAY,0); printf(" expected=%f\n", expected); 
            color(CMD_WHITE,0);  
        }
*/
    }

    ~Profiler() 
    {
        if( skip ) return;
        
        float elapsed = timer.elapsed();
        bool exceeded = elapsed > expected;
        if( verbose || exceeded ) {
            _PRINT_PATH(function, filename, line); 
            color(CMD_LIGHTGRAY,0); 
            printf("elapsed ");
            color(exceeded ? CMD_LIGHTRED : CMD_GREEN, 0); printf("%f\n", elapsed); 
            color(CMD_WHITE,0);  
        }
    }    
};
