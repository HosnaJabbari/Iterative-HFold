#ifndef RESULT_HEADER
#define RESULT_HEADER

#include <stdio.h>  
#include <assert.h> 

class Result{
    public:
        Result(char* sequence,char* structure, char* final_structure, double final_energy, int method_chosen);
        ~Result();
};

#endif
