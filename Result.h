#ifndef RESULT_HEADER
#define RESULT_HEADER

#include <stdio.h>  
#include <assert.h> 
#include <stdlib.h>
#include <string.h>

class Result{
    public:
        //constructor
        Result(char* sequence,char* restricted, char* final_structure, double final_energy, int method_chosen);
        //destructor
        ~Result();

        //getter
        char* get_sequence();
        char* get_restricted();
        char* get_final_structure();
        double get_final_energy();
        int get_method_chosen();
    private:
        char* sequence;
        char* restricted;
        char* final_structure;
        double final_energy;
        int method_chosen;
};

#endif
