#ifndef RESULT_HEADER
#define RESULT_HEADER

#include <stdio.h>  
#include <assert.h> 
#include <stdlib.h>
#include <string>

class Result{
    public:
        //constructor
        Result(std::string sequence,std::string restricted,double restricted_energy, std::string final_structure, double final_energy, int method_chosen);
        //destructor
        ~Result();

        //getter
        std::string get_sequence();
        std::string get_restricted();
        std::string get_final_structure();
        double get_restricted_energy();
        double get_final_energy();
        int get_method_chosen();

        struct Result_comp{
		bool operator ()(Result &x, Result &y) const {
			if(x.get_final_energy() < y.get_final_energy()) return true;

            if(x.get_restricted_energy() < y.get_restricted_energy()) return true;

            return false;
		}
		} result_comp;
        
    private:
        std::string sequence;
        std::string restricted;
        double restricted_energy;
        std::string final_structure;
        double final_energy;
        int method_chosen;
};



#endif
