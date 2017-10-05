#include "Result.h"

//constructor
Result::Result(char* sequence,char* restricted, char* final_structure, double final_energy, int method_chosen){
    this->sequence = (char*) malloc(sizeof(char)*(strlen(sequence)+1));
    strcpy(this->sequence,sequence);

    this->restricted = (char*) malloc(sizeof(char)*(strlen(sequence)+1));
    strcpy(this->restricted,restricted);

    this->final_structure = (char*) malloc(sizeof(char)*(strlen(sequence)+1));
    strcpy(this->final_structure,final_structure);

    this->final_energy = final_energy;
    this->method_chosen = method_chosen;
}

//destructor
Result::~Result(){
    free(this->sequence);
    free(this->restricted);
    free(this->final_structure);
}

char* Result::get_sequence(){
    return this->sequence;
}
char* Result::get_restricted(){
    return this->restricted;
}
char* Result::get_final_structure(){
    return this->final_structure;
}
double Result::get_final_energy(){
    return this->final_energy;
}
int Result::get_method_chosen(){
    return this->method_chosen;
}



//getters