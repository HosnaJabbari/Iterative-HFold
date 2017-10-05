#include "Hotspot.h"

//constructor
Hotspot::Hotspot(int left_inner_index, int right_inner_index, int length){
    assert(left_inner_index >= 0);
    assert(right_inner_index >= left_inner_index);
    this->left_inner_index = left_inner_index;
    this->right_inner_index = right_inner_index;
    this->left_outer_index = left_inner_index;
    this->right_outer_index = right_inner_index;
    this->length = length;
    this->energy = 0.0; //todo kevin set it to hairpin of left_inner_index,right_inner_index?
    this->size = 1;
    this->structure = (char*) malloc(sizeof(char) * (length+1));

}

//destructor
Hotspot::~Hotspot(){
    free(this->structure);
}

//getters
int Hotspot::get_left_inner_index(){
    return this->left_inner_index;
}
int Hotspot::get_right_inner_index(){
    return this->right_inner_index;
}
int Hotspot::get_left_outer_index(){
    return this->left_outer_index;
}
int Hotspot::get_right_outer_index(){
    return this->right_outer_index;
}
double Hotspot::get_energy(){
    return this->energy;
}
int Hotspot::get_size(){
    return this->size;
}

char* Hotspot::get_structure(){
    return this->structure;
}

//setters
void Hotspot::set_energy(double energy){
    this->energy = energy;
}

void Hotspot::set_structure(){
    memset(this->structure,'\0',this->length+1);
    for(int i =0; i < this->length; i++){
        if(i >= this->left_outer_index && i <= this->left_inner_index){
            this->structure[i] = '(';
        }else if(i <= this->right_outer_index && i >= this->right_inner_index){
            this->structure[i] = ')';
        }else{
            this->structure[i] = '_';
        }
    }
    this->structure[this->length] = '\0';
}

void Hotspot::set_default_structure(){
    memset(this->structure,'\0',this->length);
    for(int i =0; i < this->length; i++){
        this->structure[i] = '_';
    }
    this->structure[this->length] = '\0';
}

void Hotspot::move_left_outer_index(){
    this->left_outer_index -= 1;
}

void Hotspot::move_right_outer_index(){
    this->right_outer_index += 1;
}

void Hotspot::increment_size(){
    this->size += 1;
}

bool Hotspot::is_invalid_energy(){
    return (this->energy >= 0.0);
}