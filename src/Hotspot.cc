#include "hotspot.hh"

//getters
int Hotspot::get_left_inner_index() const{
    return this->left_inner_index;
}
int Hotspot::get_right_inner_index() const{
    return this->right_inner_index;
}
int Hotspot::get_left_outer_index() const{
    return this->left_outer_index;
}
int Hotspot::get_right_outer_index() const{
    return this->right_outer_index;
}
double Hotspot::get_energy() const{
    return this->energy;
}
int Hotspot::get_size() const{
    return this->size;
}
int Hotspot::get_length() const{
    return this->n;
}

std::string Hotspot::get_structure() const{
    return this->structure;
}

//setters
void Hotspot::set_energy(double energy){
    this->energy = energy;
}

void Hotspot::set_structure(){
    for(int i =1; i <= n; ++i){
        if(i >= left_outer_index && i <= left_inner_index){
            structure[i] = '(';
        }else if(i <= right_outer_index && i >= right_inner_index){
            structure[i] = ')';
        }else{
            structure[i] = '.';
        }
    }
    structure = structure.substr(1,n);
}


void Hotspot::set_structure(std::string structure){
    this->structure = structure;
}


void Hotspot::set_default_structure(){
    // memset(this->structure,'\0',this->length);
    for(int i =1; i <= this->n; i++){
        this->structure[i] = '.';
    }
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