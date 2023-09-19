#ifndef HOTSPOT_HEADER
#define HOTSPOT_HEADER

#include <stdio.h>  
#include <assert.h> 
#include <stdlib.h>
#include <string>
#include <algorithm>

class Hotspot{
    public:
       
        //constructor
        Hotspot(int left_inner_index, int right_inner_index, int length){
            assert(left_inner_index >= 0);
            assert(right_inner_index >= left_inner_index);
            this->left_inner_index = left_inner_index;
            this->right_inner_index = right_inner_index;
            this->left_outer_index = left_inner_index;
            this->right_outer_index = right_inner_index;
            this->length = length;
            this->energy = 0.0; 
            this->size = 1;
            std::string temp(length,'_');
            structure = temp;

        }
        //Copy Constuctor
        Hotspot(const Hotspot &hotspot){
            this->left_inner_index = hotspot.get_left_inner_index();
            this->right_inner_index = hotspot.get_right_inner_index();
            this->left_outer_index = hotspot.get_left_outer_index();
            this->right_outer_index = hotspot.get_right_outer_index();
            this->length = hotspot.get_length();
            this->energy = hotspot.get_energy(); 
            this->size = hotspot.get_size();
            structure = hotspot.get_structure();
        }
        //destructor
        ~Hotspot(){
        }

        //getters
        int get_left_inner_index() const;
        int get_right_inner_index() const;
        int get_left_outer_index() const;
        int get_right_outer_index() const;
        int get_length() const;
        double get_energy() const;
        int get_size() const;
        std::string get_structure() const;

        //setters
        void set_energy(double energy);
        void set_structure();
        void set_structure(std::string structure);
        void set_default_structure();

        void move_left_outer_index();
        void move_right_outer_index();
        void increment_size();
        bool is_invalid_energy();
    
    private:
        int left_inner_index;
        int right_inner_index;
        int left_outer_index;
        int right_outer_index;
        double energy;
        int size;
        std::string structure;
        int length;
};

#endif
