#ifndef HOTSPOT_HEADER
#define HOTSPOT_HEADER

#include <stdio.h>  
#include <assert.h> 
#include <stdlib.h>
#include <string.h>
#include <algorithm>

class Hotspot{
    public:
        //constructor
        Hotspot(int left_inner_index, int right_inner_index, int length);
        //destructor
        ~Hotspot();

        //getters
        int get_left_inner_index();
        int get_right_inner_index();
        int get_left_outer_index();
        int get_right_outer_index();
        double get_energy();
        int get_size();
        char* get_structure();

        //setters
        void set_energy(double energy);
        void set_structure();
        void set_default_structure();

        void move_left_outer_index();
        void move_right_outer_index();
        void increment_size();
        bool is_invalid_energy();
        /*
        is_complete_hotspot();
        increment_size();
        set_complete_hotspot();
        set_inner_begin_index(cur_interval->i);
        set_inner_end_index(cur_interval->j);
        set_energy(cur_interval->i,cur_interval->j, V);
        is_strong_stack();
        is_valid_energy(cur_interval->i,cur_interval->j);
        get_energy();
        */
    private:
        int left_inner_index;
        int right_inner_index;
        int left_outer_index;
        int right_outer_index;
        double energy;
        int size;
        char* structure;
        int length;
};

#endif
