#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.hh"
#include "constants.hh"
#include <stdio.h>
#include <string.h>
#include "s_energy_matrix.hh"

class VM_final;
class V_final;
class pseudo_loop{

public:
	// constructor
	pseudo_loop(std::string seq, std::string restricted, s_energy_matrix *V, short *S, short *S1, vrna_param_t *params);

	// destructor
	~pseudo_loop();

    void compute_energies(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

    // energy_t get_energy(cand_pos_t i, cand_pos_t j);
	// in order to be able to check the border values consistantly
	// I am adding these get functions
	energy_t get_WI(cand_pos_t i, cand_pos_t j);
	energy_t get_WIP(cand_pos_t i, cand_pos_t j);

	energy_t get_VP(cand_pos_t i, cand_pos_t j);
	energy_t get_VPL(cand_pos_t i, cand_pos_t j);
	energy_t get_VPR(cand_pos_t i, cand_pos_t j);
	energy_t get_WMB(cand_pos_t i, cand_pos_t j);
	energy_t get_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree);

	energy_t get_WMBP(cand_pos_t i, cand_pos_t j);
	energy_t get_WMBW(cand_pos_t i, cand_pos_t j);

    void back_track(std::string structure, minimum_fold *f, seq_interval *cur_interval, sparse_tree &tree);

    void set_stack_interval(seq_interval *stack_interval);
    seq_interval *get_stack_interval(){return stack_interval;}
    std::string get_structure(){return structure;}
    minimum_fold *get_minimum_fold(){return f;}

private:

	cand_pos_t n;
	std::string res;
	std::string seq;

    s_energy_matrix *V;		        // the V object

	seq_interval *stack_interval;
	std::string structure;
	minimum_fold *f;
	vrna_param_t *params_;


	//Hosna
    std::vector<energy_t> WI;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot)
    std::vector<energy_t> VP;				// the loop corresponding to the pseudoknotted region of WMB
	std::vector<energy_t> VPL;				// the loop corresponding to the pseudoknotted region of WMB
    std::vector<energy_t> VPR;				// the loop corresponding to the pseudoknotted region of WMB
    std::vector<energy_t> WMB;				// the main loop for pseudoloops and bands
	std::vector<energy_t> WMBP; 				// the main loop to calculate WMB
	std::vector<energy_t> WMBW;
	std::vector<energy_t> WIP;				// the loop corresponding to WI'
    std::vector<energy_t> BE;				// the loop corresponding to BE
    std::vector<cand_pos_t> index;				// the array to keep the index of two dimensional arrays like WI and weakly_closed

	short *S_;
	short *S1_;

    // function to allocate space for the arrays
    void allocate_space();

    void compute_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: This function is supposed to fill in the WI array

	void compute_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the VP array

	// Computes the non-redundant recurrence from CParty (replaces VPP from original)
	void compute_VPL(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	void compute_VPR(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

	void compute_WMB(cand_pos_t i,cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the WMB array

	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)
	void compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// this is the helper recurrence to fill the WMB array

	// Computes the non-redundant recurrence from CParty (replaces WMBP case 2 from original)
	void compute_WMBW(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

	void compute_WIP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the WIP array

	void compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree);
	// Hosna: this function is supposed to fill the BE array

	// Hosna Feb 8th, 2007:
	// I have to calculate the e_stP in a separate function
	energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
	energy_t get_e_intP(cand_pos_t i,cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
	energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params);

  	// Hosna: Feb 19th 2007
  	// used for backtracking
  	void insert_node (cand_pos_t i, cand_pos_t j, char type);//, seq_interval *stack_interval);

};
#endif /*PSEUDO_LOOP_H_*/
