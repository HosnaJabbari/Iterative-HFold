#ifndef W_FINAL_H_
#define W_FINAL_H_

#include "VM_final.h"
#include "V_final.h"
#include "pseudo_loop.h"
#include "s_min_folding.h"
#include "h_common.h"


class W_final: public s_min_folding{
	public:
		W_final(char *seq, char *res);
        // constructor for the restricted mfe case

        ~W_final ();
        // The destructor

        double hfold ();
        // PRE:  the init_data function has been called;
        //       the space for structure has been allocate
        // POST: fold sequence, return the MFE structure in structure, and return the MFE

		double hfold_pkonly ();
		// PRE:  the init_data function has been called;
		//       the space for structure has been allocate
		// POST: fold sequence, return the MFE structure in structure, and return the MFE

        void return_structure (char *structure) ;
        // writes the predicted MFE structure into structure


    protected:
    	// Hosna: June 18th, 2007:
        // this pointer is the main part of the Hierarchical fold program
        // and corresponds to WMB recurrence
        pseudo_loop *WMB;
        // pointer to the final V matrix
        V_final *v;

        // pointer to the final VM matrix
        VM_final *vm;

        void space_allocation();

        // allocate the necessary memory
        double fold_sequence_restricted ();

        void backtrack_restricted (seq_interval *cur_interval, str_features *fres);
        // backtrack, the restricted case

		void backtrack_restricted_pkonly (seq_interval *cur_interval, str_features *fres);
		// backtrack, the restricted case with pk only base pairs

        void compute_W_restricted (int j, str_features *fres);
        // fill the W array, the restricted case

		// Hosna, April 3, 2012
		void compute_W_restricted_pkonly(int j, str_features *fres);
		// fill the W array, with addition of just pseudoknotted base pairs to the original structure

        int compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch);

		// Hosna, April 3, 2012
		int compute_W_br2_restricted_pkonly (int j, str_features *fres, int &must_choose_this_branch);

        int compute_W_br3_restricted (int j, str_features *fres);

        void print_result();
        // PRE:  The matrix V has been calculated and the results written in f
        // POST: Prints details of each elementary structure

        //s_min_folding *s_w;
        //int seq_length;
        //int *intseq;
//        int *index;

};

#endif /*W_FINAL_H_*/
