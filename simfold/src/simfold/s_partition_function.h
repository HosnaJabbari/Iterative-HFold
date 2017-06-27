/***************************************************************************
                          s_partition_function.h  -  description
                             -------------------
    begin                : Wed Mar 15 2006
    copyright            : (C) 2006 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PARTITION_FUNCTION_H
#define PARTITION_FUNCTION_H

//#include "s_stacked_pair.h"
#include "s_hairpin_loop.h"
#include "s_internal_loop.h"
#include "s_multi_loop.h"
#include "s_multi_loop_sub.h"



class s_partition_function
{
    public:

//        friend class s_stacked_pair;
//        friend class s_internal_loop;
//        friend class s_multi_loop;

        s_partition_function (char *seq, int ignore_dangles=1, int compute_g_dangles=0, char *restricted=NULL);
        // The constructor                

        ~s_partition_function ();
        // The destructor

        PFTYPE compute_partition_function ();
  
        void compute_base_pair_probabilities ();
        // Nov 9, 2006. Computes base pair probabilities
        // PRE: the partition function arrays have been filled

        void print_base_pair_probabilities (PFTYPE threshold = 0.0);

        PFTYPE get_probability (int i, int j) { return p[index[i]+j-i]; }
        
        void PS_dot_plot(char *wastlfile);        
        
        PFTYPE getOneoverRT () { return oneoverRT; }
        
        // exhaustive, for verification
        PFTYPE compute_partition_function_exhaustively ();
        PFTYPE Z;        
        PFTYPE Zexhaustive;
        int verify_partition_function();
        // PRE: two functions to compute the partition function were called:
        //      compute_partition_function();
        //      compute_partition_function_exhaustively();
        // POST: returns 1 if there's error, 0 otherwise      
        
        void verify_recurrences();
        
        void compute_logZ_gradient ();
        // this actually computes -RT*logZ_gradient, which is the expected counts
        
        void compute_logZ_gradient_exhaustively ();
        // this actually computes -RT*logZ_gradient, which is the expected counts
        
        int correct_gradient ();
        void copy_gradient (PFTYPE *grad);
        void copy_gradient_exhaustively (PFTYPE *grad);
        int correct_gradient_nan ();
        void print_gradient();        
        // return 1 if no derivative is nan
        void print_u();
        // print the u array

        void verify_simple_vs_complicated_bpp (s_partition_function *pf);
        // assume this one has ignore_dangles = 1 and otherpf has ignore_dangles=0 
        // assume that the input parameters for the object with ignore_dangles=0 has dangling ends 0
        // then the p base pair probabilities should be the same
        // this function verifies that

        void verify_simple_vs_complicated_recurrences (s_partition_function *otherpf);
        // assume this one has ignore_dangles = 1 and otherpf has ignore_dangles=0
        // assume that the input parameters for the object with ignore_dangles=0 has dangling ends 0
        // then the equivalent data structures should have the same values
        // this function verifies that

        void compute_logZ_gradient_finite_differences ();

        // arrays to compute the partition function
        PFTYPE *up;    // (i,j) paired
        PFTYPE *upm;    // (i,j) closes a multi-loop
        
        // arrays for exterior loops
        PFTYPE *u_ip_jp;
        PFTYPE *u_ip_ju;
        PFTYPE *u_iu_jp;
        PFTYPE *u_iu_ju;
        PFTYPE *s1_jp;
        PFTYPE *s1_ju;

        // arrays for multi-loop                
        PFTYPE *u1_ip_jp;
        PFTYPE *u1_ip_ju_jm1p;
        PFTYPE *u1_ip_ju;
        PFTYPE *u1_iu_jp;
        PFTYPE *u1_iu_ju_jm1p;
        PFTYPE *u1_iu_ju;
        PFTYPE *s2_jp;
        PFTYPE *s2_ju;
        PFTYPE *s3_jp;
        PFTYPE *s3_ju_jm1p;
        PFTYPE *s3_ju;

              
        PFTYPE *GlogZ;    // gradient of logZ


        int num_internal_in_up;
        int num_internal_in_p;
              
    private:

        char *csequence;
        char *restricted;           // a restricted string, with ( for paired, . for unpaired, _ for unrestricted
        int *ptable_restricted;         // corresponding to the restricted string
        
        int *sequence;             // the entire sequence for which we compute the energy
        int seqlen;                 // sequence length
        int *index;                // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
        int ignore_dangles;        // if 1, do not consider dangling ends at all - gradient will be 0.
        int compute_gradient_dangles;   // if 1, compute the gradient for the dangling ends, otherwise don't
        PFTYPE oneoverRT;
        PFTYPE expOneoverRT;
        
                
        // arrays for base pair probabilities        
        PFTYPE *p;    // base pair probabilities

        // the following are needed if the dangling ends are not considered at all (equivalent with the dangling ends being 0)
        PFTYPE *u;
        PFTYPE *s1;
        PFTYPE *u1;
        PFTYPE *s2;
        PFTYPE *s3;
        PFTYPE *pm;    // used if ignore_dangles
        PFTYPE *pm1;    // used if ignore_dangles
        PFTYPE *pm2;    // used if ignore_dangles

        // the next 6 are used if the dangling ends are considered.
        PFTYPE *pmnod3_needmidd3;
        PFTYPE *pmnod3_noneedmidd3;
        PFTYPE *pmd3_needmidd3;
        PFTYPE *pmd3_noneedmidd3;
        PFTYPE *pm1nod3_needendd3;
        PFTYPE *pm1d3_needendd3;
        
        // the next 4 are necessary to compute the derivative wrt multi_free_base_penalty
        PFTYPE *pm2d5_needmidd5;
        PFTYPE *pm2d5_noneedmidd5;
        PFTYPE *pm2nod5_needmidd5;
        PFTYPE *pm2nod5_noneedmidd5;
        
        void initialize_arrays ();
        
        void compute_u_ip_jp (int i, int j);
        void compute_u_ip_ju (int i, int j);
        void compute_u_iu_jp (int i, int j);
        void compute_u_iu_ju (int i, int j);
        void compute_s1_jp (int i, int j);
        void compute_s1_ju (int i, int j);
        
        void compute_up (int i, int j);
        void compute_upm (int i, int j);
        void compute_upm_nodangles (int i, int j);
                
        void compute_u1_ip_jp (int i, int j);
        void compute_u1_ip_ju_jm1p (int i, int j);
        void compute_u1_ip_ju (int i, int j);
        void compute_u1_iu_jp (int i, int j);
        void compute_u1_iu_ju_jm1p (int i, int j);
        void compute_u1_iu_ju (int i, int j);
        
        void compute_s2_jp (int i, int j);
        void compute_s2_ju (int i, int j);
        
        void compute_s3_jp (int i, int j);
        void compute_s3_ju_jm1p (int i, int j);
        void compute_s3_ju (int i, int j);
        
        PFTYPE exp_AUpenalty (int i, int j);
        PFTYPE exp_dangle5 (int i, int j, int k);
        PFTYPE exp_dangle3 (int i, int j, int k);
        PFTYPE fd3 (int jplus1, int h, int l);
        void compute_p (int h, int l);

        void compute_pmnod3_needmidd3 (int h, int l);   
        void compute_pmnod3_noneedmidd3 (int h, int l);   
        void compute_pmd3_needmidd3 (int h, int l);
        void compute_pmd3_noneedmidd3 (int h, int l);
        void compute_pm1nod3_needendd3 (int h, int l);
        void compute_pm1d3_needendd3 (int h, int l);
        
        void compute_pm2d5_needmidd5 (int h, int j);
        void compute_pm2d5_noneedmidd5 (int h, int j);
        void compute_pm2nod5_needmidd5 (int h, int j);
        void compute_pm2nod5_noneedmidd5 (int h, int j);                


        // the next functions are called only if ignore_dangles is true
        void compute_pm (int h, int l);    // used if ignore_dangles
        void compute_pm1 (int h, int l);    // used if ignore_dangles
        void compute_pm2 (int h, int l);    // used if ignore_dangles
        void compute_u (int i, int j);
        void compute_u1 (int i, int j);
        void compute_s1 (int i, int j);
        void compute_s2 (int i, int j);
        void compute_s3 (int i, int j);
        
        PFTYPE *pexhaustive;
        PFTYPE *uexhaustive;
        PFTYPE *upexhaustive;
        int has_base_pair (int i, int j, char *structure);
        int identical_structure (int i, int j, char *structure1, char *structure2);
        PFTYPE exp_free_energy_partial (int i, int j, char *sequence, char *structure, int removeAU);
                
        PFTYPE *GlogZexhaustive;

        PFTYPE *GlogZ_finite_differences;
        PFTYPE eAU; // exp AUpenalty, so that we only compute it once. - seems faster
        PFTYPE EXPA;
        PFTYPE EXPB1;
        PFTYPE EXPB2;
        PFTYPE *EXPC;
        PFTYPE edangle3[NUCL][NUCL][NUCL];  // fill these arrays from the beginning, so that we don't call a function - it's actually a looot faster (about 5 times on length 200)
        PFTYPE edangle5[NUCL][NUCL][NUCL];
        int num_substructures;  // the number of substructures, as computed by simfold suboptimal
        

//         // I'm using the terminology from Ding & Lawrence NAR 2003, "A statistical sampling algorithm for RNA secondary structure prediction"
//         double *u;     // partition function on region i,j (can have anything)
//         //double *unod;  // partition function on region i,j (can have anything), but which does not consider the dangling ends at all
//         double *unod5;  // partition function on region i,j (can have anything), but which does not consider the dangling ends at the most 5' end
//         double *unod3;  // partition function on region i,j (can have anything), but which does not consider the dangling ends at the most 3' end
//         double *unod35;  // partition function on region i,j (can have anything), but which does not consider the dangling ends at the most 3' end and at the most 5' end
// 
//         double *s1;    // helper
//         //double *s1nod;    // helper corresponding to unod
//         double *s1nod5;    // helper corresponding to unod5
//         double *s1nod3;    // helper corresponding to unod3
//         double *s1nod35;    // helper corresponding to unod35
// 
//         // arrays for multi-loops
// 
//         double *s2;    // helper
//         double *s2nod3;    // helper
//         double *u1;    // partial multi-loop
//         double *u1nod5;    // partial multi-loop, no dangling ends to the 5' end
//         double *u1nod3;    // partial multi-loop, no dangling ends to the 3' end
//         double *u1nod35;    // partial multi-loop, no dangling ends to the 5' and 3' ends
//
//         double *s3;    // helper
//         double *s3nod5;    // helper, no dangling end to the 5' end
//         double *s3nod3;    // helper, no dangling end to the 3' end
//         double *s3nod35;    // helper, no dangling end to the 5' and 3' ends
// 
//         void compute_u (int i, int j);
//         //void compute_unod (int i, int j);
//         void compute_unod3 (int i, int j);
//         void compute_unod5 (int i, int j);
//         void compute_unod35 (int i, int j);
// 
//         void compute_u1 (int i, int j);
//         void compute_u1nod5 (int i, int j);
//         void compute_u1nod3 (int i, int j);
//         void compute_u1nod35 (int i, int j);
// 
//         void compute_s1 (int h, int j);
//         //void compute_s1nod (int h, int j);
//         void compute_s1nod5 (int h, int j);
//         void compute_s1nod3 (int h, int j);
//         void compute_s1nod35 (int h, int j);
//         void compute_s2 (int h, int j);        
//         void compute_s2nod3 (int h, int j);        
//         void compute_s3 (int h, int j);
//         void compute_s3nod5 (int h, int j);
//         void compute_s3nod3 (int h, int j);
//         void compute_s3nod35 (int h, int j);
        
};



#endif


