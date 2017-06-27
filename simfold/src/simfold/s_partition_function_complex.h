/***************************************************************************
                          s_partition_function_complex.h  -  description
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

#ifndef PARTITION_FUNCTION_COMPLEX_H
#define PARTITION_FUNCTION_COMPLEX_H

#include <complex>

// TODO: should this be long double?
typedef PARAMTYPE Real;
//typedef double Real;
using std::complex;
//typedef complex<Real> Complex;
typedef complex<double> Complex;

//#include "s_stacked_pair.h"
#include "s_hairpin_loop.h"
#include "s_internal_loop.h"
#include "s_multi_loop.h"
#include "s_multi_loop_sub.h"
#include "s_partition_function.h"

class s_partition_function_complex
{
    public:

//        friend class s_stacked_pair;
//        friend class s_internal_loop;
//        friend class s_multi_loop;

        //Complex function(Complex val);
        //s_partition_function_complex ();
        

        
        //#undef EXPA
        //#undef EXPB
        //#undef EXPC
        
        //#define EXPA   (exp (misc_complex.multi_offset * oneoverRT))
        //#define EXPB(X)   (exp ((((Complex)(X))*misc_complex.multi_helix_penalty) * oneoverRT))
        //#define EXPC(X)   (exp ((((Complex)(X))*misc_complex.multi_free_base_penalty) * oneoverRT))
        #define AU_penalty_complex(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?misc_complex.terminal_AU_penalty:0)

        
        s_partition_function_complex (char *seq, int ignore_dangles=1, int compute_g_dangles=0, int params_are_double=0);
        // The constructor                

        ~s_partition_function_complex ();
        // The destructor

        Complex compute_partition_function ();
        Complex Z;  // complex partition function



        // ignore next ones for now
  
        void compute_base_pair_probabilities ();
        // Nov 9, 2006. Computes base pair probabilities
        // PRE: the partition function arrays have been filled

        void print_base_pair_probabilities (PFTYPE threshold = 0);

        void print_partition_function();
        void print_u();
        
        
        Complex getOneoverRT () { return oneoverRT; }
        

        void compute_logZ_gradient ();
        void copy_gradient_numerical (PFTYPE *grad);
        void compute_logZ_gradient_finite_differences();
        void copy_gradient (PFTYPE *grad);
        int correct_gradient_nan ();
        void print_gradient();
        // return 1 if no derivative is nan

        void validate_partition_function (s_partition_function *pf);
        void compute_logZ_gradient_numerical ();
        void compute_logZ_gradient_numerical (int index);
        // only compute the partial derivative wrt index
        void compute_logZ_second_derivatives_numerical ();
        void compute_logZ_second_derivatives_numerical (int index1);
        // compute the second derivatives of logZ wrt theta_index1
        // this way it is easier to parallelize
        // unfortunately we can't really evaluate this
        void compute_logZ_second_derivatives_finite_differences ();
        
        int validate_gradient_numerical (s_partition_function *pf);
        int validate_gradient_dp (s_partition_function *pf);
        int validate_gradient_analytical_vs_numerical ();

        Complex simple_function (int which, Complex addition);
        void test_simple_function ();

        void print_Hessian ();
        // print the second derivative matrix to the screen
        void print_Hessian (char *filename);
        // print the second derivative matrix into the file filename
        void print_Hessian (int index1, char *filename);
        // print the second derivative matrix into the file filename
        // only print the row index1
        
              
    private:

        char *csequence;
        int *sequence;             // the entire sequence for which we compute the energy
        int seqlen;                 // sequence length
        int *index;                // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
        int ignore_dangles;        // if 1, do not consider dangling ends at all - gradient will be 0.
        int compute_gradient_dangles;   // if 1, compute the gradient for the dangling ends, otherwise don't                
        Complex oneoverRT;
        Complex RT;

        Complex EXPC (int i);
        
        int params_are_double;  // sometimes use more precision for better resolution

        // the following are needed if the dangling ends are not considered at all (equivalent with the dangling ends being 0)
        Complex *u;
        Complex *s1;
        Complex *u1;
        Complex *s2;
        Complex *s3;
        Complex *pm;    // used if ignore_dangles
        Complex *pm1;    // used if ignore_dangles
        Complex *pm2;    // used if ignore_dangles

        
        // arrays to compute the partition function
        Complex *up;    // (i,j) paired
        Complex *upm;    // (i,j) closes a multi-loop
        
        // arrays for exterior loops
        Complex *u_ip_jp;
        Complex *u_ip_ju;
        Complex *u_iu_jp;
        Complex *u_iu_ju;
        Complex *s1_jp;
        Complex *s1_ju;

        // arrays for multi-loop                
        Complex *u1_ip_jp;
        Complex *u1_ip_ju_jm1p;
        Complex *u1_ip_ju;
        Complex *u1_iu_jp;
        Complex *u1_iu_ju_jm1p;
        Complex *u1_iu_ju;
        Complex *s2_jp;
        Complex *s2_ju;
        Complex *s3_jp;
        Complex *s3_ju_jm1p;
        Complex *s3_ju;
                
        // arrays for base pair probabilities        
        Complex *p;    // base pair probabilities
        Complex *pmnod3_needmidd3;
        Complex *pmnod3_noneedmidd3;
        Complex *pmd3_needmidd3;
        Complex *pmd3_noneedmidd3;
        Complex *pm1nod3_needendd3;
        Complex *pm1d3_needendd3;
        
        // the next 4 are necessary to compute the derivative wrt multi_free_base_penalty
        Complex *pm2d5_needmidd5;
        Complex *pm2d5_noneedmidd5;
        Complex *pm2nod5_needmidd5;
        Complex *pm2nod5_noneedmidd5;
        
        void initialize_arrays ();
        
        void fill_data_structures_complex ();
        // Mirela: Mar 12, 2007
        // fill the complex parameter data structures from the global ones

        //void fill_data_structures_complex_double ();
        
        void compute_u_ip_jp (int i, int j);
        void compute_u_ip_ju (int i, int j);
        void compute_u_iu_jp (int i, int j);
        void compute_u_iu_ju (int i, int j);
        void compute_s1_jp (int i, int j);
        void compute_s1_ju (int i, int j);
        
        void compute_up (int i, int j);
        void compute_upm (int i, int j);
                
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

        // the next functions are called only if ignore_dangles is true
        void compute_pm (int h, int l);    // used if ignore_dangles
        void compute_pm1 (int h, int l);    // used if ignore_dangles
        void compute_pm2 (int h, int l);    // used if ignore_dangles
        void compute_u (int i, int j);
        void compute_u1 (int i, int j);
        void compute_s1 (int i, int j);
        void compute_s2 (int i, int j);
        void compute_s3 (int i, int j);
        void compute_upm_nodangles (int i, int j);
        
        
        Complex exp_AUpenalty (int i, int j);
        Complex exp_dangle5 (int i, int j, int k);
        Complex exp_dangle3 (int i, int j, int k);

        // the following are added in this complex version
        Complex asymmetry_penalty_complex (int size1, int size2);
        Complex penalty_by_size_complex (int size, char type);  
        Complex get_hairpin_energy (int i, int j, int* sequence, char *csequence);
        Complex get_stacked_energy (int i, int j, int *sequence);
        Complex get_internal_energy (int i, int j, int ip, int jp, int *sequence);
         
        Complex fd3 (int jplus1, int h, int l);
        void compute_p (int h, int l);
//         void compute_pm (int h, int l);    // used if ignore_dangles
//         void compute_pm1 (int h, int l);    // used if ignore_dangles
//         void compute_pm2 (int h, int l);    // used if ignore_dangles
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
        
        int has_base_pair (int i, int j, char *structure);
        int identical_structure (int i, int j, char *structure1, char *structure2);
        double exp_free_energy_partial (int i, int j, char *sequence, char *structure, int removeAU);
        
        Complex *GlogZ;    // gradient of logZ
        PFTYPE *GlogZ_numerical;    // gradient of logZ
        PFTYPE *GlogZ_finite_differences;    // gradient of logZ, as computed using the finite differences method
        PFTYPE **HlogZ_numerical;   // Hessian of logZ
        PFTYPE **HlogZ_finite_differences;   // Hessian of logZ


        // we need to use complex parameters

        //Complex eAU; // exp AUpenalty, so that we only compute it once. - seems faster
//         Complex EXPA;
//         Complex EXPB1;
//         Complex EXPB2;
        //Complex *EXPC;
        Complex edangle3[NUCL][NUCL][NUCL];  // fill these arrays from the beginning, so that we don't call a function - it's actually a looot faster (about 5 times on length 200)
        Complex edangle5[NUCL][NUCL][NUCL];

        Complex **parameters_complex;    // an array of pointers to all the parameters in the model, so that I don't take them case by case 
        void initialize_parameters_complex ();
        // make the parameters_complex array to point to the right parameters;

        Complex get_special_internal_complex (int i, int j, int ip, int jp);
        // Return the energy obtained when we consider 6 additional parameters for internal loop 3x3 and larger, 
        //  as described in Chen_Turner_2006b.
        // the arguments are positions in sequence      
                  
        Complex get_int22_complex_MODEL_EXTENDED (int ii, int jj, int kk, int ll, int mm, int nn, int oo, int pp);
        // get the energy for int22, when MODEL is EXTENDED

        Complex get_int11_complex_MODEL_EXTENDED (int ii, int jj, int kk, int ll, int mm, int nn);
        // get the energy for int11, when MODEL is EXTENDED

        Complex get_int21_complex_MODEL_EXTENDED (int ii, int jj, int kk, int ll, int mm, int nn, int oo);
        // return the energy for int21, when MODEL is EXTENDED

        
        typedef struct
        {
            // Extrapolation for large loops based on polymer theory
            // internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30)
            double param_greater30;    // we keep this fixed for learning   // let's not make this complex yet
            // helix
            Complex terminal_AU_penalty;
            // hairpin data from miscloop file
            Complex hairpin_GGG;          // bonus for GGG hairpin
            Complex hairpin_c1;           // c hairpin slope
            Complex hairpin_c2;           // c hairpin intercept
            Complex hairpin_c3;           // c hairpin of 3
            // internal loops
            Complex asymmetry_penalty_max_correction;
            Complex asymmetry_penalty_array[4];    //for param learning, only [0] and [1] are variable, [2] and [3] are fixed
            int gail_rule;  // this is not a complex parameter
            // multi-branched loops
            Complex multi_offset;
            Complex multi_helix_penalty;
            Complex multi_free_base_penalty;
            Complex intermolecular_initiation;
            // 3 parameters which replace the tstacki table - added on Dec 20, 2006
            Complex internal_AU_closure;   // does not need terminal_AU_penalty to be added to it
            Complex internal_GA_AG_mismatch;
            Complex internal_UU_mismatch;
            Complex internal22_delta_same_size;
            Complex internal22_delta_different_size;
            Complex internal22_delta_1stable_1unstable;
            Complex internal22_delta_AC;
            Complex internal22_match;   // if it has a CG or AU match in the middle
            Complex internal21_match;   // if it has a CG or AU or GU match in the middle
            Complex internal21_AU_closure;   // does not need terminal_AU_penalty to be added to it
            Complex internal11_basic_mismatch;
            Complex internal11_GG_mismatch;
            
//#if (MODEL == EXTENDED)
            // 4 params to replace the tstackh table if MODEL==EXTENDED and parsi_tstackh is 1
            Complex hairpin_AU_closure;   // does not need terminal_AU_penalty to be added to it
            Complex hairpin_AG_mismatch;
            Complex hairpin_GA_mismatch;
            Complex hairpin_UU_mismatch;
            
            // use GA and AG separately, and add GG, according to Schroeder_Turner_2000
            Complex internal_AG_mismatch;
            Complex internal_GA_mismatch;
            Complex internal_GG_mismatch;
    
            // 6 more parameters for internal loops 3x3 or larger, following Chen_Turner_2006b
            // I only added the ones that would not be included in tstacki
            // These parameters are to be used when !parsi_special
            Complex internal_special_3GA;     // 5'-YGGA/GAAR-3' or 5'-GGAR/YGAA-3', in loops 3x3 and larger
            Complex internal_special_2GA;     // 5'-GA/GA-3' next to a closing base pair, or 5'-GG/AA-3' next to a closing base pair, for 3x3, 3x4, 4x4, and 4x5 loops; 
                                //  ALSO 5'-RGGA/GAAY-3' or 5'-GGAY/YGAA-3' for 3x5, 3x6 and 4x6 loops.
                                // internal_2GA is USED only if internal_3GA was not used
            Complex internal_special_2xGA_GC;  // 5'-GANGC/GANGC-3' in 3x3 loops
            Complex internal_special_midGA;       // middle GA adjacent to RY in 3x3 loops
                                    // internal_midGA is USED only if none of internal_3GA and internal_2GA is used.
            Complex internal_special_UG_AG;       // once or twice, for each 5'-UG/AG-3' at the terminus of loops 3x3 or larger
            Complex internal_special_GU_A;        // first mismatch is GA, and U is 3' of G, for loops 3x3
            
            Complex internal11_AU_closure;
            Complex internal11_GU_closure;
            Complex internal11_AG_mismatch;
            Complex internal11_UU_mismatch;
            Complex internal11_5YRR_5YRR;
            Complex internal11_5RYY_5RYY;
            Complex internal11_5YYR_5YYR;
            Complex internal11_5YRY_5RYR;
            Complex internal11_5RRY_5RYY;    
            
            Complex internal21_initiation;
            Complex internal21_GU_closure;
            Complex internal21_AG_mismatch;    // applied once per loop, not applied to 5'RA/3'YG loops
            Complex internal21_GG_mismatch;    // applied once per loop
            Complex internal21_UU_mismatch;    // applied once per loop
            Complex internal22mid_group1;      // group 1 according to Christiansen_Znosko_2008
    // That is: a U · U pair adjacent to an R · R pair, a G · A  or A · G pair adjacent to a Y · Y pair, or  any combination of A · C, U · C, C · U,  C · C, C · A, or A · A pairs
            Complex internal22mid_group2;      // group 2 according to Christiansen_Znosko_2008
    // That is: any combination of adjacent G · A and A · G pairs or two U · U pairs
            Complex internal22mid_group3;      // group 3 according to Christiansen_Znosko_2008
    // That is: a U · U pair adjacent to a Y · Y (not U · U),  C · A, or A · C pair
            Complex internal22mid_group4;      // group 4 according to Christiansen_Znosko_2008
    // That is: a G · G pair not adjacent to a U · U pair
    
    // 2 more parameters for the case parsi_int22 == 1
            Complex internal22_AU_closure;     // as suggested by Christiansen_Znosko_2008
            Complex internal22_GU_closure;     // as suggested by Christiansen_Znosko_2008
            
//#endif
            
        } miscinfo_Complex;

        // info from tloop.dat
        typedef struct
        {
            //char *seq;  // pointer to the hairpin_tloop data structure, because no needed is necessary to the seq
            Complex energy;
        } hairpin_tloop_Complex;

        // energies information
        Complex stack_complex [NUCL] [NUCL] [NUCL] [NUCL];
        Complex tstackh_complex [NUCL] [NUCL] [NUCL] [NUCL];
        Complex tstacki_complex [NUCL] [NUCL] [NUCL] [NUCL];
        Complex int11_complex   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
        Complex int21_complex   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
        Complex int22_complex   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
        Complex dangle_top_complex  [NUCL] [NUCL] [NUCL];
        Complex dangle_bot_complex  [NUCL] [NUCL] [NUCL];
        Complex internal_penalty_by_size_complex [MAXLOOP+1];
        Complex bulge_penalty_by_size_complex [MAXLOOP+1];
        Complex hairpin_penalty_by_size_complex [MAXLOOP+1];
        miscinfo_Complex misc_complex;
        hairpin_tloop_Complex triloop_complex[MAXTRILOOPNO];
        //Complex nb_triloops_complex;
        hairpin_tloop_Complex tloop_complex[MAXTLOOPNO];
        //Complex nb_tloops;

//#if (MODEL == EXTENDED)        
        hairpin_tloop_Complex special_hl_complex[MAX_SPECIAL_LOOP_NO];
        Complex int11_experimental_addition_complex   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
        Complex int21_experimental_addition_complex   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
        Complex int22_experimental_addition_complex   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];        
        Complex internal_asymmetry_initiation_complex;
        Complex internal_asymmetry_slope_complex;
        Complex internal_asymmetry_complex [MAXLOOP+1];
        Complex bulgeA_complex;
        Complex bulgeC_complex;
        Complex bulgeG_complex;
        Complex bulgeU_complex;
        Complex bulge1_complex[NUCL] [NUCL] [NUCL] [NUCL] [NUCL];     // bulge of size 1     [i][j][k][ip][jp], where k=i+1, ip=k+1, j=jp+1
        
//#endif       
        
};



#endif


