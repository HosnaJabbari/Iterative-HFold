/***************************************************************************
                          s_specific_functions.cpp  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2002 by Mirela Andronescu
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


// This file contains functions specific to simfold

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_stacked_pair.h"
#include "s_specific_functions.h"
#include "s_min_folding.h"
#include "s_sub_folding.h"
#include "s_partition_function.h"


PARAMTYPE s_dangling_energy (int *sequence, char *structure, int i1, int i2, int i3, int i4)
//      (   )...(   )
//      i1..i2..i3..i4
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling energy between i2 and i3
// Feb 28, 2008. We might have a situation like this: <   >...(   ) or like this: (   )...<   >.
//  In that case, only add the parameter dangling onto the () pair, if it's at least 2 unpaired bases away
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    d_top = 0;
    d_bot = 0;

    // on Sep 5, 2008, I removed the minimization, so that 
    d_top = dangle_top[sequence[i2]] [sequence[i1]] [sequence[i2+1]];
    d_bot = dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]];

    //d_top = MIN (0, dangle_top[sequence[i2]] [sequence[i1]] [sequence[i2+1]]);
    //d_bot = MIN (0, dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]]);
    if (structure[i2] == '>' && structure[i3] == '(')   // pseudoknot, ignore dangling end dangling on it
    {
        if (i3 <= i2+2)     // >.( or >(   ignore completely
            energy = 0;
        else                // >...(
            energy = d_bot;
    }    
    else if (structure[i2] == ')' && structure[i3] == '<')   // pseudoknot, ignore dangling end dangling on it
    {
        if (i3 <= i2+2)     // ).< or )<   ignore completely
            energy = 0;
        else                // )...<
            energy = d_top;
    }
    else if (structure[i2] == '>' && structure[i3] == '<')  // case >..<  ignore completely
    {
        energy = 0;
    }                          
    else if (i2+1 == i3-1)     // see which is smaller
    {
        if (simple_dangling_ends)
            energy = d_top;
        else
            energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i2+1 < i3-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}


PARAMTYPE s_dangling_energy_left (int *sequence, char *structure, int i1, int i2, int i3, int i4)
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i1 and i3
// Feb 28, 2008. We might have a situation like this:  (....<    >   ). In that case, only add the 3' dangling end.
// If it's (.<    > ), don't add any
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    d_top = 0;
    d_bot = 0;            

    d_top = dangle_top[sequence[i1]] [sequence[i2]] [sequence[i1+1]];
    d_bot = dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]];
    
//     d_top = MIN (0, dangle_top[sequence[i1]] [sequence[i2]] [sequence[i1+1]]);
//     d_bot = MIN (0, dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]]);

    if (structure[i3] == '<')   // pseudoknot inside, ignore dangling end dangling on it
    {
        if (i3 <= i1+2)     // (< or (.<, ignore completely
            energy = 0;
        else                // (....<
            energy = d_top;
    }    
    else if (i1+1 == i3-1)     // see which is smaller
    {
        if (simple_dangling_ends)
            energy = d_top;
        else    
            energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i1+1 < i3-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}


PARAMTYPE s_dangling_energy_right (int *sequence, char *structure, int i1, int i2, int i3, int i4)
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2
// Feb 28, 2008. We might have a situation like this:  (    <    >...)
//  In that case, only add the 5' dangling end if it's at least two unpaired bases away
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    d_top = 0;

    d_bot = 0;

    d_top = dangle_top[sequence[i4]] [sequence[i3]] [sequence[i4+1]];
    d_bot = dangle_bot[sequence[i1]] [sequence[i2]] [sequence[i2-1]];
//     d_top = MIN (0, dangle_top[sequence[i4]] [sequence[i3]] [sequence[i4+1]]);
//     d_bot = MIN (0, dangle_bot[sequence[i1]] [sequence[i2]] [sequence[i2-1]]);

    if (structure[i4] == '>')   // pseudoknot inside, ignore dangling end dangling on it
    {
        if (i2 <= i4+2)     // >.) or >)   ignore completely
            energy = 0;
        else                // >...)
            energy = d_bot;
    }                              
    else if (i4+1 == i2-1)     // see which is smaller
    {
        if (simple_dangling_ends)
            energy = d_top;
        else    
            energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i4+1 < i2-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}


PARAMTYPE s_dangling_enthalpy (int *sequence, int i1, int i2, int i3, int i4)
// not in use
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling enthalpy between i2 and i3
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;

    d_top = MIN (0, enthalpy_dangle_top[sequence[i2]]
                      [sequence[i1]]
                      [sequence[i2+1]]);
    d_bot = MIN (0, enthalpy_dangle_bot[sequence[i4]]
                      [sequence[i3]]
                      [sequence[i3-1]]);
    if (i2+1 == i3-1)     // see which is smaller
    {
        energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i2+1 < i3-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}





PARAMTYPE s_dangling_enthalpy_left (int *sequence, int i1, int i2, int i3, int i4)
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling enthalpy between i1 and i3
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    d_top = 0;
    d_bot = 0;

    // this will be used in multi-loops.
    // add the dangle_top, even if it is positive
    d_top = MIN (0, enthalpy_dangle_top[sequence[i1]]
                          [sequence[i2]]
                          [sequence[i1+1]]);
    // in the other parts of the multi-loop, the dangles are added only if they are negative
    d_bot = MIN (0, enthalpy_dangle_bot[sequence[i4]]
                          [sequence[i3]]
                          [sequence[i3-1]]);
    if (i1+1 == i3-1)     // see which is smaller
    {
        energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i1+1 < i3-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}


PARAMTYPE s_dangling_enthalpy_right (int *sequence, int i1, int i2, int i3, int i4)
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling enthalpy between i4 and i2
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    d_top = 0;
    d_bot = 0;

    d_top = MIN (0, enthalpy_dangle_top[sequence[i4]]
                          [sequence[i3]]
                          [sequence[i4+1]]);
    d_bot = MIN (0, enthalpy_dangle_bot[sequence[i1]]
                          [sequence[i2]]
                          [sequence[i2-1]]);
    if (i4+1 == i2-1)     // see which is smaller
    {
        energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i4+1 < i2-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}



double free_energy_simfold (char *sequence, char *structure)
// PRE:  sequence and structure are given as input
// POST: returns the free energy of sequence folded into structure
//       Worst case complexity = n^2 (n = length of sequence)
{
	printf("in s_specific_functions.cpp -> free_energy_simfold \n");
    int i;
    str_features *f;
    double energy;
    int *int_sequence;
    PARAMTYPE en;
    int nb_nucleotides;
    
    nb_nucleotides = strlen(sequence); 

    if ((f = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");   
    // detect the structure features
    detect_structure_features (structure, f);
    
    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);
    en = s_calculate_energy (int_sequence, sequence, structure, f, NULL);
    energy = en/ 100.0;
    delete [] f;
    delete [] int_sequence;
    return energy;
}


double free_energy_simfold_restricted (char *sequence, char *structure, char *restricted)
// PRE:  sequence and structure are given as input
// POST: returns the free energy of sequence folded into structure
//       Worst case complexity = n^2 (n = length of sequence)
{
    int i;
    str_features *f;
    double energy;
    int *int_sequence;
    PARAMTYPE en;
    int nb_nucleotides;
        
    nb_nucleotides = strlen(sequence); 
    
    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");   
    // detect the structure features
    detect_structure_features (restricted, fres);            
    

    if ((f = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");   
    // detect the structure features
    detect_structure_features (structure, f);
    
    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);
    en = s_calculate_energy (int_sequence, sequence, structure, f, fres);
    energy = en/ 100.0;
    delete [] fres;
    delete [] f;
    delete [] int_sequence;
    return energy;
}


double enthalpy_simfold (char *sequence, char *structure)
// PRE:  sequence and structure are given as input
// POST: returns the enthalpy of sequence folded into structure
//       Worst case complexity = n^2 (n = length of sequence)
{
    int i;
    str_features *f;

    double energy;
    int *int_sequence;
    PARAMTYPE en;
    int nb_nucleotides;

    nb_nucleotides = strlen(sequence);
    if ((f = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (structure, f);

    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);
    en = s_calculate_enthalpy (int_sequence, sequence, f);
    energy = en/ 100.0;
    delete [] f;
    delete [] int_sequence;
    return energy;
}



PARAMTYPE s_calculate_energy (int *sequence, char *csequence, char *structure, str_features *f, str_features *fres)
// PRE:  the structure features have been determined
// POST: return the free energy (as integer)
// restricted is null if no restrictions
// restricted to be used for hairpin loops only
{
    int i; 
    PARAMTYPE energy, en, AUpen;
    PARAMTYPE dang;
    PARAMTYPE misc_energy;
    int h,l, nb_nucleotides;    
    static int cannot_add_dangling[MAXSLEN];
    int* ptable_restricted = NULL;
    
    nb_nucleotides = strlen (csequence);    
    for (i=0; i < nb_nucleotides; i++) cannot_add_dangling[i] = 0;

    if (fres != NULL)
    {
        ptable_restricted = new int[nb_nucleotides];
        for (i=0; i < nb_nucleotides; i++)
            ptable_restricted[i] = fres[i].pair;
    }

    energy = 0;
    AUpen = 0;
    
    for (i=0; i < nb_nucleotides; i++)
    {    
        // add some AU_penalties
        if (i==0 && f[i].pair > i)
        {
            AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
            if (debug)
                printf ("%d - AUpen1 \t- add energy %6d\n", i, AUpen);
            energy += AUpen;            
        }            
        else if ( i > 0 && f[i].pair > i && f[i-1].pair < i-1 &&
             f[i-1].pair != -1 && !cannot_add_dangling[i])
            //  )(  
        {            
            AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
            if (debug)
                printf ("%d - AUpen2 \t- add energy %6d\n", i, AUpen);
            energy += AUpen;            
        }            
    
        // add dangling energies and AU_penalties
        if (f[i].pair == -1 && !cannot_add_dangling[i])
        {
            if ((i == 0 || (i > 0 && f[i-1].pair == -1)) &&
                 i < nb_nucleotides-1 && f[i+1].pair > i+1)
                // .( or ..(
            {
                if (no_dangling_ends)
                    dang = 0;
                else
                {
                    // removed MIN (0, dangle) on Sep 5, 2008
                    dang = dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]];
                    //dang = MIN (0, dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]]);
                }
                AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
                if (debug)
                {
                    printf ("%d - dangle1 \t- add energy %6d\n", i, dang);
                    printf ("%d - AUpen3 \t- add energy %6d\n", i, AUpen);
                }                    
                energy += dang + AUpen;
            }
            else if ((i == nb_nucleotides-1 || 
                     (i < nb_nucleotides-1 && f[i+1].pair == -1)) &&
                     i > 0 && f[i-1].pair > -1 && f[i-1].pair < i-1)
                // ). or )..
            {
                if (no_dangling_ends)
                    dang = 0;
                else
                {
                    // removed MIN (0, dangle) on Sep 5, 2008
                    dang = dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]];
                    //dang = MIN (0, dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]]);
                }
                if (debug)                
                    printf ("%d - dangle2 \t- add energy %6d\n", i, dang);
                energy += dang;
            }
            else if (i < nb_nucleotides-1 && f[i+1].pair > i+1 && f[i-1].pair < i-1 && f[i-1].pair != -1)
               // ).( 
            {
                if (no_dangling_ends)
                    dang = 0;
                else
                {
                    // removed MIN (0, dangle) on Sep 5, 2008
                    dang = s_dangling_energy (sequence, structure, f[i-1].pair, i-1, i+1, f[i+1].pair);
                    //dang = MIN (0, s_dangling_energy (sequence, structure, f[i-1].pair, i-1, i+1, f[i+1].pair));
                }
                AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
                if (debug)
                {              
                    printf ("%d - dangle1 \t- add energy %6d\n", i, dang);
                    printf ("%d - AUpen4 \t- add energy %6d\n", i, AUpen);
                }    
                energy += dang + AUpen;
            }
            else
            {
                continue;
            }
        }
        
        if (f[i].pair < i)
        {
            continue;
        }       

        if (f[i].type == STACK)
        {
            en = s_stacked_pair::get_energy (i, f[i].pair, sequence);
            if (debug)            
                printf ("%d stack \t- add energy %6d\n", i, en);        
            energy += en;
        }
        else if (f[i].type == HAIRP)
        {
            en = s_hairpin_loop::get_energy (i, f[i].pair, sequence, csequence, ptable_restricted);
            // Mirela: added restriction in May 7th, 2005
            // Mar 24, 2008: I don't think I need this, it's dealt with inside s_hairpin_loop
            if (fres != NULL)
                if (fres[i].pair == f[i].pair && f[i].pair-i-1 < TURN)
                    en = 0;
            if (debug)
                printf ("%d hairpin \t- add energy %6d\n", i, en);
            energy += en;
        }

        else if (f[i].type == INTER)
        {
            int ip, jp;
            ip = f[i].bri[0];
            jp = f[f[i].bri[0]].pair;
            cannot_add_dangling[ip-1] = 1;
            cannot_add_dangling[jp+1] = 1;
            en = s_internal_loop::get_energy (i, f[i].pair, ip, jp, sequence, ptable_restricted);
            if (debug)
                printf ("%d internal \t- add energy %6d\n", i, en);        
            energy += en;
        }
        else  // (f[i].type == MULTI)
        {
            dang = 0;
            misc_energy = 0;
            AUpen = 0;
            int special;
            special = 0;
            // add the energies/enthalpies for free bases

            for (l=i+1; l < f[i].bri[0]; l++)
                misc_energy += misc.multi_free_base_penalty;
            for (h=0; h < f[i].num_branches-1; h++)
            {
                for (l = f[f[i].bri[h]].pair + 1; l < f[i].bri[h+1]; l++)
                    misc_energy += misc.multi_free_base_penalty;
            }
            for (l = f[f[i].bri[f[i].num_branches-1]].pair + 1; l < f[i].pair; l++)
                misc_energy += misc.multi_free_base_penalty;
                
            misc_energy += misc.multi_offset;
            misc_energy += misc.multi_helix_penalty * (f[i].num_branches + 1);

            // add AU_penalties for multi-loop
            AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);
            for (h=0; h < f[i].num_branches; h++)
                AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);
        
            // add dangling energies for multi-loop
            if (no_dangling_ends)
                dang = 0;
            else
            {
                dang += s_dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
                for (l=0; l < f[i].num_branches - 1; l++)
                    dang += s_dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
                dang += s_dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);
            }
            
            // add "no-dangling" restriction                                    
            for (l=0; l < f[i].num_branches; l++)
            {
                cannot_add_dangling [f[i].bri[l] -1] = 1;
                cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
            }
            if (debug)
            {
                printf ("%d - multi m\t- add energy %6d\n", i, misc_energy);            
                printf ("%d - multi d\t- add energy %6d\n", i, dang);
                printf ("%d - multi AU\t- add energy %6d\n", i, AUpen);
            }                
            energy += misc_energy + dang + AUpen;                           
        }        
    }
    if (fres != NULL)
    {
        delete [] ptable_restricted;
    }

    return energy;
}




PARAMTYPE s_calculate_enthalpy (int *sequence, char *csequence, str_features *f)
// PRE:  the structure features have been determined
// POST: return the enthalpy (as integer)
// to check!!!
{
    int i;
    PARAMTYPE energy, en, AUpen;
    PARAMTYPE dang;
    PARAMTYPE misc_energy;
    int h,l, nb_nucleotides;
    static int cannot_add_dangling[MAXSLEN];


    nb_nucleotides = strlen (csequence);
    for (i=0; i < nb_nucleotides; i++) cannot_add_dangling[i] = 0;

    energy = 0;
    AUpen = 0;

    for (i=0; i < nb_nucleotides; i++)
    {

        // add some AU_penalties
        if (i==0 && f[i].pair > i)
        {
            AUpen = AU_penalty_enthalpy (sequence[i], sequence[f[i].pair]);
            if (debug)
                printf ("%d - AUpen1 \t- add energy %6d\n", i, AUpen);
            energy += AUpen;
        }
        else if ( i > 0 && f[i].pair > i && f[i-1].pair < i-1 &&
             f[i-1].pair != -1 && !cannot_add_dangling[i])
            //  )(
        {
            AUpen = AU_penalty_enthalpy (sequence[i], sequence[f[i].pair]);
            if (debug)
                printf ("%d - AUpen2 \t- add energy %6d\n", i, AUpen);
            energy += AUpen;
        }

        // add dangling energies and AU_penalties
        if (f[i].pair == -1 && !cannot_add_dangling[i])
        {
            if ((i == 0 || (i > 0 && f[i-1].pair == -1)) &&
                 i < nb_nucleotides-1 && f[i+1].pair > i+1)
                // .( or ..(
            {
                dang = MIN (0, enthalpy_dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]]);
                AUpen = AU_penalty_enthalpy (sequence[i+1], sequence[f[i+1].pair]);
                if (debug)
                {
                    printf ("%d - dangle1 \t- add energy %6d\n", i, dang);
                    printf ("%d - AUpen3 \t- add energy %6d\n", i, AUpen);
                }
                energy += dang + AUpen;
            }
            else if ((i == nb_nucleotides-1 || 
                     (i < nb_nucleotides-1 && f[i+1].pair == -1 )) &&
                     i > 0 && f[i-1].pair > -1 && f[i-1].pair < i-1)
                // ). or )..
            {
                dang = MIN (0, enthalpy_dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]]);
                if (debug)
                    printf ("%d - dangle2 \t- add energy %6d\n", i, dang);
                energy += dang;
            }
            else if (i < nb_nucleotides-1 && f[i+1].pair > i+1 && f[i-1].pair < i-1 && f[i-1].pair != -1)
               // ).(
            {
                dang = MIN (0, s_dangling_enthalpy (sequence, f[i-1].pair, i-1, i+1, f[i+1].pair));
                AUpen = AU_penalty_enthalpy (sequence[i+1], sequence[f[i+1].pair]);
                if (debug)
                {
                    printf ("%d - dangle1 \t- add energy %6d\n", i, dang);
                    printf ("%d - AUpen4 \t- add energy %6d\n", i, AUpen);
                }
                energy += dang + AUpen;
            }
            else
            {
                continue;
            }
        }

        if (f[i].pair < i)
        {
            continue;
        }

        if (f[i].type == STACK)
        {
            en = s_stacked_pair::get_enthalpy (i, f[i].pair, sequence);
            if (debug)
                printf ("%d stack \t- add energy %6d\n", i, en);
            energy += en;
        }
        else if (f[i].type == HAIRP)
        {
            en = s_hairpin_loop::get_enthalpy (i, f[i].pair, sequence, csequence);
            if (debug)
                printf ("%d hairpin \t- add energy %6d\n", i, en);
            energy += en;
        }
        else if (f[i].type == INTER)
        {
            int ip, jp;
            ip = f[i].bri[0];
            jp = f[f[i].bri[0]].pair;
            cannot_add_dangling[ip-1] = 1;
            cannot_add_dangling[jp+1] = 1;
            en = s_internal_loop::get_enthalpy (i, f[i].pair, ip, jp, sequence);
            if (debug)
                printf ("%d internal \t- add energy %6d\n", i, en);
            energy += en;
        }
        else  // (f[i].type == MULTI)
        {
            dang = 0;
            misc_energy = 0;
            AUpen = 0;
            int special;
            special = 0;
            // add the energies/enthalpies for free bases
            // first find out if it is a regular multi-loop or a special multi-loop
            l = i;

//                printf ("Regular ML\n");
                for (l=i+1; l < f[i].bri[0]; l++)
                    misc_energy += enthalpy_misc.multi_free_base_penalty;
                for (h=0; h < f[i].num_branches-1; h++)
                {
                    for (l = f[f[i].bri[h]].pair + 1; l < f[i].bri[h+1]; l++)
                        misc_energy += enthalpy_misc.multi_free_base_penalty;
                }
                for (l = f[f[i].bri[f[i].num_branches-1]].pair + 1; l < f[i].pair; l++)
                    misc_energy += enthalpy_misc.multi_free_base_penalty;


                misc_energy += enthalpy_misc.multi_offset;
                misc_energy += enthalpy_misc.multi_helix_penalty * (f[i].num_branches + 1);
            
            // add AU_penalties for multi-loop
            AUpen += AU_penalty_enthalpy (sequence[i], sequence[f[i].pair]);
            for (h=0; h < f[i].num_branches; h++)
                AUpen += AU_penalty_enthalpy (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);

            // add dangling energies for multi-loop
            dang += s_dangling_enthalpy_left (sequence, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
            for (l=0; l < f[i].num_branches - 1; l++)
                dang += s_dangling_enthalpy (sequence, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
            dang += s_dangling_enthalpy_right (sequence, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);

            // add "no-dangling" restriction
            for (l=0; l < f[i].num_branches; l++)
            {
                cannot_add_dangling [f[i].bri[l] -1] = 1;
                cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
            }
            if (debug)
            {
                printf ("%d - multi m\t- add energy %6d\n", i, misc_energy);
                printf ("%d - multi d\t- add energy %6d\n", i, dang);
                printf ("%d - multi AU\t- add energy %6d\n", i, AUpen);
            }
            energy += misc_energy + dang + AUpen;
        }
    }
    return energy;
}



double simfold (char *sequence, char *structure)
// PRE:  the init_data function has been called;
//       the space for structure has been allocated
// POST: fold sequence, return the MFE structure in structure, and return the MFE
{
    double min_energy;
    s_min_folding *min_fold = new s_min_folding (sequence);
    min_energy = min_fold->s_simfold();
    min_fold->return_structure (structure);
    delete min_fold;
    return min_energy;
}    


double simfold_loss_augmented (char *sequence, char *known_structure, char *pred_structure)
// like simfold, but loss-augmented
// Returns arg min_y (dG(x,y) - loss(y,y*))
//  where x is sequence and y* is known_structure            
{
    double min_energy;
    s_min_folding *min_fold = new s_min_folding (sequence);
    int length = strlen (sequence);
    known_pairings = new int[length];
    pred_pairings = new int[length];
    detect_original_pairs (known_structure, known_pairings);
    min_energy = min_fold->s_simfold();
    min_fold->return_structure (pred_structure);
    delete [] known_pairings;
    delete [] pred_pairings;
    delete min_fold;
    return min_energy;
}    


double simfold_restricted (char *sequence, char *restricted, char *structure)
// PRE:  the init_data function has been called;
//       the space for structure has been allocated
// POST: fold sequence, return the MFE structure in structure, and return the MFE
{
    double min_energy;
    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted();
    min_fold->return_structure (structure);
    delete min_fold;
    // check if restricted is included in structure
    // ignore structures with infinite free energy - actually not necessarily
    //if (min_energy > 0)
    //    return INF;
    for (int i=0; i < strlen (sequence); i++)
    {
        if ((restricted[i] == '(' || restricted[i] == ')' || restricted[i] == '.') &&
            (restricted[i] != structure[i]))
        {
            printf ("ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
            printf ("  %s\n  %s\n  %s\t%.2lf\n", sequence, restricted, structure, min_energy);
            printf ("ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
            exit(1);
        }    
    }
    // now check if the free energy obtained with simfold_restricted is correct
    double correct_energy = free_energy_simfold_restricted (sequence, structure, restricted);
    if (fabs (correct_energy-min_energy) > 1.0)
    {
        printf ("ERROR!!! The dp energy is different from the energy calculated at the end!!\n");
        printf ("%s\n%s\n%s\n correct_energy=%.2lf, energy=%.2lf\n", sequence, restricted, structure, correct_energy, min_energy);
        printf ("ERROR!!! The dp energy is different from the energy calculated at the end!!\n");
        exit(1);
    }
    return min_energy;
}


int simfold_unordered_suboptimals_range (char *sequence, double max_energy, char structures[][MAXSLEN], double energies[], int num_subopt)
// Compute all suboptimal structures (limited by MAXSUBSTR), which have free energy <= max_energy
//     the free energies include the 5' and 3' dangling energies in all cases
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy;
    int energy_range;

    s_min_folding *min_fold = new s_min_folding (sequence);
    min_energy = min_fold->s_simfold();
    min_fold->return_structure (structure);
    delete min_fold;
                      
    if (min_energy >= max_energy)
    {
        strcpy (structures[0], structure);
        energies[0] = min_energy;    
        actual_num_str = 1;
    }
    else
    {
        energy_range = (int) ((max_energy - min_energy)*100);    
        //printf ("R1: %s %d\nR2: %s %d\n\n", sequence1, strlen(sequence1), sequence2, strlen(sequence2));
        s_sub_folding* sub_fold = new s_sub_folding(sequence, energy_range);
        sub_fold->set_limit(num_subopt);
        sub_fold->s_simfold (enthalpy);
        actual_num_str = sub_fold->return_structures(structures, energies);
        delete sub_fold;
    }    
    /*
    printf ("%s     %.2lf   %.2lf\n", sequence, max_energy, min_energy);
    for (int i=0; i < actual_num_str; i++)
    {
        double en = free_energy_simfold (sequence, structures[i]);
        if (en <= max_energy)
            printf ("%s    %.2lf\n", structures[i], energies[i]);
    }
    */
    return actual_num_str;      
}


int simfold_restricted_unordered_suboptimals_range (char *sequence, char *restricted, double max_energy, char structures[][MAXSLEN], double energies[])
// Compute all (restricted) suboptimal structures (limited by MAXSUBSTR), which have free energy <= max_energy
//     the free energies include the 5' and 3' dangling energies in all cases
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy;
    int energy_range;

    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted();
    min_fold->return_structure (structure);
    delete min_fold;

    if (min_energy >= max_energy)
    {
        strcpy (structures[0], structure);
        energies[0] = min_energy;    
        actual_num_str = 1;
    }
    else
    {
        energy_range = (int) ((max_energy - min_energy)*100);    
        //printf ("R1: %s %d\nR2: %s %d\n\n", sequence1, strlen(sequence1), sequence2, strlen(sequence2));
        s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, energy_range);
        sub_fold->set_limit(MAXSUBSTR);
        sub_fold->s_simfold_restricted (enthalpy);
        actual_num_str = sub_fold->return_structures(structures, energies);
        delete sub_fold;
    }    
    /*
    printf ("%s     %.2lf   %.2lf\n", sequence, max_energy, min_energy);
    for (int i=0; i < actual_num_str; i++)
    {
        double en = free_energy_simfold (sequence, structures[i]);
        if (en <= max_energy)
            printf ("%s    %.2lf\n", structures[i], energies[i]);
    }
    */
    return actual_num_str;      
}


int simfold_unordered_suboptimals (char *sequence, int number, char structures[][MAXSLEN], double energies[])
// compute "number" suboptimal structures, with dangling energies added everywhere
// return the number of suboptimal structures (between 1 and number)
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy;

    s_min_folding *min_fold = new s_min_folding (sequence);
    min_energy = min_fold->s_simfold();
    min_fold->return_structure (structure);

    delete min_fold;
      
    if (number == 1)
      {
        strcpy (structures[0], structure);
        energies[0] = min_energy;
      }
    else
      {    
        //printf ("R1: %s %d\nR2: %s %d\n\n", sequence1, strlen(sequence1), sequence2, strlen(sequence2));
        s_sub_folding* sub_fold = new s_sub_folding(sequence, 5000);
        sub_fold->set_limit(number);
        sub_fold->s_simfold (enthalpy);
        actual_num_str = sub_fold->return_structures(structures, energies);
        //        printf ("Actual unordered: %d\n", actual_num_str);
        delete sub_fold;
      }
    return actual_num_str;
      
}



int simfold_ordered_suboptimals (char *sequence, int number, char structures[][MAXSLEN], double energies[])
// compute "number" suboptimal structures, with dangling energies added everywhere
// in fact, compute 2*number structures, then compute the good free energy, reorder, and get the best "number" structures
// return the number of suboptimal structures (between 1 and number)
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy;
    int i;

    int positions[MAXSUBSTR];
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];

    if (number > MAXSUBSTR/2)
    {
        printf ("Desired number of suboptimal structures should be at most %d\n", MAXSUBSTR/2);
        exit(1);
    }

    s_min_folding *min_fold = new s_min_folding (sequence);
    min_energy = min_fold->s_simfold();
    min_fold->return_structure (structure);
	printf("in s_specific_functions.cpp: structure=%s \n",structure);

    delete min_fold;
      
    if (number == 1)
      {
        strcpy (structures[0], structure);
        energies[0] = min_energy;
      }
    else
      {    
        //printf ("R1: %s %d\nR2: %s %d\n\n", sequence1, strlen(sequence1), sequence2, strlen(sequence2));
        s_sub_folding* sub_fold = new s_sub_folding(sequence, 50000);        
        sub_fold->set_limit(2*number);
        sub_fold->s_simfold (enthalpy);
        actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
        printf("in s_specific_functions.cpp: sub_fold-> return_structures() DONE! \n");
		  
		printf("in s_specific_functions.cpp and actual_num_str = %d\n", actual_num_str);
        // added on Oct 16, 2007. The good free energies have to be recomputed here! For some reason it wasn't.
        for (i=0; i < actual_num_str; i++)
        {
            tmp_energies[i] = free_energy_simfold (sequence, tmp_structures[i]);
            //printf ("Unordered S %d: %s   %.2lf\n", i, tmp_structures[i], tmp_energies[i]);
        }
		printf ("in s_specific_functions.cpp: Actual unordered: %d\n", actual_num_str);
		  
        delete sub_fold;

        get_sorted_positions (actual_num_str, tmp_energies, positions);
		printf ("in s_specific_functions.cpp: get_sorted_positions() DONE!\n");  
		
        int MFE_there, ii;
        MFE_there = 0;
        ii = 0;
		  
	
        printf("in s_specific_functions.cpp, tmp_energies[positions[0]]=%d\n",tmp_energies[positions[0]]);
		  
        if (min_energy != tmp_energies[positions[0]])
          {
            //printf ("MFE1 %s    %.2lf \n =======\n", i+1, structure, energy);
            strcpy (structures[ii], structure);
            energies[ii] = min_energy;
            ii++;
          }
        else 
          {
            i = 0;
            while (tmp_energies[positions[i]] == min_energy)
              {
                if (strcmp (structure, tmp_structures[positions[i]]) == 0)
                  {
                    MFE_there = 1;
                    //printf ("MFE is subseq: ordered=%d, original=%d\n", i, positions[i]);
                    break;
                  }
                i++;
              }
            if (!MFE_there)
              {
                //printf ("MFE2 %s    %.2lf \n =======\n", i+1, structure, energy);
                strcpy (structures[ii], structure);
                energies[ii] = min_energy;
                ii++;
              }
          }

        // don't return more structures than number
        if (actual_num_str > number) actual_num_str = number;
        for (i=0; i < actual_num_str; i++)
          {
            //printf ("%2d  %s    %.2lf \n =======\n", i+1, tmp_structures[positions[i]], tmp_energies[positions[i]]);
            strcpy (structures[ii], tmp_structures[positions[i]]);
            energies[ii] = tmp_energies[positions[i]];
            ii++;
          }                
      }
    return actual_num_str;      
}



int simfold_restricted_unordered_suboptimals (char *sequence, char *restricted, int number, char structures[][MAXSLEN], double energies[])
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy;

    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted ();
    min_fold->return_structure (structure);

    delete min_fold;
    printf ("Structure: %s\n", structure);
      
    if (number == 1)
      {
        strcpy (structures[0], structure);
        energies[0] = min_energy;
        actual_num_str=1;
      }
    else
      {    
        //printf ("R1: %s %d\nR2: %s %d\n\n", sequence1, strlen(sequence1), sequence2, strlen(sequence2));
        s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, 5000);
        sub_fold->set_limit(number);
        sub_fold->s_simfold_restricted (enthalpy);
        actual_num_str = sub_fold->return_structures(structures, energies);
        //printf ("Actual unordered: %d\n", actual_num_str);
        delete sub_fold;
      }
    return actual_num_str;      
}



int simfold_restricted_all_suboptimals (char *sequence, char *restricted, char structures[][MAXSLEN], double energies[])
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy;

    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted ();
    min_fold->return_structure (structure);

    delete min_fold;      
    s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, -(int)(min_energy*100.0));
    sub_fold->set_limit(MAXSUBSTR);
    sub_fold->s_simfold_restricted (enthalpy);
    actual_num_str = sub_fold->return_structures(structures, energies);
    //printf ("Actual unordered: %d\n", actual_num_str);
    delete sub_fold;
    return actual_num_str;      
}

int simfold_restricted_all_mfe_structures (char *sequence, char *restricted, char structures[][MAXSLEN], double energies[])
// return MAXSUBSTR suboptimal structures, which have free energy <= 0. Note that the suboptimal structure computation adds both 5' and 3' dangling ends in any case. The suboptimal structures are not reordered in this function.
{
    int actual_num_str;
    char structure[MAXSLEN];
    double min_energy, enthalpy, energy;
    char tmp_structures[10][MAXSLEN];
    double tmp_energies[10];
    
    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted ();
    min_fold->return_structure (structure);
    delete min_fold;   
    //printf ("Res mfe:%s\t%.2lf\n", structure, min_energy);
    // generate the first 10, and recompute free energy
    s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, 1000);
    sub_fold->set_limit(10);
    sub_fold->s_simfold_restricted (enthalpy);
    actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
    int j = 0;
    for (int i=0; i < actual_num_str; i++)
    {        
        energy = free_energy_simfold_restricted (sequence, tmp_structures[i], restricted);
        if (energy == min_energy)
        {
            strcpy (structures[j], tmp_structures[i]);
            energies[j] = tmp_energies[i];
            j++;
        }
        //printf ("%d\t%s\t%.2lf\t%.2lf\n", i, tmp_structures[i], tmp_energies[i], energy);
    }
    if (j==0)    // if nothing was added, add the mfe structure
    {    
        strcpy (structures[j], structure);        
        energies[j] = min_energy;
        j = 1;        
    }
    //printf ("Actual unordered: %d\n", actual_num_str);
    delete sub_fold;
    return j;      
}


PFTYPE simfold_partition_function_smart (char *sequence, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
{
    PFTYPE pf;
    int i,j;
    s_partition_function *part = new s_partition_function (sequence, ignore_dangles, compute_gradient_dangles);
    pf = part->compute_partition_function();
    //part->print_u();
    //part->compute_base_pair_probabilities();
    //part->PS_dot_plot("dot.ps");
    delete part;
    //printf ("ignore_dangles=%d, pf=%Lg\n", ignore_dangles, pf);
    return pf;
}


void simfold_gradient_smart (char *sequence, PFTYPE *grad, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
{
    PFTYPE pf;
    int i,j;
    s_partition_function *part = new s_partition_function (sequence, ignore_dangles, compute_gradient_dangles);
    pf = part->compute_partition_function();   
    part->compute_base_pair_probabilities();
    part->compute_logZ_gradient();
    part->copy_gradient(grad);
    //part->PS_dot_plot("dot.ps");
    delete part;
}

PFTYPE simfold_f_and_gradient_smart (char *sequence, char *restricted, PFTYPE *grad, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
{
    PFTYPE pf;
    int i,j;
    s_partition_function *part = new s_partition_function (sequence, ignore_dangles, compute_gradient_dangles, restricted);
    
    //pf = part->compute_partition_function_exhaustively();   
    //part->compute_logZ_gradient_exhaustively();
    //part->copy_gradient(grad);
    
    pf = part->compute_partition_function();
    part->compute_base_pair_probabilities();
    part->compute_logZ_gradient();
    part->copy_gradient(grad);

    /*
    for (i=260; i <= 307; i++)
    {
        if (grad[i] != 0)
        {
            printf ("grad[%d]=%g, %s\n", i, grad[i], sequence);
        }
    }
    */
    //part->PS_dot_plot("dot.ps");
    delete part;
    return pf;
}

void simfold_partition_function_both (char *sequence)
{
    //long double pf, pfexhaust;
    int i,j;
    s_partition_function *part = new s_partition_function (sequence);
    part->compute_partition_function();   
    part->compute_base_pair_probabilities();

    part->compute_partition_function_exhaustively();   
    part->verify_partition_function ();
    //part->PS_dot_plot("dot.ps");
    delete part;
}

PFTYPE simfold_partition_function_approximately (char *sequence)
{
    char structure[MAXSLEN];
    double enthalpy, energy;
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];
    double min_energy, max_energy;
    int actual_num_str;
    PFTYPE Z, Zexact;
    int i,ii,j, k;
    double R, temp, beta;
    //R = 0.00198717;
    //temp = 310.15;
    beta = 1000.0/(1.98717*310.15);
    //beta = 1; 

    simple_dangling_ends = 1;
    int seqlen = strlen(sequence);
    double papp[seqlen][seqlen];
    for(i=0; i < seqlen; i++)
    {
        for (j=i+1; i < seqlen; i++)
            papp[i][j] = 0;
    }    
    Zexact = simfold_partition_function_smart (sequence);
    int ptable[MAXSLEN];
    PFTYPE strprob;
        
    s_min_folding *min_fold = new s_min_folding (sequence);
    min_energy = min_fold->s_simfold ();
    min_fold->return_structure (structure);
    delete min_fold;      
    
    //s_sub_folding* sub_fold = new s_sub_folding(sequence, -(int)(min_energy*100.0));
    s_sub_folding* sub_fold = new s_sub_folding(sequence, 10000);
    sub_fold->set_limit (MAXSUBSTR);
    sub_fold->s_simfold (enthalpy);
    actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
    delete sub_fold;
    
    max_energy = tmp_energies[actual_num_str-1];
    Z = 0.0;
    for (i=0; i < actual_num_str; i++)
    {
        // recompute the free energy, i.e. with the correct dangling ends
        energy = free_energy_simfold (sequence, tmp_structures[i]);
        strprob = (PFTYPE)(exp(-1.0 * energy * beta)/Zexact);
        detect_original_pairs (tmp_structures[i], ptable); 
        for (ii=0; ii < seqlen; ii++)
        {
            if (ptable[ii] > ii)    
            {
                //printf ("ptable[%d]=%d\n", ii, ptable[ii]);
                papp[ii][ptable[ii]] += strprob;
            }
        }   
        
        printf ("Substr %d: %s\ten=%.2lf \tprob=%Le\n", i, tmp_structures[i], energy, strprob);
        Z += exp (-1.0 * energy * beta);
    }
    for(i=0; i < seqlen; i++)
    {
        for (j=i+1; j < seqlen; j++)
            if (papp[i][j] > 0)
                printf ("papp[%d][%d] = %g\n", i, j, sqrt(papp[i][j]));
    }    

    /*    
    // try to get the exhaustive up[i,j]
    double up = 0;
    int l;
    for (i=0; i < seqlen; i++)
    {
        for (j=i+TURN+1; j < seqlen; j++)
        {    
            up = 0;
            if (can_pair(nuc_to_int(sequence[i]), nuc_to_int(sequence[j])))
            {
                // compute exhaustive up[i,j]
                for (k=0; k < actual_num_str; k++)                
                {
                    if (has_base_pair(i,j, tmp_structures[k]))
                    {
                        int identical = 0;
                        for (l=0; l < k; l++)
                        {
                            if (identical_structure (i, j, tmp_structures[l], tmp_structures[k]))
                                identical = 1;
                        }
                        if (!identical)
                            up += exp_free_energy_partial (i, j, sequence, tmp_structures[k]);
                        //else
                        //    printf ("Str %s identical with other for i=%d, j=%d\n", tmp_structures[k], i, j);
                    }                        
                }
                printf ("upexhaust(%d,%d) = %g\n", i, j, up);
            }
        }
    }
    */
    return Z;    
}



// for now, the following functions exist only if the parameters are [long] double 

#ifdef INCLUDE_FUNCTIONS
#include "s_partition_function_complex.h"

PFTYPE simfold_partition_function_smart_numerical (char *sequence, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
{
    Complex pf;
    int i,j;
    s_partition_function_complex *part = new s_partition_function_complex (sequence, ignore_dangles, compute_gradient_dangles);
    pf = part->compute_partition_function();
    //part->print_u();
    //part->compute_base_pair_probabilities();
    //part->PS_dot_plot("dot.ps");
    delete part;
    //printf ("ignore_dangles=%d, pf=%Lg\n", ignore_dangles, pf);
    return pf.real();
}

PFTYPE simfold_partition_function_smart_double (char *sequence, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
{
    PFTYPE pf;
    int i,j;
    s_partition_function_complex *part = new s_partition_function_complex (sequence, ignore_dangles, compute_gradient_dangles, 1);
    pf = part->compute_partition_function().real();
    //part->print_u();
    //part->compute_base_pair_probabilities();
    //part->PS_dot_plot("dot.ps");
    delete part;
    //printf ("ignore_dangles=%d, pf=%Lg\n", ignore_dangles, pf);
    return pf;
}

void simfold_gradient_numerical (char *sequence, PFTYPE *grad, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
{
    Complex pf;
    int i,j;
    s_partition_function_complex *part = new s_partition_function_complex (sequence, ignore_dangles, compute_gradient_dangles);
    pf = part->compute_partition_function();   
    part->compute_logZ_gradient_numerical();
    part->copy_gradient_numerical(grad);
    //part->PS_dot_plot("dot.ps");
    delete part;
}


PFTYPE simfold_f_and_gradient_smart_double (char *sequence, PFTYPE *grad, int ignore_dangles, int compute_gradient_dangles)
// unrestricted for now
// Sep 11, 2007: all parameters are double
{
    PFTYPE pf;
    int i,j;
    s_partition_function_complex *part = new s_partition_function_complex (sequence, ignore_dangles, compute_gradient_dangles, 1);
    Complex pfc = part->compute_partition_function();
    pf = pfc.real();
    //printf ("partition function = %lf\n", pf);
    part->compute_base_pair_probabilities();
    part->compute_logZ_gradient();
    part->copy_gradient(grad);
    if (!part->correct_gradient_nan())
    {
        printf ("Nan gradients\n");
        part->print_gradient();
    }
    /*
    for (i=260; i <= 307; i++)
    {
        if (grad[i] != 0)
        {
            printf ("grad[%d]=%g, %s\n", i, grad[i], sequence);
        }
    }
    */
    //part->PS_dot_plot("dot.ps");
    delete part;
    return pf;
}


PFTYPE simfold_f_and_gradient_smart_numerical (char *sequence, PFTYPE *grad, int ignore_dangles,
    int compute_gradient_dangles, int which_param)    
// unrestricted for now
// if which_param is -1, compute gradient for all params, otherwise just for which_param
{
    Complex pf;
    int i,j;
    s_partition_function_complex *part = new s_partition_function_complex (sequence, ignore_dangles, compute_gradient_dangles);
    pf = part->compute_partition_function();
    if (which_param == -1)
        part->compute_logZ_gradient_numerical();
    else
        part->compute_logZ_gradient_numerical(which_param);
    part->copy_gradient_numerical(grad);
    delete part;
    return pf.real();
}


#endif



