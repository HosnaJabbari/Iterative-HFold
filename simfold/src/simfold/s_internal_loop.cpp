/***************************************************************************
                          s_internal_loop.cpp  -  description
                             -------------------
    begin                : Fri Apr 12 2002
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

// a class for internal loop related functions
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "s_internal_loop.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"                             
#include "params.h"


s_internal_loop::s_internal_loop (int *seq, int length)
// The constructor
{
    seqlen = length;
    sequence = seq;
    this->V = NULL;
}


s_internal_loop::~s_internal_loop ()
// The destructor
{
}


PARAMTYPE s_internal_loop::compute_energy (int i, int j)
// computes the MFE of the structure closed by an internal loop closed at (i,j)
{

    // TODO
    //return 0;
    //printf ("\n1\n");
    
    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en, specialen; 
    int branch1, branch2, l;
    ttmp = INF;
    mmin = INF;
    //specialen = 0;
    ip_jp_energy = 0;
    i_j_energy = 0;
    
    for (ip = i+1; ip <= MIN(j-2-TURN,i+MAXLOOP+1) ; ip++)  // j-2-TURN
    {
		// Hosna, August 28, 2012
		// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
		// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
		// (____(_____((__(_______)__))_______________________________)____)
		// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
		// in this example int(5,59,11,27) is falsely missed and is equal to INF
		// So I am changing it to the be jp=ip+1+TURN; jp<j; jp++ instead
        //minq = MAX (j-i+ip-MAXLOOP-2, ip+1+TURN);    // ip+1+TURN);
		minq = ip+1 + TURN;
		
        for (jp = minq; jp < j; jp++)
        {        
            // could replace the code from now to the end of the function, by get_energy_str
            // it's very redundant (duplicated code) as it is now, but it's faster.
            if (sequence[ip]+sequence[jp] == 3 ||
                sequence[ip]+sequence[jp] == 5)        
            {
            
                branch1 = ip-i-1;
                branch2 = j-jp-1;
                if (branch1 == 0 && branch2 == 0)
                    continue;

                if (branch1 != 0 || branch2 != 0)
                {
                    // check if it is a bulge loop of size 1
                    // check if it is int11 or int21 or int22
                    // 1x1 internal loop
                    if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
                    {
                        // int11[i][j][i+1][j-1][ip][jp]
                        en = int11 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]];
                        ttmp = en + V->get_energy (ip, jp);
                    }
                    // 1x2 internal loop
                    else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
                    {
                        // int21[i][j][i+1][j-1][ip][jp][jp+1]
                        en = int21 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]]
                                   [sequence[jp+1]];
                        ttmp = en + V->get_energy (ip, jp);
                    }
                    // 1x2 internal loop
                    else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
                    {
                        // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
                        en = int21 [sequence[jp]][sequence[ip]]
                                   [sequence[j-1]][sequence[ip-1]]
                                   [sequence[j]][sequence[i]]
                                   [sequence[i+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        //printf ("1: ttmp=%g, mmin=%g\n", ttmp, mmin);
                    }
                    // 2x2 internal loop
                    else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
                    {
                        // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
                        en = int22 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]]
                                   [sequence[ip-1]][sequence[jp+1]];
                        ttmp = en + V->get_energy (ip, jp);
                    }
                    else
                    {
                        // this case is not int11, int21, int22
                        
                        // check if it is a bulge
                        if (branch1 == 0 || branch2 == 0)
                        {
                            l = branch1+branch2;
                            penalty_size = penalty_by_size (l, 'B');
                            if (l == 1)
                            {
                                // bulge of size 1
                                // stack[i][j][i+1][j-1]
                                if (parsi_bulge1 == T99)
                                {
                                    en = stack [sequence[i]][sequence[j]]
                                            [sequence[ip]][sequence[jp]];                                
                                }
                                else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
                                {                                
                                    int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
                                    if (branch1 == 1)
                                    {
                                        i2 = sequence[i];
                                        j2 = sequence[j];
                                        k2 = sequence[i+1];
                                        ip2 = sequence[ip];
                                        jp2 = sequence[jp];
                                    }
                                    else        // it's upside down
                                    {
                                        i2 = sequence[jp];
                                        j2 = sequence[ip];
                                        k2 = sequence[j-1];
                                        ip2 = sequence[j];
                                        jp2 = sequence[i];
                                    }
                                    en = bulge1[i2][j2][k2][ip2][jp2];
                                    penalty_size = 0;   // we don't add it for case 1
                                }
                                ttmp = en + penalty_size + V->get_energy (ip, jp);
                            }
                            else
                            {
                                // bulge of size bigger than 1
                                // check if (i,j) and (ip,jp) can pair
                                ttmp = penalty_size +
                                    AU_penalty (sequence[i],sequence[j]) +
                                    AU_penalty (sequence[ip], sequence[jp]) +
                                    V->get_energy (ip, jp);
                            }
                        }
                        // it is a regular internal loop (not a bulge)
                        else
                        {
                            l = branch1+branch2;
                            penalty_size = penalty_by_size (l, 'I');
                            asym_penalty = asymmetry_penalty (branch1, branch2);

                            if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                            // If gail_rule is set to 1 in miscloop file,
                            // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                            {
//#if (MODEL == SIMPLE)
                                // In the simple model I only use 3 parameters for tstacki, 
                                //  So tstacki[i][j][[0][0] never comes up, I just ignore it 
                                //i_j_energy  =  tstacki[sequence[i]][sequence[j]] [0][0];
                                //ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [0][0];
//#elif (MODEL == EXTENDED)
                                if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                                      ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                                {
                                    i_j_energy = misc.internal_AU_closure;
                                }
                                if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                                      ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                                {
                                    ip_jp_energy = misc.internal_AU_closure;
                                }
//#endif                                
                            }
                            else
                            {
                                i_j_energy   = tstacki[sequence[i]][sequence[j]]
                                                      [sequence[i+1]][sequence[j-1]];
                                ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
                                                      [sequence[jp+1]][sequence[ip-1]];
                                //specialen = special_energy_internal (sequence, i,j,ip,jp);
                                i_j_energy += special_energy_internal (sequence, i,j,ip,jp);
                            }
                            ttmp = i_j_energy + ip_jp_energy + penalty_size + // specialen +
                                         asym_penalty + V->get_energy (ip, jp);

                        }
                    }
                }
            }
            
            // add the loss
            if (pred_pairings != NULL && ttmp < INF/2)
            {
                pred_pairings[i] = j;     pred_pairings[j] = i;
                pred_pairings[ip] = jp;   pred_pairings[jp] = ip;
                for (int kk=i+1; kk < ip; kk++) pred_pairings[kk] = -1;
                for (int kk=jp+1; kk < j; kk++) pred_pairings[kk] = -1;
                ttmp = ttmp - loss (i,ip-1) - loss (jp+1,j);
            }            

            if (ttmp < mmin)
            {
                mmin = ttmp;
            }       
        }
    }
    return mmin;
}


PARAMTYPE s_internal_loop::compute_energy_restricted (int i, int j, str_features *fres)
// computes the MFE of the structure closed by a restricted internal loop closed by (i,j)
{

    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;

    for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there
    {
        // Hosna, August 28, 2012
		// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
		// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
		// (____(_____((__(_______)__))_______________________________)____)
		// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
		// in this example int(5,59,11,27) is falsely missed and is equal to INF
		// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
        //minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
		minq = ip+1;		
        for (jp = minq; jp < j; jp++)
        {        
            if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
                continue;
            //ttmp = get_energy_str (i, j, ip, jp);
			// Hosna, March 26, 2012
			// changed to accommodate non-canonical base pairs in the restricted structure
			ttmp = get_energy_str_restricted (i, j, ip, jp, fres);

            if (ttmp < mmin)
            {
                mmin = ttmp;
            }
        }
    }
    return mmin;
}

//TODO: PMO
PARAMTYPE s_internal_loop::compute_energy_restricted_pmo (int i, int j, str_features *fres)
// computes the MFE of the structure closed by a restricted internal loop closed by (i,j)
{
    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;

    for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there
    {
        // Hosna, August 28, 2012
		// TODO: cannot understand why we have the following calculations, as it makes the following case be missed!
		// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
		// (____(_____((__(_______)__))_______________________________)____)
		// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
		// in this example int(5,59,11,27) is falsely missed and is equal to INF
		// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
        //minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
		minq = ip+1;		
        for (jp = minq; jp < j; jp++)
        {        
            if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
                continue;
            //ttmp = get_energy_str (i, j, ip, jp);
			// Hosna, March 26, 2012
			// changed to accommodate non-canonical base pairs in the restricted structure
			ttmp = get_energy_str_restricted_pmo (i, j, ip, jp, fres);
			
            if (ttmp < mmin) {
                mmin = ttmp;
            }
        }
    }

    return mmin;
}

PARAMTYPE s_internal_loop::compute_energy_restricted_pkonly (int i, int j, str_features *fres)
// computes the MFE of the structure closed by a restricted internal loop closed by (i,j)
{
	
    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;
	
	// Hosna, August 31, 2012
	// The following restriction misses the long restricted loops, so I am chaning it
    //for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there
	for (ip = i+1; ip <= j-2 ; ip++)  // the -TURN shouldn't be there
    {
		// Hosna, August 28, 2012
		// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
		// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
		// (____(_____((__(_______)__))_______________________________)____)
		// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
		// in this example int(5,59,11,27) is falsely missed and is equal to INF
		// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
        //minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
		minq = ip+1;
        for (jp = minq; jp < j; jp++)
        {        
            if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres)){
				continue;
			}
            //ttmp = get_energy_str (i, j, ip, jp);
			// pkonly version
			ttmp = (fres[ip].pair == jp && fres[jp].pair == ip)? get_energy_str_restricted (i, j, ip, jp, fres): INF;
			
            if (ttmp < mmin)
            {
                mmin = ttmp;
            }
        }
    }
    return mmin;
}

//TODO: PK
PARAMTYPE s_internal_loop::compute_energy_restricted_pkonly_pmo (int i, int j, str_features *fres)
// computes the MFE of the structure closed by a restricted internal loop closed by (i,j)
{
	
    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;
	
	// Hosna, August 31, 2012
	// The following restriction misses the long restricted loops, so I am chaning it
    //for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there
	for (ip = i+1; ip <= j-2 ; ip++)  // the -TURN shouldn't be there
    {
		// Hosna, August 28, 2012
		// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
		// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
		// (____(_____((__(_______)__))_______________________________)____)
		// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
		// in this example int(5,59,11,27) is falsely missed and is equal to INF
		// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
        //minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
		minq = ip+1;
        for (jp = minq; jp < j; jp++) {        
			if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres)){
				continue;
			}
            //ttmp = get_energy_str (i, j, ip, jp);
			// pkonly version
			ttmp = (fres[ip].pair == jp && fres[jp].pair == ip)? get_energy_str_restricted_pmo (i, j, ip, jp, fres): INF;
            if (ttmp < mmin) {
                mmin = ttmp;
            }
        }
    }
    return mmin;
}

PARAMTYPE s_internal_loop::get_energy_str_restricted (int i, int j, int ip, int jp, str_features *fres)
// returns the free energy of the structure closed by the internal loop (i,j,ip,jp)
// This function is just most of what is inside the second for loop of compute_energy

// Hosna, March 26, 2012
// cases that need to be considered in this function are:
// 1) i.j and ip.jp are enforced in the restricted structure and all are non-canonical
// 2) i.j is enforced in the restricted structure and we have i.j non-canonical but ip.jp is canonical
// 3) ip.jp is enforced in the restricted structure and we have i.j canonical but ip.jp is non-canonical
{
    // TODO
    //return 0;
    //printf ("\n2\n");
    
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;
    i_j_energy = 0; ip_jp_energy = 0;
    
	if ((sequence[ip]+sequence[jp] == 3 || sequence[ip]+sequence[jp] == 5) && can_pair(sequence[i],sequence[j])) // normal case        
	{
		
		branch1 = ip-i-1;
		branch2 = j-jp-1;
		
		if (branch1 != 0 || branch2 != 0)
		{
			// check if it is a bulge loop of size 1
			// check if it is int11 or int21 or int22
			if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
			{
				// int11[i][j][i+1][j-1][ip][jp]
				en = int11 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]];
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
			{
				// int21[i][j][i+1][j-1][ip][jp][jp+1]
				en = int21 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]]
				[sequence[jp+1]];
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
			{
				// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
				en = int21 [sequence[jp]][sequence[ip]]
				[sequence[j-1]][sequence[ip-1]]
				[sequence[j]][sequence[i]]
				[sequence[i+1]];
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
			{
				// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
				en = int22 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]]
				[sequence[ip-1]][sequence[jp+1]];
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else
			{
				// this case is not int11, int21, int22
				
				// check if it is a bulge
				if (branch1 == 0 || branch2 == 0)
				{
					l = branch1+branch2;
					penalty_size = penalty_by_size (l, 'B');
					if (l == 1)
					{
						// bulge of size 1
						// stack[i][j][i+1][j-1]
						if (parsi_bulge1 == T99)
						{                                    
							en = stack [sequence[i]][sequence[j]]
							[sequence[ip]][sequence[jp]];                                
						}
						else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
						{                                
							int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
							if (branch1 == 1)
							{
								i2 = sequence[i];
								j2 = sequence[j];
								k2 = sequence[i+1];
								ip2 = sequence[ip];
								jp2 = sequence[jp];
							}
							else        // it's upside down
							{
								i2 = sequence[jp];
								j2 = sequence[ip];
								k2 = sequence[j-1];
								ip2 = sequence[j];
								jp2 = sequence[i];
							}
							en = bulge1[i2][j2][k2][ip2][jp2];
							penalty_size = 0;   // we don't add it for case 1
						}
						
						ttmp = en + penalty_size + V->get_energy (ip, jp);
						if (ttmp < mmin)
						{
							mmin = ttmp;
						}
					}
					else
					{
						// bulge of size bigger than 1
						// check if (i,j) and (ip,jp) can pair
						ttmp = penalty_size +
						AU_penalty (sequence[i],sequence[j]) +
						AU_penalty (sequence[ip], sequence[jp]) +
						V->get_energy (ip, jp);
						if (ttmp < mmin)
						{
							mmin = ttmp;
						}
					}
				}
				// it is an internal loop (not a bulge)
				else
				{
					l = branch1+branch2;
					penalty_size = penalty_by_size (l, 'I');
					asym_penalty = asymmetry_penalty (branch1, branch2);
					
					if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
						// If gail_rule is set to 1 in miscloop file,
						// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
					{
						//#if (MODEL == SIMPLE)
						// In the simple model I only use 3 parameters for tstacki, 
						//  So tstacki[i][j][[0][0] never comes up, I just ignore it 
						//i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
						//ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
						//#elif (MODEL == EXTENDED)
						if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
							((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
						{
							i_j_energy = misc.internal_AU_closure;
						}
						if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
							((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
						{
							ip_jp_energy = misc.internal_AU_closure;
						}                                
						//#endif
					}
					else
					{
						i_j_energy   = tstacki[sequence[i]][sequence[j]]
						[sequence[i+1]][sequence[j-1]];
						
						ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
						[sequence[jp+1]][sequence[ip-1]];
						
						i_j_energy += special_energy_internal (sequence, i,j,ip,jp);
						
					}
					ttmp = i_j_energy + ip_jp_energy + penalty_size +
					asym_penalty + V->get_energy (ip, jp);
					
					
					if (ttmp < mmin)
					{
						mmin = ttmp;
					}
				}
			}
		}
	}  else if((fres[ip].pair==jp && fres[jp].pair==ip && can_pair(sequence[i],sequence[j])) || 
			   (fres[ip].pair==jp && fres[jp].pair==ip && fres[i].pair==j && fres[j].pair==i) ||
			   (fres[i].pair==j && fres[j].pair==i && can_pair(sequence[ip],sequence[jp]))){ 		
		branch1 = ip-i-1;
		branch2 = j-jp-1;
		
		if (branch1 != 0 || branch2 != 0)
		{
			// check if it is a bulge loop of size 1
			// check if it is int11 or int21 or int22
			if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
			{
				// int11[i][j][i+1][j-1][ip][jp]
				en = int11 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]];
				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
			{
				// int21[i][j][i+1][j-1][ip][jp][jp+1]
				en = int21 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]]
				[sequence[jp+1]];
				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
			{
				// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
				en = int21 [sequence[jp]][sequence[ip]]
				[sequence[j-1]][sequence[ip-1]]
				[sequence[j]][sequence[i]]
				[sequence[i+1]];
				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
			{
				// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
				en = int22 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]]
				[sequence[ip-1]][sequence[jp+1]];
				
				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else
			{
				// this case is not int11, int21, int22
				
				// check if it is a bulge
				if (branch1 == 0 || branch2 == 0)
				{
					l = branch1+branch2;
					penalty_size = penalty_by_size (l, 'B');
					if (l == 1)
					{
						// bulge of size 1
						// stack[i][j][i+1][j-1]
						if (parsi_bulge1 == T99)
						{                                    
							en = stack [sequence[i]][sequence[j]]
							[sequence[ip]][sequence[jp]]; 
							en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
						}
						else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
						{                                
							int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
							if (branch1 == 1)
							{
								i2 = sequence[i];
								j2 = sequence[j];
								k2 = sequence[i+1];
								ip2 = sequence[ip];
								jp2 = sequence[jp];
							}
							else        // it's upside down
							{
								i2 = sequence[jp];
								j2 = sequence[ip];
								k2 = sequence[j-1];
								ip2 = sequence[j];
								jp2 = sequence[i];
							}
							en = bulge1[i2][j2][k2][ip2][jp2];
							en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
							penalty_size = 0;   // we don't add it for case 1
						}
						ttmp = en + penalty_size + V->get_energy (ip, jp);
						
						if (ttmp < mmin)
						{
							mmin = ttmp;
						}
					}
					else
					{
						// bulge of size bigger than 1
						// check if (i,j) and (ip,jp) can pair
						ttmp = penalty_size +
						AU_penalty (sequence[i],sequence[j]) +
						AU_penalty (sequence[ip], sequence[jp]) +
						V->get_energy (ip, jp);
						if (ttmp < mmin)
						{
							mmin = ttmp;
						}
					}
				}
				// it is an internal loop (not a bulge)
				
				// Hosna, March 26, 2012
				// the following parts may need modification to accommodate non-canonical base pairs in restricted structure
				else
				{
					l = branch1+branch2;
					penalty_size = penalty_by_size (l, 'I');
					asym_penalty = asymmetry_penalty (branch1, branch2);
					
					if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
						// If gail_rule is set to 1 in miscloop file,
						// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
					{
						//#if (MODEL == SIMPLE)
						// In the simple model I only use 3 parameters for tstacki, 
						//  So tstacki[i][j][[0][0] never comes up, I just ignore it 
						//i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
						//ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
						//#elif (MODEL == EXTENDED)
						if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
							((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
						{
							i_j_energy = misc.internal_AU_closure;
						}
						if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
							((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
						{
							ip_jp_energy = misc.internal_AU_closure;
						}                                
						//#endif
					}
					else // Hosna, June 4, 2012, fixed the non-canonical base pairing problem as folllows
					{
						if (can_pair(sequence[i],sequence[j])){
							i_j_energy   = tstacki[sequence[i]][sequence[j]]
												[sequence[i+1]][sequence[j-1]];
						}else{
							i_j_energy=0;
						}
						if(can_pair(sequence[ip],sequence[jp])){
							ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
												[sequence[jp+1]][sequence[ip-1]];
						}else{
							ip_jp_energy =0;
						}
						
						i_j_energy += special_energy_internal (sequence, i,j,ip,jp);
					}
					ttmp = i_j_energy + ip_jp_energy + penalty_size +
					asym_penalty + V->get_energy (ip, jp);
					
					if (ttmp < mmin)
					{
						mmin = ttmp;
					}
				}
			}
		}
		
		
	}
    if (mmin < INF/2)
    {
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = j;     pred_pairings[j] = i;
            pred_pairings[ip] = jp;   pred_pairings[jp] = ip;
            for (int kk=i+1; kk < ip; kk++) pred_pairings[kk] = -1;
            for (int kk=jp+1; kk < j; kk++) pred_pairings[kk] = -1;
            mmin = mmin - loss (i,ip-1) - loss (jp+1,j);
        }            
		
        return mmin;
    }
    return INF;            
}

PARAMTYPE s_internal_loop::get_energy_str_restricted_pmo (int i, int j, int ip, int jp, str_features *fres)
// returns the free energy of the structure closed by the internal loop (i,j,ip,jp)
// This function is just most of what is inside the second for loop of compute_energy

// Hosna, March 26, 2012
// cases that need to be considered in this function are:
// 1) i.j and ip.jp are enforced in the restricted structure and all are non-canonical
// 2) i.j is enforced in the restricted structure and we have i.j non-canonical but ip.jp is canonical
// 3) ip.jp is enforced in the restricted structure and we have i.j canonical but ip.jp is non-canonical
{
    // TODO
    //return 0;
    //printf ("\n2\n");
    
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
	double pmo_percentage, rna_percentage, pmo_percentage_p, rna_percentage_p, pmo_percentage_avg, rna_percentage_avg;
    mmin = INF;
    i_j_energy = 0; ip_jp_energy = 0;
    
	// Compute percente usages
	get_pmo_usage_percentages(i, j, &pmo_percentage, &rna_percentage);
	get_pmo_usage_percentages(ip, jp, &pmo_percentage_p, &rna_percentage_p);
	pmo_percentage_avg = (pmo_percentage + pmo_percentage_p)/2;
	rna_percentage_avg = (rna_percentage + rna_percentage_p)/2;

	pmo_percentage_avg = pmo_percentage;
	rna_percentage_avg = rna_percentage;

	if (sequence[i] == 4 || sequence[j] == 4 || sequence[ip] == 4 || sequence[jp] == 4)
		return 0;
	if (sequence[i+1] == 4 || sequence[j+1] == 4 || sequence[ip+1] == 4 || sequence[jp+1] == 4 || sequence[i-1] == 4 || sequence[j-1] == 4 || sequence[ip-1] == 4 || sequence[jp-1] == 4)
		return INF;

	if ((sequence[ip]+sequence[jp] == 3 || sequence[ip]+sequence[jp] == 5) && can_pair(sequence[i],sequence[j])) // normal case        
	{
		
		branch1 = ip-i-1;
		branch2 = j-jp-1;
		
		if (branch1 != 0 || branch2 != 0)
		{
			// check if it is a bulge loop of size 1
			// check if it is int11 or int21 or int22
			if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
			{
				// int11[i][j][i+1][j-1][ip][jp]
				en = (PARAMTYPE) round(pmo_percentage_avg*int11_pmo[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]] + 
									rna_percentage_avg*int11[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]);
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
			{
				// int21[i][j][i+1][j-1][ip][jp][jp+1]
				en = (PARAMTYPE) round(pmo_percentage_avg*int21_pmo[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[jp+1]] + 
									rna_percentage_avg*int21[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[jp+1]]);
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
			{
				// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
				en = (PARAMTYPE) round(pmo_percentage_avg*int21_pmo[sequence[jp]]
															[sequence[ip]]
															[sequence[j-1]]
															[sequence[ip-1]]
															[sequence[j]]
															[sequence[i]]
															[sequence[i+1]] + 
									rna_percentage_avg*int21[sequence[jp]]
															[sequence[ip]]
															[sequence[j-1]]
															[sequence[ip-1]]
															[sequence[j]]
															[sequence[i]]
															[sequence[i+1]]);
				
				ttmp = en + V->get_energy (ip, jp);			

				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
			{
				// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
				en = (PARAMTYPE) round(pmo_percentage_avg*int22_pmo[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[ip-1]]
															[sequence[jp+1]] + 
									rna_percentage_avg*int22[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[ip-1]]
															[sequence[jp+1]]);
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin)
				{
					mmin = ttmp;
				}
			}
			else
			{
				// this case is not int11, int21, int22
				
				// check if it is a bulge
				if (branch1 == 0 || branch2 == 0)
				{
					l = branch1+branch2;
					penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'B') + rna_percentage_avg*penalty_by_size(l, 'B'));
					if (l == 1)
					{
						// bulge of size 1
						// stack[i][j][i+1][j-1]
						if (parsi_bulge1 == T99)
						{                                    
							en = (PARAMTYPE) round(pmo_percentage_avg*stack_pmo[sequence[i]]
																		[sequence[j]]
																		[sequence[ip]]
																		[sequence[jp]] + 
												rna_percentage_avg*stack[sequence[i]]
																		[sequence[j]]
																		[sequence[ip]]
																		[sequence[jp]]);
						}
						else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
						{                                
							int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
							if (branch1 == 1)
							{
								i2 = sequence[i];
								j2 = sequence[j];
								k2 = sequence[i+1];
								ip2 = sequence[ip];
								jp2 = sequence[jp];
							}
							else        // it's upside down
							{
								i2 = sequence[jp];
								j2 = sequence[ip];
								k2 = sequence[j-1];
								ip2 = sequence[j];
								jp2 = sequence[i];
							}
							en = (PARAMTYPE) round(pmo_percentage_avg*bulge1_pmo[i2][j2][k2][ip2][jp2] + rna_percentage_avg*bulge1[i2][j2][k2][ip2][jp2]);
							penalty_size = 0;   // we don't add it for case 1
						}
						
						ttmp = en + penalty_size + V->get_energy (ip, jp);
						if (ttmp < mmin) {
							mmin = ttmp;
						}
					}
					else
					{
						// bulge of size bigger than 1
						// check if (i,j) and (ip,jp) can pair
						ttmp = penalty_size + 
								(PARAMTYPE) round(pmo_percentage*AU_penalty_pmo(sequence[i],sequence[j]) + rna_percentage*AU_penalty(sequence[i],sequence[j])) + 
								(PARAMTYPE) round(pmo_percentage_p*AU_penalty_pmo(sequence[ip], sequence[jp]) + rna_percentage_p*AU_penalty(sequence[ip], sequence[jp])) + 
								V->get_energy (ip, jp);
						
						if (ttmp < mmin) {
							mmin = ttmp;
						}
					}
				}
				// it is an internal loop (not a bulge)
				else
				{
					l = branch1+branch2;
					penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'I') + rna_percentage_avg*penalty_by_size(l, 'I'));
					asym_penalty = (PARAMTYPE) round(pmo_percentage_avg*asymmetry_penalty_pmo(branch1, branch2) + rna_percentage_avg*asymmetry_penalty(branch1, branch2));
					asym_penalty = asymmetry_penalty(branch1, branch2);
					if ((branch1 == 1 || branch2 == 1) && (misc.gail_rule || misc_pmo.gail_rule))
						// If gail_rule is set to 1 in miscloop file,
						// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
					{
						//#if (MODEL == SIMPLE)
						// In the simple model I only use 3 parameters for tstacki, 
						//  So tstacki[i][j][[0][0] never comes up, I just ignore it 
						//i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
						//ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
						//#elif (MODEL == EXTENDED)
						if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
							((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
						{
							i_j_energy = (PARAMTYPE) round(pmo_percentage*misc_pmo.internal_AU_closure + rna_percentage*misc.internal_AU_closure);
						}
						if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
							((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
						{
							ip_jp_energy = (PARAMTYPE) round(pmo_percentage_p*misc_pmo.internal_AU_closure + rna_percentage_p*misc.internal_AU_closure);
						}                                
						//#endif
					}
					else
					{
						i_j_energy  = (PARAMTYPE) round(pmo_percentage*tstacki_pmo[sequence[i]]
																			[sequence[j]]
																			[sequence[i+1]]
																			[sequence[j-1]] + 
														rna_percentage*tstacki[sequence[i]]
																			[sequence[j]]
																			[sequence[i+1]]
																			[sequence[j-1]]);
						
						ip_jp_energy = (PARAMTYPE) round(pmo_percentage_p*tstacki_pmo[sequence[jp]]
																				[sequence[ip]]
																				[sequence[jp+1]]
																				[sequence[ip-1]] + 
														rna_percentage_p*tstacki[sequence[jp]]
																				[sequence[ip]]
																				[sequence[jp+1]]
																				[sequence[ip-1]]);
						
						i_j_energy += (PARAMTYPE) round(pmo_percentage_avg*special_energy_internal_pmo(sequence, i,j,ip,jp) + rna_percentage_avg*special_energy_internal(sequence, i,j,ip,jp));
						
					}
					ttmp = i_j_energy + ip_jp_energy + penalty_size +
					asym_penalty + V->get_energy (ip, jp);
					
					
					if (ttmp < mmin) {
						mmin = ttmp;
					}
				}
			}
		}
	}  else if((fres[ip].pair==jp && fres[jp].pair==ip && can_pair(sequence[i],sequence[j])) || 
			   (fres[ip].pair==jp && fres[jp].pair==ip && fres[i].pair==j && fres[j].pair==i) ||
			   (fres[i].pair==j && fres[j].pair==i && can_pair(sequence[ip],sequence[jp]))){ 		
		branch1 = ip-i-1;
		branch2 = j-jp-1;
		
		if (branch1 != 0 || branch2 != 0)
		{
			// check if it is a bulge loop of size 1
			// check if it is int11 or int21 or int22
			if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
			{
				// int11[i][j][i+1][j-1][ip][jp]
				en = (PARAMTYPE) round(pmo_percentage_avg*int11_pmo[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]] + 
									rna_percentage_avg*int11[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]);
				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
			{
				// int21[i][j][i+1][j-1][ip][jp][jp+1]
				en = (PARAMTYPE) round(pmo_percentage_avg*int21_pmo[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[jp+1]] + 
									rna_percentage_avg*int21[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[jp+1]]);

				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
			{
				// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
				en = (PARAMTYPE) round(pmo_percentage_avg*int21_pmo[sequence[jp]]
															[sequence[ip]]
															[sequence[j-1]]
															[sequence[ip-1]]
															[sequence[j]]
															[sequence[i]]
															[sequence[i+1]] + 
									rna_percentage_avg*int21[sequence[jp]]
															[sequence[ip]]
															[sequence[j-1]]
															[sequence[ip-1]]
															[sequence[j]]
															[sequence[i]]
															[sequence[i+1]]);

				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
			{
				// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
				en = (PARAMTYPE) round(pmo_percentage_avg*int22_pmo[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[ip-1]]
															[sequence[jp+1]] + 
									rna_percentage_avg*int22[sequence[i]]
															[sequence[j]]
															[sequence[i+1]]
															[sequence[j-1]]
															[sequence[ip]]
															[sequence[jp]]
															[sequence[ip-1]]
															[sequence[jp+1]]);
				
				en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
				
				ttmp = en + V->get_energy (ip, jp);
				
				if (ttmp < mmin) {
					mmin = ttmp;
				}
			}
			else
			{
				// this case is not int11, int21, int22
				
				// check if it is a bulge
				if (branch1 == 0 || branch2 == 0)
				{
					l = branch1+branch2;
					penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'B') + rna_percentage_avg*penalty_by_size(l, 'B'));
					if (l == 1)
					{
						// bulge of size 1
						// stack[i][j][i+1][j-1]
						if (parsi_bulge1 == T99)
						{                                    
							en = (PARAMTYPE) round(pmo_percentage_avg*stack_pmo[sequence[i]]
																		[sequence[j]]
																		[sequence[ip]]
																		[sequence[jp]] + 
												rna_percentage_avg*stack[sequence[i]]
																		[sequence[j]]
																		[sequence[ip]]
																		[sequence[jp]]); 

							en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
						}
						else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
						{                                
							int i2, j2, k2, ip2, jp2;
							penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'I') + rna_percentage_avg*penalty_by_size(l, 'I'));
							asym_penalty = (PARAMTYPE) round(pmo_percentage_avg*asymmetry_penalty_pmo(branch1, branch2) + rna_percentage_avg*asymmetry_penalty(branch1, branch2));//  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
							if (branch1 == 1)
							{
								i2 = sequence[i];
								j2 = sequence[j];
								k2 = sequence[i+1];
								ip2 = sequence[ip];
								jp2 = sequence[jp];
							}
							else        // it's upside down
							{
								i2 = sequence[jp];
								j2 = sequence[ip];
								k2 = sequence[j-1];
								ip2 = sequence[j];
								jp2 = sequence[i];
							}
							en = (PARAMTYPE) round(pmo_percentage_avg*bulge1_pmo[i2][j2][k2][ip2][jp2] + rna_percentage_avg*bulge1[i2][j2][k2][ip2][jp2]);
							en = MIN(0,en); // Hosna, March 26, 2012, added to accommodate non-cannonical base pairs in restricted structure
							penalty_size = 0;   // we don't add it for case 1
						}
						ttmp = en + penalty_size + V->get_energy (ip, jp);
						
						if (ttmp < mmin) {
							mmin = ttmp;
						}
					}
					else
					{
						// bulge of size bigger than 1
						// check if (i,j) and (ip,jp) can pair
						ttmp = penalty_size + 
							(PARAMTYPE) round(pmo_percentage*AU_penalty_pmo(sequence[i],sequence[j]) + rna_percentage*AU_penalty(sequence[i],sequence[j])) + 
							(PARAMTYPE) round(pmo_percentage_p*AU_penalty_pmo(sequence[ip], sequence[jp]) + rna_percentage_p*AU_penalty(sequence[ip], sequence[jp])) + 
							V->get_energy (ip, jp);
						
						if (ttmp < mmin) {
							mmin = ttmp;
						}
					}
				}
				// it is an internal loop (not a bulge)
				
				// Hosna, March 26, 2012
				// the following parts may need modification to accommodate non-canonical base pairs in restricted structure
				else
				{
					l = branch1+branch2;
					penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'I') + rna_percentage_avg*penalty_by_size(l, 'I'));
					asym_penalty = (PARAMTYPE) round(pmo_percentage_avg*asymmetry_penalty_pmo(branch1, branch2) + rna_percentage_avg*asymmetry_penalty(branch1, branch2));
					
					if ((branch1 == 1 || branch2 == 1) && (misc.gail_rule || misc_pmo.gail_rule))
						// If gail_rule is set to 1 in miscloop file,
						// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
					{
						//#if (MODEL == SIMPLE)
						// In the simple model I only use 3 parameters for tstacki, 
						//  So tstacki[i][j][[0][0] never comes up, I just ignore it 
						//i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
						//ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
						//#elif (MODEL == EXTENDED)
						if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
							((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
						{
							i_j_energy = (PARAMTYPE) round(pmo_percentage*misc_pmo.internal_AU_closure + rna_percentage*misc.internal_AU_closure);
						}
						if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
							((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
						{
							ip_jp_energy = (PARAMTYPE) round(pmo_percentage_p*misc_pmo.internal_AU_closure + rna_percentage_p*misc.internal_AU_closure);
						}                                
						//#endif
					}
					else // Hosna, June 4, 2012, fixed the non-canonical base pairing problem as folllows
					{
						if (can_pair(sequence[i],sequence[j])){
							i_j_energy = (PARAMTYPE) round(pmo_percentage*tstacki_pmo[sequence[i]]
																				[sequence[j]]
																				[sequence[i+1]]
																				[sequence[j-1]] +
															rna_percentage*tstacki[sequence[i]]
																				[sequence[j]]
																				[sequence[i+1]]
																				[sequence[j-1]]);

						}else{
							i_j_energy=0;
						}
						if(can_pair(sequence[ip],sequence[jp])){
							ip_jp_energy = (PARAMTYPE) round(pmo_percentage_p*tstacki_pmo[sequence[jp]]
																					[sequence[ip]]
																					[sequence[jp+1]]
																					[sequence[ip-1]] + 
															rna_percentage_p*tstacki[sequence[jp]]
																					[sequence[ip]]
																					[sequence[jp+1]]
																					[sequence[ip-1]]);
						}else{
							ip_jp_energy =0;
						}
						
						i_j_energy += (PARAMTYPE) round(pmo_percentage_avg*special_energy_internal_pmo(sequence, i,j,ip,jp) + rna_percentage_avg*special_energy_internal(sequence, i,j,ip,jp));
					}
					ttmp = i_j_energy + ip_jp_energy + penalty_size +
					asym_penalty + V->get_energy (ip, jp);
					
					if (ttmp < mmin) {
						mmin = ttmp;
					}
				}
			}
		}
		
		
	}
    if (mmin < INF/2)
    {
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = j;     
			pred_pairings[j] = i;
            pred_pairings[ip] = jp;   
			pred_pairings[jp] = ip;
            for (int kk=i+1; kk < ip; kk++) { pred_pairings[kk] = -1; }
            for (int kk=jp+1; kk < j; kk++) { pred_pairings[kk] = -1; }
            mmin = mmin - loss (i,ip-1) - loss (jp+1,j);
        }               

        return mmin;
    }
    return INF;            
}


PARAMTYPE s_internal_loop::get_energy_str (int i, int j, int ip, int jp)
// returns the free energy of the structure closed by the internal loop (i,j,ip,jp)
// This function is just most of what is inside the second for loop of compute_energy
{
    // TODO
    //return 0;
    //printf ("\n2\n");
    
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;
    i_j_energy = 0; ip_jp_energy = 0;
    
            if (sequence[ip]+sequence[jp] == 3 || sequence[ip]+sequence[jp] == 5 )        
            {
            
                branch1 = ip-i-1;
                branch2 = j-jp-1;

                if (branch1 != 0 || branch2 != 0)
                {
                    // check if it is a bulge loop of size 1
                    // check if it is int11 or int21 or int22
                    if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
                    {
                        // int11[i][j][i+1][j-1][ip][jp]
                        en = int11 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]];
						
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin)
                        {
                            mmin = ttmp;
                        }
                    }
                    else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
                    {
                        // int21[i][j][i+1][j-1][ip][jp][jp+1]
                        en = int21 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]]
                                   [sequence[jp+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin)
                        {
                            mmin = ttmp;
                        }
                    }
                    else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
                    {
                        // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
                        en = int21 [sequence[jp]][sequence[ip]]
                                   [sequence[j-1]][sequence[ip-1]]
                                   [sequence[j]][sequence[i]]
                                   [sequence[i+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin)
                        {
                            mmin = ttmp;
                        }
                    }
                    else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
                    {
                        // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
                        en = int22 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]]
                                   [sequence[ip-1]][sequence[jp+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin)
                        {
                            mmin = ttmp;
                        }
                    }
                    else
                    {
                        // this case is not int11, int21, int22
                        
                        // check if it is a bulge
                        if (branch1 == 0 || branch2 == 0)
                        {
                            l = branch1+branch2;
                            penalty_size = penalty_by_size (l, 'B');
                            if (l == 1)
                            {
                                // bulge of size 1
                                // stack[i][j][i+1][j-1]
                                if (parsi_bulge1 == T99)
                                {                                    
                                    en = stack [sequence[i]][sequence[j]]
                                            [sequence[ip]][sequence[jp]];                                
                                }
                                else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
                                {                                
                                    int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
                                    if (branch1 == 1)
                                    {
                                        i2 = sequence[i];
                                        j2 = sequence[j];
                                        k2 = sequence[i+1];
                                        ip2 = sequence[ip];
                                        jp2 = sequence[jp];
                                    }
                                    else        // it's upside down
                                    {
                                        i2 = sequence[jp];
                                        j2 = sequence[ip];
                                        k2 = sequence[j-1];
                                        ip2 = sequence[j];
                                        jp2 = sequence[i];
                                    }
                                    en = bulge1[i2][j2][k2][ip2][jp2];
                                    penalty_size = 0;   // we don't add it for case 1
                                }
                                ttmp = en + penalty_size + V->get_energy (ip, jp);
                                if (ttmp < mmin)
                                {
                                    mmin = ttmp;
                                }
                            }
                            else
                            {
                                // bulge of size bigger than 1
                                // check if (i,j) and (ip,jp) can pair
                                ttmp = penalty_size +
                                    AU_penalty (sequence[i],sequence[j]) +
                                    AU_penalty (sequence[ip], sequence[jp]) +
                                    V->get_energy (ip, jp);
                                if (ttmp < mmin)
                                {
                                    mmin = ttmp;
                                }
                            }
                        }
                        // it is an internal loop (not a bulge)
                        else
                        {
                            l = branch1+branch2;
                            penalty_size = penalty_by_size (l, 'I');
                            asym_penalty = asymmetry_penalty (branch1, branch2);

                            if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                            // If gail_rule is set to 1 in miscloop file,
                            // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                            {
//#if (MODEL == SIMPLE)
                                // In the simple model I only use 3 parameters for tstacki, 
                                //  So tstacki[i][j][[0][0] never comes up, I just ignore it 
                                //i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
                                //ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
//#elif (MODEL == EXTENDED)
                                if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                                      ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                                {
                                    i_j_energy = misc.internal_AU_closure;
                                }
                                if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                                      ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                                {
                                    ip_jp_energy = misc.internal_AU_closure;
                                }                                
//#endif
                            }
                            else
                            {
                                i_j_energy   = tstacki[sequence[i]][sequence[j]]
                                                      [sequence[i+1]][sequence[j-1]];
                                ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
                                                      [sequence[jp+1]][sequence[ip-1]];
                                i_j_energy += special_energy_internal (sequence, i,j,ip,jp);
                            }
                            ttmp = i_j_energy + ip_jp_energy + penalty_size +
                                         asym_penalty + V->get_energy (ip, jp);
                            if (ttmp < mmin)
                            {
                                mmin = ttmp;
                            }
                        }
                    }
                }
            }    
    if (mmin < INF/2)
    {
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = j;     pred_pairings[j] = i;
            pred_pairings[ip] = jp;   pred_pairings[jp] = ip;
            for (int kk=i+1; kk < ip; kk++) pred_pairings[kk] = -1;
            for (int kk=jp+1; kk < j; kk++) pred_pairings[kk] = -1;
            mmin = mmin - loss (i,ip-1) - loss (jp+1,j);
        }            

        return mmin;
    }
    return INF;            
}


// Not used. Tried to see if the gradient in s_partition_function gets computed faster
//  if you make a different function for each situation.
// It doesn't seem to be the case.
PARAMTYPE s_internal_loop::get_energy_00 (int i, int j, int ip, int jp, int *sequence)
{
    PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
    int branch1, branch2, l;

    branch1 = ip-i-1;
    branch2 = j-jp-1;
    l = branch1+branch2;
    penalty_size = penalty_by_size (l, 'I');
    asym_penalty = asymmetry_penalty (branch1, branch2);

    i_j_energy  =  tstacki[sequence[i]][sequence[j]]
                            [0][0];
    ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
                            [0][0];
    energy = i_j_energy + ip_jp_energy + penalty_size +
                            asym_penalty;
    return energy;
}                            
                            

PARAMTYPE s_internal_loop::get_energy (int i, int j, int ip, int jp, int *sequence, int *ptable)
// returns the free energy of the internal loop closed at (i,j,ip,jp)
// Also called by the partition function
{
    // TODO
    //return 100;
    //printf ("\n3\n");

    PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
    int branch1, branch2, l;

    if (exists_restricted_ptable (i,ip,ptable) || exists_restricted_ptable (jp,j,ptable))
        return INF;

    branch1 = ip-i-1;
    branch2 = j-jp-1;

    if (branch1 != 0 || branch2 != 0)
    {
        // check if it is a bulge loop of size 1
        // check if it is int11 or int21 or int22
        if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
        {
                        // int11[i][j][i+1][j-1][ip][jp]
            energy = int11 [sequence[i]][sequence[j]]
                           [sequence[i+1]][sequence[j-1]]
                           [sequence[ip]][sequence[jp]];
            check_int11_parameters (sequence[i], sequence[j], sequence[i+1], sequence[j-1], sequence[ip], sequence[jp]);
        }
        else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
        {
            // int21[i][j][i+1][j-1][ip][jp][jp+1]
            energy = int21 [sequence[i]][sequence[j]]
                           [sequence[i+1]][sequence[j-1]]
                           [sequence[ip]][sequence[jp]]
                           [sequence[jp+1]];
        }
        else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
        {
            // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
            energy = int21 [sequence[jp]][sequence[ip]]
                           [sequence[j-1]][sequence[ip-1]]
                           [sequence[j]][sequence[i]]
                           [sequence[i+1]];
        }
        else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
        {
            // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
            energy = int22 [sequence[i]][sequence[j]]
                           [sequence[i+1]][sequence[j-1]]
                           [sequence[ip]][sequence[jp]]
                           [sequence[ip-1]][sequence[jp+1]];
            //printf ("IN GET_ENERGY int22(%d,%d,%d,%d,%d,%d,%d,%d): %g\n", sequence[i],sequence[j],sequence[i+1],sequence[j-1],sequence[ip],sequence[jp],sequence[ip-1],sequence[jp+1],energy);
        }
        else
        {
            // this case is not int11, int21, int22

            // check if it is a bulge
            if (branch1 == 0 || branch2 == 0)
            {
                l = branch1+branch2;
                penalty_size = penalty_by_size (l, 'B');
                if (l == 1)
                {
                    // bulge of size 1
                    // stack[i][j][i+1][j-1]
                    if (parsi_bulge1 == T99)
                    {
                        energy = stack [sequence[i]][sequence[j]]
                                    [sequence[ip]][sequence[jp]] +
                                penalty_size;
                    }
                    else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
                    {
                        int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
                        if (branch1 == 1)
                        {
                            i2 = sequence[i];
                            j2 = sequence[j];
                            k2 = sequence[i+1];
                            ip2 = sequence[ip];
                            jp2 = sequence[jp];
                        }
                        else        // it's upside down
                        {
                            i2 = sequence[jp];
                            j2 = sequence[ip];
                            k2 = sequence[j-1];
                            ip2 = sequence[j];
                            jp2 = sequence[i];
                        }
                        energy = bulge1[i2][j2][k2][ip2][jp2];
                    }
                }
                else
                {
                    // bulge of size bigger than 1
                    // check if (i,j) and (ip,jp) can pair
                    //printf ("REAL len = %d, pen = %Lg\n", l, penalty_size);
                    energy = penalty_size +
                               AU_penalty (sequence[i],sequence[j]) +
                               AU_penalty (sequence[ip], sequence[jp]);
                }
            }
            // it is an internal loop (not a bulge)
            else
            {
                l = branch1+branch2;
                penalty_size = penalty_by_size (l, 'I');
                asym_penalty = asymmetry_penalty (branch1, branch2);
                //if (l == 5  && abs(branch1-branch2) == 3)
                //    printf ("REAL pen = %Lg, asym = %Lf\n", penalty_size, asym_penalty);

                
                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                // If gail_rule is set to 1 in miscloop file,
                // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                {
                    i_j_energy = 0;
                    ip_jp_energy = 0;

//#if (MODEL == SIMPLE)
                    // In the simple model I only use 3 parameters for tstacki, 
                    //  So tstacki[i][j][[0][0] never comes up, I just ignore it                    
                    //i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
                    //ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
//#elif (MODEL == EXTENDED)                    
                    if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                          ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                    {
                        i_j_energy = misc.internal_AU_closure;
                    }
                    if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                          ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                    {
                        ip_jp_energy = misc.internal_AU_closure;
                    }    
//#endif                                    
                    //printf ("IN GET_ENERGY, energy is %lf\n", i_j_energy + ip_jp_energy);
                }
                else
                {
                    i_j_energy   = tstacki[sequence[i]][sequence[j]]
                                          [sequence[i+1]][sequence[j-1]];
                    ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
                                          [sequence[jp+1]][sequence[ip-1]];
                    i_j_energy += special_energy_internal (sequence, i,j,ip,jp);
                }
                energy = i_j_energy + ip_jp_energy + penalty_size +
                                         asym_penalty;
            }
        }
    }
    //printf ("REAL:   i=%d, j=%d, ip=%d, jp=%d, energy=%g\n", i, j, ip, jp, energy);
    return energy;
}

PARAMTYPE s_internal_loop::get_energy_pmo (int i, int j, int ip, int jp, int *sequence, int *ptable)
// returns the free energy of the internal loop closed at (i,j,ip,jp)
// Also called by the partition function
{
    // TODO
    //return 100;
    //printf ("\n3\n");

    PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
    int branch1, branch2, l;
	double pmo_percentage, rna_percentage, pmo_percentage_p, rna_percentage_p, pmo_percentage_avg, rna_percentage_avg;
    
	// Compute percent usages
	get_pmo_usage_percentages(i, j, &pmo_percentage, &rna_percentage);
	get_pmo_usage_percentages(ip, jp, &pmo_percentage_p, &rna_percentage_p);
	pmo_percentage_avg = (pmo_percentage + pmo_percentage_p)/2;
	rna_percentage_avg = (rna_percentage + rna_percentage_p)/2;

	pmo_percentage_avg = pmo_percentage;
	rna_percentage_avg = rna_percentage;
	pmo_percentage_p = pmo_percentage;
	rna_percentage_p = rna_percentage;

	if (sequence[i] == 4 || sequence[j] == 4 || sequence[ip] == 4 || sequence[jp] == 4 || sequence[i+1] == 4 || sequence[j+1] == 4 || sequence[ip+1] == 4 || sequence[jp+1] == 4 || sequence[i-1] == 4 || sequence[j-1] == 4 || sequence[ip-1] == 4 || sequence[jp-1] == 4)
		return INF;

    if (exists_restricted_ptable (i,ip,ptable) || exists_restricted_ptable (jp,j,ptable))
        return INF;

    branch1 = ip-i-1;
    branch2 = j-jp-1;

    if (branch1 != 0 || branch2 != 0)
    {
        // check if it is a bulge loop of size 1
        // check if it is int11 or int21 or int22
        if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
        {
                        // int11[i][j][i+1][j-1][ip][jp]
            energy = (PARAMTYPE) round(pmo_percentage_avg*int11_pmo[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]][sequence[ip]][sequence[jp]] + rna_percentage_avg*int11[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]][sequence[ip]][sequence[jp]]);
			
            check_int11_pmo_parameters (sequence[i], sequence[j], sequence[i+1], sequence[j-1], sequence[ip], sequence[jp]);
        }
        else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
        {
            // int21[i][j][i+1][j-1][ip][jp][jp+1]
            energy = (PARAMTYPE) round(pmo_percentage_avg*int21_pmo[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]][sequence[ip]][sequence[jp]][sequence[jp+1]] + rna_percentage_avg*int21[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]][sequence[ip]][sequence[jp]][sequence[jp+1]]);
        }
        else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
        {
            // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
            energy = (PARAMTYPE) round(pmo_percentage_avg*int21_pmo[sequence[jp]][sequence[ip]][sequence[j-1]][sequence[ip-1]][sequence[j]][sequence[i]][sequence[i+1]] + rna_percentage_avg*int21[sequence[jp]][sequence[ip]][sequence[j-1]][sequence[ip-1]][sequence[j]][sequence[i]][sequence[i+1]]);
        }
        else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
        {
            // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
            energy = (PARAMTYPE) round(pmo_percentage_avg*int22_pmo[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]][sequence[ip]][sequence[jp]][sequence[ip-1]][sequence[jp+1]] + rna_percentage_avg*int22[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]][sequence[ip]][sequence[jp]][sequence[ip-1]][sequence[jp+1]]);
        }
        else
        {
            // this case is not int11, int21, int22

            // check if it is a bulge
            if (branch1 == 0 || branch2 == 0)
            {
                l = branch1+branch2;
                penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'B') + rna_percentage_avg*penalty_by_size(l, 'B'));
                if (l == 1)
                {
                    // bulge of size 1
                    // stack[i][j][i+1][j-1]
                    if (parsi_bulge1 == T99)
                    {
                        energy = (PARAMTYPE) round(pmo_percentage_avg*stack_pmo[sequence[i]][sequence[j]][sequence[ip]][sequence[jp]] + rna_percentage_avg*stack[sequence[i]][sequence[j]][sequence[ip]][sequence[jp]]) + penalty_size;
                    }
                    else if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)
                    {
                        int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
                        if (branch1 == 1)
                        {
                            i2 = sequence[i];
                            j2 = sequence[j];
                            k2 = sequence[i+1];
                            ip2 = sequence[ip];
                            jp2 = sequence[jp];
                        }
                        else        // it's upside down
                        {
                            i2 = sequence[jp];
                            j2 = sequence[ip];
                            k2 = sequence[j-1];
                            ip2 = sequence[j];
                            jp2 = sequence[i];
                        }
                        energy = (PARAMTYPE) round(pmo_percentage_avg*bulge1_pmo[i2][j2][k2][ip2][jp2] + rna_percentage_avg*bulge1[i2][j2][k2][ip2][jp2]);
                    }
                }
                else
                {
                    // bulge of size bigger than 1
                    // check if (i,j) and (ip,jp) can pair
                    //printf ("REAL len = %d, pen = %Lg\n", l, penalty_size);
                    energy = penalty_size + (PARAMTYPE) round(pmo_percentage*AU_penalty_pmo(sequence[i],sequence[j]) + rna_percentage*AU_penalty(sequence[i],sequence[j])) + (PARAMTYPE) round(pmo_percentage_p*AU_penalty_pmo(sequence[ip], sequence[jp]) + rna_percentage_p*AU_penalty(sequence[ip], sequence[jp]));
                }
            }
            // it is an internal loop (not a bulge)
            else
            {
                l = branch1+branch2;
                penalty_size = (PARAMTYPE) round(pmo_percentage_avg*penalty_by_size_pmo(l, 'I') + rna_percentage_avg*penalty_by_size(l, 'I'));
                asym_penalty = (PARAMTYPE) round(pmo_percentage_avg*asymmetry_penalty_pmo(branch1, branch2) + rna_percentage_avg*asymmetry_penalty(branch1, branch2));
                //if (l == 5  && abs(branch1-branch2) == 3)
                //    printf ("REAL pen = %Lg, asym = %Lf\n", penalty_size, asym_penalty);

                
                if ((branch1 == 1 || branch2 == 1) && (misc.gail_rule || misc_pmo.gail_rule))
                // If gail_rule is set to 1 in miscloop file,
                // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                {
                    i_j_energy = 0;
                    ip_jp_energy = 0;

//#if (MODEL == SIMPLE)
                    // In the simple model I only use 3 parameters for tstacki, 
                    //  So tstacki[i][j][[0][0] never comes up, I just ignore it                    
                    //i_j_energy  =  tstacki[sequence[i]][sequence[j]][0][0];
                    //ip_jp_energy = tstacki[sequence[jp]][sequence[ip]][0][0];
//#elif (MODEL == EXTENDED)                    
                    if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                          ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                    {
                        i_j_energy = (PARAMTYPE) round(pmo_percentage*misc_pmo.internal_AU_closure + rna_percentage*misc.internal_AU_closure);
                    }
                    if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                          ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                    {
                        ip_jp_energy = (PARAMTYPE) round(pmo_percentage_p*misc_pmo.internal_AU_closure + rna_percentage_p*misc.internal_AU_closure);
                    }    
//#endif                                    
                    //printf ("IN GET_ENERGY, energy is %lf\n", i_j_energy + ip_jp_energy);
                }
                else
                {
                    i_j_energy   = (PARAMTYPE) round(pmo_percentage*tstacki_pmo[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]] + rna_percentage*tstacki[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]]);
                    ip_jp_energy = (PARAMTYPE) round(pmo_percentage_p*tstacki_pmo[sequence[jp]][sequence[ip]][sequence[jp+1]][sequence[ip-1]] + rna_percentage_p*tstacki[sequence[jp]][sequence[ip]][sequence[jp+1]][sequence[ip-1]]);
                    i_j_energy += (PARAMTYPE) round(pmo_percentage_avg*special_energy_internal_pmo(sequence, i,j,ip,jp) + rna_percentage_avg*special_energy_internal(sequence, i,j,ip,jp));
                }
                energy = i_j_energy + ip_jp_energy + penalty_size + asym_penalty;
            }
        }
    }
    //printf ("REAL:   i=%d, j=%d, ip=%d, jp=%d, energy=%g\n", i, j, ip, jp, energy);
    return energy;
}

// NOT UPDATED FOR MODEL=EXTENDED
PARAMTYPE s_internal_loop::get_enthalpy (int i, int j, int ip, int jp, int *sequence)
// returns the enthalpy of the internal loop closed at (i,j,ip,jp)
{
    PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
    int branch1, branch2, l;

    branch1 = ip-i-1;
    branch2 = j-jp-1;

    if (branch1 != 0 || branch2 != 0)
    {
        // check if it is a bulge loop of size 1
        // check if it is int11 or int21 or int22
        if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
        {
                        // int11[i][j][i+1][j-1][ip][jp]
            energy = enthalpy_int11 [sequence[i]][sequence[j]]
                           [sequence[i+1]][sequence[j-1]]
                           [sequence[ip]][sequence[jp]];
        }
        else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
        {
            // int21[i][j][i+1][j-1][ip][jp][jp+1]
            energy = enthalpy_int21 [sequence[i]][sequence[j]]
                           [sequence[i+1]][sequence[j-1]]
                           [sequence[ip]][sequence[jp]]
                           [sequence[jp+1]];
        }
        else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
        {
            // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
            energy = enthalpy_int21 [sequence[jp]][sequence[ip]]
                           [sequence[j-1]][sequence[ip-1]]
                           [sequence[j]][sequence[i]]
                           [sequence[i+1]];
        }
        else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
        {
            // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
            energy = enthalpy_int22 [sequence[i]][sequence[j]]
                           [sequence[i+1]][sequence[j-1]]
                           [sequence[ip]][sequence[jp]]
                           [sequence[ip-1]][sequence[jp+1]];
        }
        else
        {
            // this case is not int11, int21, int22

            // check if it is a bulge
            if (branch1 == 0 || branch2 == 0)
            {
                l = branch1+branch2;
                penalty_size = penalty_by_size_enthalpy (l, 'B');
                if (l == 1)
                {
                    // bulge of size 1
                    // stack[i][j][i+1][j-1]
                    energy = enthalpy_stack [sequence[i]][sequence[j]]
                                   [sequence[ip]][sequence[jp]] +
                               penalty_size;
                }
                else
                {
                    // bulge of size bigger than 1
                    // check if (i,j) and (ip,jp) can pair
                    energy = penalty_size +
                               AU_penalty_enthalpy (sequence[i],sequence[j]) +
                               AU_penalty_enthalpy (sequence[ip], sequence[jp]);
                }
            }
            // it is an internal loop (not a bulge)
            else
            {
                l = branch1+branch2;
                penalty_size = penalty_by_size_enthalpy (l, 'I');
                asym_penalty = asymmetry_penalty_enthalpy_pmo (branch1, branch2);

                if ((branch1 == 1 || branch2 == 1) && enthalpy_misc.gail_rule)
                // If gail_rule is set to 1 in miscloop file,
                // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                {
                    i_j_energy  =  enthalpy_tstacki[sequence[i]][sequence[j]]
                                          [0][0];
                    ip_jp_energy = enthalpy_tstacki[sequence[jp]][sequence[ip]]
                                          [0][0];
                }
                else
                {
                    i_j_energy   = enthalpy_tstacki[sequence[i]][sequence[j]]
                                          [sequence[i+1]][sequence[j-1]];
                    ip_jp_energy = enthalpy_tstacki[sequence[jp]][sequence[ip]]
                                          [sequence[jp+1]][sequence[ip-1]];
                    
                    //i_j_energy += special_enthalpy_internal (i,j,ip,jp);
                }
                energy = i_j_energy + ip_jp_energy + penalty_size +
                                         asym_penalty;
            }
        }
    }
    return energy;
}


PARAMTYPE count_int21_MODEL_EXTENDED (double *counter, int ii, int jj, int kk, int ll, int mm, int nn, int oo)
// do the counts for int21, when MODEL is EXTENDED
// return the energy, for checking purposes
{
    PARAMTYPE energy = 0;
    if (parsi_int21 == T99)     return 0;
    char type[100];
    int index;

    // SHOULD NOT apply rule 1 if !parsi_int21, they are independent parameters
//     if (parsi_int21)
//     {
//         apply_rule_1 (kk, ll, kk, ll);
//         apply_rule_1 (kk, oo, kk, oo);
//     }


    // The 6 parameters are counted only if it's parsi_int21 or it's !parsi_int21 and is part of experimental_addition
    if (parsi_int21 == PARSI || creating_model ||
        (parsi_int21 == LAVISH && int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] < INF))
    {
        // exactly the same as above
        // try to follow the model proposed by Badhwar_Znosko_2007
        // the initiation appears in all cases
        index = structure_type_index ("misc.internal21_initiation");
        counter[index]++;
        energy += misc.internal21_initiation;
        //printf ("Adding misc.internal21_initiation = %g\n", misc.internal21_initiation);
        // look for AU closure
        if ((ii==A && jj==U) || (ii==U && jj==A))
        {
            index = structure_type_index ("misc.internal21_AU_closure");
            counter[index]++;
            energy += misc.internal21_AU_closure;
            //printf ("Adding misc.internal21_AU_closure = %g\n", misc.internal21_AU_closure);
        }
        if ((mm==A && nn==U) || (mm==U && nn==A))
        {
            index = structure_type_index ("misc.internal21_AU_closure");
            counter[index]++;
            energy += misc.internal21_AU_closure;
            //printf ("Adding misc.internal21_AU_closure = %g\n", misc.internal21_AU_closure);
        }
        // look for GU closure
        if ((ii==G && jj==U) || (ii==U && jj==G))
        {
            index = structure_type_index ("misc.internal21_GU_closure");
            counter[index]++;
            energy += misc.internal21_GU_closure;
            //printf ("Adding misc.internal21_GU_closure = %g\n", misc.internal21_GU_closure);
        }
        if ((mm==G && nn==U) || (mm==U && nn==G))
        {
            index = structure_type_index ("misc.internal21_GU_closure");
            counter[index]++;
            energy += misc.internal21_GU_closure;
            //printf ("Adding misc.internal21_GU_closure = %g\n", misc.internal21_GU_closure);
        }
        // look for AG mismatch - but not applied to 5'RA/3'YG loops
        if ((kk==A && ll==G &&   ii!=A && ii!=G   &&   jj!=U && jj!=C) ||
            (kk==G && ll==A) ||
            (kk==G && oo==A &&   mm!=U && mm!=C   &&   nn!=A && nn!=G) ||
            (kk==A && oo==G))
        {
            index = structure_type_index ("misc.internal21_AG_mismatch");
            counter[index]++;
            energy += misc.internal21_AG_mismatch;
            //printf ("Adding misc.internal21_AG_mismatch = %g\n", misc.internal21_AG_mismatch);
        }        
        // look for GG mismatch
        if (kk==G && (ll==G || oo==G))
        {
            index = structure_type_index ("misc.internal21_GG_mismatch");
            counter[index]++;
            energy += misc.internal21_GG_mismatch;
            //printf ("Adding misc.internal21_GG_mismatch = %g\n", misc.internal21_GG_mismatch);
        }
        // look for UU mismatch
        if (kk==U && (ll==U || oo==U))
        {
            index = structure_type_index ("misc.internal21_UU_mismatch");
            counter[index]++;
            energy += misc.internal21_UU_mismatch;
            //printf ("Adding misc.internal21_UU_mismatch = %g\n", misc.internal21_UU_mismatch);
        }
        // IN THIS CASE, int21 is just the additional value, on top of the above values, in order to make them be what the experiments say, and not an approximation
        // HMM - MAYBE THIS IS NOT SUCH A GOOD IDEA ACTUALLY, BUT THAT'S WHAT'S RECOMMENDED BY THE OPTICAL MELTING PAPERS
        // TODO: separate into experimental_addition and not
    
        if (parsi_int21 == LAVISH)
        {
            //sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
            if (creating_model)
                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
            else
                sprintf (type, "int21_experimental_addition[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
            index = structure_type_index (type);
            counter[index]++;
            //energy +=  int21[ii][jj][kk][ll][mm][nn][oo];
            
            if (creating_model)
                // Note I'm not ADDING to energy here, but writing everything in it
                energy =  int21[ii][jj][kk][ll][mm][nn][oo];
            else
                energy +=  int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo];
                        
        }
    }
    else    // it is parsi_int21, and  int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] is infinity
    {
        sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
        index = structure_type_index (type);
        counter[index]++;
        energy +=  int21[ii][jj][kk][ll][mm][nn][oo];
        
    }
    return energy;
}


PARAMTYPE count_int11_MODEL_EXTENDED (double *counter, int ii, int jj, int kk, int ll, int mm, int nn)
// do the counts for int11, when MODEL is EXTENDED
{
    PARAMTYPE energy = 0;
    if (parsi_int11 == T99)     return 0;
    char type[100];
    int index;

    // SHOULD NOT apply rule 1 if !parsi_int11, they are independent parameters
//     if (parsi_int11)    // I don't think it matters if I add rule 1 or not, but add it for consistency
//     {        
//         // Apply rule 1
//         apply_rule_1 (kk, ll, kk, ll);
//         // Applying rule 1 might not get the order of ii, jj, kk, ll, mm, nn to be in the first symmetric part, which is part of the feature set
//     }

    // The 10 parameters are counted only if it's parsi_int11 or it's !parsi_int11 and is part of experimental_addition
    if (parsi_int11 == PARSI || creating_model ||
        ((parsi_int11 == LAVISH || parsi_int11 == HLI) && int11_experimental_addition[ii][jj][kk][ll][mm][nn] < INF))
    {    
        //printf ("In COUNT, int11_expadd(%d,%d,%d,%d,%d,%d) = %Lg\n", ii,jj,kk,ll,mm,nn, int11_experimental_addition[ii][jj][kk][ll][mm][nn]);
        //printf ("parsi=%d, creating_model=%d\n", parsi_int11, creating_model);
        // try to follow the model proposed by Davis_Znosko_2007
        // look for AU closure
        if ((ii==A && jj==U) || (ii==U && jj==A))
        {
            index = structure_type_index ("misc.internal11_AU_closure");
            counter[index]++;
            energy += misc.internal11_AU_closure;
            //printf ("In COUNT11 add misc.internal11_AU_closure = %lf\n", misc.internal11_AU_closure);        
        }
        if ((mm==A && nn==U) || (mm==U && nn==A))
        {
            index = structure_type_index ("misc.internal11_AU_closure");
            counter[index]++;
            energy += misc.internal11_AU_closure;
            //printf ("In COUNT11 add misc.internal11_AU_closure = %lf\n", misc.internal11_AU_closure);
        }
        // look for GU closure
        if ((ii==G && jj==U) || (ii==U && jj==G))
        {
            index = structure_type_index ("misc.internal11_GU_closure");
            counter[index]++;
            energy += misc.internal11_GU_closure;
            //printf ("In COUNT11 add misc.internal11_GU_closure = %lf\n", misc.internal11_GU_closure);
        }
        if ((mm==G && nn==U) || (mm==U && nn==G))
        {
            index = structure_type_index ("misc.internal11_GU_closure");
            counter[index]++;
            energy += misc.internal11_GU_closure;
            //printf ("In COUNT11 add misc.internal11_GU_closure = %lf\n", misc.internal11_GU_closure);
        }
        // look for AG mismatch
        if ((kk==A && ll==G) || (kk==G && ll==A))
        {
            index = structure_type_index ("misc.internal11_AG_mismatch");
            counter[index]++;
            energy += misc.internal11_AG_mismatch;
            //printf ("In COUNT11 add misc.internal11_AG_mismatch = %lf\n", misc.internal11_AG_mismatch);
        }    
        // look for GG mismatch
        if (kk==G && ll==G)
        {
            index = structure_type_index ("misc.internal11_GG_mismatch");
            counter[index]++;
            energy += misc.internal11_GG_mismatch;
            //printf ("In COUNT11 add misc.internal11_GG_mismatch = %lf\n", misc.internal11_GG_mismatch);
        }
        // look for UU mismatch
        if (kk==U && ll==U)
        {
            index = structure_type_index ("misc.internal11_UU_mismatch");
            counter[index]++;
            energy += misc.internal11_UU_mismatch;
            //printf ("In COUNT11 add misc.internal11_UU_mismatch = %lf\n", misc.internal11_UU_mismatch);
        }
        // check if it is internal11_5YRR_5YRR
    
        if (isY(ii) && isR(jj) && isR(kk) && isR(ll) && isR(mm) && isY(nn))
        {
            index = structure_type_index ("misc.internal11_5YRR_5YRR");
            counter[index]++;
            energy += misc.internal11_5YRR_5YRR;
            //printf ("In COUNT11 add misc.internal11_5YRR_5YRR = %lf\n", misc.internal11_5YRR_5YRR);
        }    
        if ( isR(ii) && isY(jj) && isY(kk) && isY(ll) && isY(mm) && isR(nn) )
        {
            index = structure_type_index ("misc.internal11_5RYY_5RYY");
            counter[index]++;
            energy += misc.internal11_5RYY_5RYY;
            //printf ("In COUNT11 add misc.internal11_5RYY_5RYY = %lf\n", misc.internal11_5RYY_5RYY);
        }    
        if ( isY(ii) && isR(jj) && isY(kk) && isY(ll) && isR(mm) && isY(nn) )
        {
            index = structure_type_index ("misc.internal11_5YYR_5YYR");
            counter[index]++;
            energy += misc.internal11_5YYR_5YYR;
            //printf ("In COUNT11 add misc.internal11_5YYR_5YYR = %lf\n", misc.internal11_5YYR_5YYR);
        }    
        if ( (isY(ii) && isR(jj) && isR(kk) && isY(ll) && isY(mm) && isR(nn)) ||
            (isR(ii) && isY(jj) && isY(kk) && isR(ll) && isR(mm) && isY(nn)) )    
        {
            index = structure_type_index ("misc.internal11_5YRY_5RYR");
            counter[index]++;
            energy += misc.internal11_5YRY_5RYR;
            //printf ("In COUNT11 add misc.internal11_5YRY_5RYR = %lf\n", misc.internal11_5YRY_5RYR);
            //printf ("ii=%d. jj=%d, kk=%d, ll=%d, mm=%d, nn=%d\n", ii, jj, kk, ll, mm, nn);
        }
        if ( (isR(ii) && isY(jj) && isR(kk) && isY(ll) && isY(mm) && isR(nn)) ||
            (isR(ii) && isY(jj) && isY(kk) && isR(ll) && isY(mm) && isR(nn)) )    
        {
            index = structure_type_index ("misc.internal11_5RRY_5RYY");
            counter[index]++;
            energy += misc.internal11_5RRY_5RYY;
            //printf ("In COUNT11 add misc.internal11_5RRY_5RYY = %lf\n", misc.internal11_5RRY_5RYY);
        }    
    
        // IN THIS CASE, int21 is just the additional value, on top of the above values, in order to make them be what the experiments say, and not an approximation
        // HMM - MAYBE THIS IS NOT SUCH A GOOD IDEA ACTUALLY, BUT THAT'S WHAT'S RECOMMENDED BY THE OPTICAL MELTING PAPERS
    
        // the following is counted only if it's part of experimental_addition, otherwise it's not
        //if (int11_experimental_addition[ii][jj][kk][ll][mm][nn] < INF)
        if (parsi_int11 == LAVISH || parsi_int11 == HLI) 
        {
            // NO RULE 1
            // If creating_model, we don't have experimental additions yet in the string_params, but the index is the same            
            if (ii*100000 + jj*10000 + kk*1000 + ll*100 + mm*10 + nn <= nn*100000 + mm*10000 + ll*1000+ kk*100 + jj*10 + ii)
            {
                if (creating_model)
                    sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
                else
                    sprintf (type, "int11_experimental_addition[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
            }
            else
            {
                if (creating_model)
                    sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", nn, mm, ll, kk, jj, ii);
                else
                    sprintf (type, "int11_experimental_addition[%d][%d][%d][%d][%d][%d]", nn, mm, ll, kk, jj, ii);
            }
            //printf ("IN COUNT int11, looking for %s\n", type);
            //sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
            index = structure_type_index (type);
            counter[index]++;
            //energy +=  int11[ii][jj][kk][ll][mm][nn];
            // the following shouldn't be affected by rule 1
            if (creating_model)
                // Note I'm not ADDING to energy here, but writing everything in it
                energy =  int11[ii][jj][kk][ll][mm][nn];
            else
                energy +=  int11_experimental_addition[ii][jj][kk][ll][mm][nn];
        }
    }
    else    // !parsi_int11 and int11_experimental_addition = INF
    {
        // Applying rule 1 might not get the order of ii, jj, kk, ll, mm, nn to be in the first symmetric part, which is part of the feature set
        if (ii*100000 + jj*10000 + kk*1000 + ll*100 + mm*10 + nn <= nn*100000 + mm*10000 + ll*1000+ kk*100 + jj*10 + ii)
            sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
        else
            sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", nn, mm, ll, kk, jj, ii);
        //printf ("IN COUNT int11, looking for %s\n", type);
        index = structure_type_index (type);
        counter[index]++;
        // the following shouldn't be affected by rule 1 
        energy +=  int11[ii][jj][kk][ll][mm][nn];
    }
    return energy;
}



PARAMTYPE count_int22_MODEL_EXTENDED (double *counter, int ii, int jj, int kk, int ll, int mm, int nn, int oo, int pp)
// do the counts for int21, when MODEL is EXTENDED
{
    PARAMTYPE energy = 0;
    if (parsi_int22 == T99)     return 0;
    char type[100];
    int index;

    // SHOULD NOT apply rule 1 if !parsi_int11, they are independent parameters
//     if (parsi_int22)
//     {
//         // Apply rule 1 twice
//         apply_rule_1 (kk, ll, kk, ll);
//         apply_rule_1 (oo, pp, oo, pp);
//         // Applying rule 1 might not get the order of ii, jj, kk, ll, mm, nn, oo, pp to be in the first symmetric part, which is part of the feature set
//     }

    // The 6 parameters are counted only if it's parsi_int22 or it's !parsi_int22 and is part of experimental_addition
    if (parsi_int22 == PARSI || creating_model || 
        (parsi_int22 == LAVISH && int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] < INF))
    {
        //printf ("INT22_EXP_ADD < INF: int22_experimental_addition[%d][%d][%d][%d][%d][%d][%d][%d] = %lf\n", ii, jj, kk, ll, mm, nn, oo, pp, int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp]);
        // try to follow the model proposed by Christiansen_Znosko_2008
        if (is_int22_group_1 (kk, ll, oo, pp))
        {
            index = structure_type_index ("misc.internal22mid_group1");
            counter[index]++;
            energy += misc.internal22mid_group1;
            //printf ("Add misc.internal22mid_group1 = %lf\n", misc.internal22mid_group1);
        }
        else if (is_int22_group_2 (kk, ll, oo, pp))
        {
            index = structure_type_index ("misc.internal22mid_group2");
            counter[index]++;
            energy += misc.internal22mid_group2;
            //printf ("Add misc.internal22mid_group2 = %lf\n", misc.internal22mid_group2);
        }
        else if (is_int22_group_3 (kk, ll, oo, pp))
        {
            is_int22_group_3 (kk, ll, oo, pp);
            index = structure_type_index ("misc.internal22mid_group3");
            counter[index]++;
            energy += misc.internal22mid_group3;
            //printf ("Add misc.internal22mid_group3 = %lf\n", misc.internal22mid_group3);
        }
        else if (is_int22_group_4 (kk, ll, oo, pp))
        {
            index = structure_type_index ("misc.internal22mid_group4");
            counter[index]++;
            energy += misc.internal22mid_group4;
            //printf ("Add misc.internal22mid_group4 = %lf\n", misc.internal22mid_group4);
        }
        if ((ii == A && jj == U) || (ii == U && jj == A))
        {
            index = structure_type_index ("misc.internal22_AU_closure");
            counter[index]++;
            energy += misc.internal22_AU_closure;
            //printf ("Add misc.internal22_AU_closure = %lf\n", misc.internal22_AU_closure);
        }
        else if ((ii == G && jj == U) || (ii == U && jj == G))
        {
            index = structure_type_index ("misc.internal22_GU_closure");
            counter[index]++;
            energy += misc.internal22_GU_closure;
            //printf ("Add misc.internal22_GU_closure = %lf\n", misc.internal22_GU_closure);
        }
        if ((mm == A && nn == U) || (mm == U && nn == A))
        {
            index = structure_type_index ("misc.internal22_AU_closure");
            counter[index]++;
            energy += misc.internal22_AU_closure;
            //printf ("Add misc.internal22_AU_closure = %lf\n", misc.internal22_AU_closure);
        }
        else if ((mm == G && nn == U) || (mm == U && nn == G))
        {
            index = structure_type_index ("misc.internal22_GU_closure");
            counter[index]++;
            energy += misc.internal22_GU_closure;
            //printf ("Add misc.internal22_GU_closure = %lf\n", misc.internal22_GU_closure);
        }
    
        // the following is counted only if it's part of experimental_addition, otherwise it's not
        if (parsi_int22 == LAVISH)
        {
            // Applying rule 1 might not get the order of ii, jj, kk, ll, mm, nn to be in the first symmetric part, which is part of the feature set
            if (ii*10000000 + jj*1000000 + kk*100000 + ll*10000 + mm*1000 + nn*100 + oo*10 + pp <=
                nn*10000000 + mm*1000000 + pp*100000 + oo*10000 + jj*1000 + ii*100 + ll*10 + kk)
            {
                if (creating_model)
                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
                else
                    sprintf (type, "int22_experimental_addition[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
            }
            else
            {
                if (creating_model)
                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nn, mm, pp, oo, jj, ii, ll, kk);
                else
                    sprintf (type, "int22_experimental_addition[%d][%d][%d][%d][%d][%d][%d][%d]", nn, mm, pp, oo, jj, ii, ll, kk);
            }
            index = structure_type_index (type);
            counter[index]++;
            // the following shouldn't be affected by rule 1
            if (creating_model)
                energy = int22[ii][jj][kk][ll][mm][nn][oo][pp];
            else
                energy += int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp];
        }
    }
    else    // !parsi_int11 and int11_experimental_addition = INF
    {
        // Applying rule 1 might not get the order of ii, jj, kk, ll, mm, nn to be in the first symmetric part, which is part of the feature set
        if (ii*10000000 + jj*1000000 + kk*100000 + ll*10000 + mm*1000 + nn*100 + oo*10 + pp <=
            nn*10000000 + mm*1000000 + pp*100000 + oo*10000 + jj*1000 + ii*100 + ll*10 + kk)
            sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
        else
            sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nn, mm, pp, oo, jj, ii, ll, kk);
        //printf ("IN COUNT int11, looking for %s\n", type);
        index = structure_type_index (type);
        counter[index]++;
        // the following shouldn't be affected by rule 1 
        energy +=  int22[ii][jj][kk][ll][mm][nn][oo][pp];
    }
    return energy;
}


void s_internal_loop::count_get_energy (int i, int j, int ip, int jp, int *sequence, double *counter)
// this function is needed for parameter learning, not for folding
// fill the counter vectro accordingly
// Mirela: Nov 23, 2003
{
    PARAMTYPE energy;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
    int branch1, branch2, l;
    char type[100];
    int index;
    branch1 = ip-i-1;
    branch2 = j-jp-1;    

    energy = 0;

    if (branch1 != 0 || branch2 != 0)
    {
        // check if it is a bulge loop of size 1
        // check if it is int11 or int21 or int22
        if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
        {
                        // int11[i][j][i+1][j-1][ip][jp]
            //energy = IGINF(int11 [sequence[i]][sequence[j]]
            //               [sequence[i+1]][sequence[j-1]]
            //               [sequence[ip]][sequence[jp]]);
            
            int ii, jj, kk, ll, mm, nn;
            if (sequence[i]*100000 + sequence[j]*10000 + sequence[i+1]*1000 + sequence[j-1]*100 + sequence[ip]*10 + sequence[jp] >
                sequence[jp]*100000 + sequence[ip]*10000 + sequence[j-1]*1000+ sequence[i+1]*100 + sequence[j]*10 + sequence[i])
            {                  
                ii = sequence[jp];
                jj = sequence[ip];
                kk = sequence[j-1];
                ll = sequence[i+1];
                mm = sequence[j];
                nn = sequence[i];
            }
            else
            {
                ii = sequence[i];
                jj = sequence[j];
                kk = sequence[i+1];
                ll = sequence[j-1];
                mm = sequence[ip];
                nn = sequence[jp];
            }
                  
            if (parsi_int11 == T99)
            {
                if ( ((ii==C && jj==G) || (ii==G && jj==C)) && ((mm==C && nn==G) || (mm==G && nn==C)))
                {
                    if (!can_pair(kk,ll))
                    {
                        sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
                        energy += int11[ii][jj][kk][ll][mm][nn];
                    }
                    else
                    {
                        sprintf (type, "misc.internal11_basic_mismatch");
                        energy += misc.internal11_basic_mismatch;
                    }
                    index = structure_type_index (type);
                    counter[index]++;
                }
                else if (watson_crick(ii,jj) && watson_crick(mm,nn) && kk==U && ll==U)
                {
                    sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
                    index = structure_type_index (type);
                    counter[index]++;
                    energy += int11[ii][jj][kk][ll][mm][nn];
                }
                else
                {
                    if (kk==G && ll==G)
                    {
                        sprintf (type, "misc.internal11_GG_mismatch");
                        energy += misc.internal11_GG_mismatch;
                    }
                    else
                    {
                        sprintf (type, "misc.internal11_basic_mismatch");
                        energy += misc.internal11_basic_mismatch;
                    }
                    index = structure_type_index (type);
                    counter[index]++;
    
                    if (has_AU_penalty(ii,jj))
                    {
                        sprintf (type, "misc.internal_AU_closure");
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += misc.internal_AU_closure;
                    }
                    if (has_AU_penalty(mm,nn))
                    {
                        sprintf (type, "misc.internal_AU_closure");
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += misc.internal_AU_closure;
                    }
                }
            }
            else
            {
                energy += count_int11_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn);
            }
        }
        else if (((branch1 == 1 && branch2 == 2) || (branch1 == 2 && branch2 == 1)) && !simple_internal_energy)
        {
            // int21[i][j][i+1][j-1][ip][jp][jp+1]
//             energy = IGINF(int21 [sequence[i]][sequence[j]]
//                            [sequence[i+1]][sequence[j-1]]
//                            [sequence[ip]][sequence[jp]]
//                            [sequence[jp+1]]);
                           
            int ii, jj, kk, ll, mm, nn, oo;
            if (branch1 == 1 && branch2 == 2)
            {
                ii = sequence[i];
                jj = sequence[j];
                kk = sequence[i+1];
                ll = sequence[j-1];
                mm = sequence[ip];
                nn = sequence[jp];
                oo = sequence[jp+1];
            }
            else        //  branch1 == 2 && branch2 == 1
            {
                ii = sequence[jp];
                jj = sequence[ip];
                kk = sequence[j-1];
                ll = sequence[ip-1];
                mm = sequence[j];
                nn = sequence[i];
                oo = sequence[i+1];
            }

            if (parsi_int21 == T99)
            {
                if ((ii==C && jj==G && mm==C && nn==G) ||  // these are already filled above, except what can pair inside
                    (ii==G && jj==C && mm==G && nn==C))
                {
                    if (can_pair(kk,ll) || can_pair(kk,oo))
                    {
                        sprintf (type, "misc.internal21_match");
                        energy += misc.internal21_match;
                    }
                    else
                    {
                        sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
                        energy += int21[ii][jj][kk][ll][mm][nn][oo];
                    }
                    index = structure_type_index (type);
                    counter[index]++;
                }
                else
                {
                    if (can_pair(kk,ll) || can_pair(kk,oo))
                    {
                        sprintf (type, "misc.internal21_match");
                        index = structure_type_index (type);
                        counter[index] ++;
                        energy += misc.internal21_match;
                    }
                    else
                    {
                        sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kk, ll, C, G, oo);
                        index = structure_type_index (type);
                        counter[index] += 0.5;
                        energy += (PARAMTYPE) round(0.5 * int21[C][G][kk][ll][C][G][oo]);
                        
                        sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kk, ll, G, C, oo);
                        index = structure_type_index (type);
                        counter[index] += 0.5;
                        energy += (PARAMTYPE) round(0.5 * int21[G][C][kk][ll][G][C][oo]);
                    }        
                    if (has_AU_penalty(ii,jj))
                    {
                        sprintf (type, "misc.internal21_AU_closure");
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += misc.internal21_AU_closure;
                    }    
                    if (has_AU_penalty(mm,nn))    
                    {
                        sprintf (type, "misc.internal21_AU_closure");
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += misc.internal21_AU_closure;                      
                    }    
                }
            }
            else
            {
                energy += count_int21_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn, oo);
            }
        }
        else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
        {
            // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
//             energy = IGINF(int22 [sequence[i]][sequence[j]]
//                            [sequence[i+1]][sequence[j-1]]
//                            [sequence[ip]][sequence[jp]]
//                            [sequence[ip-1]][sequence[jp+1]]);
            int ii, jj, kk, ll, mm, nn, oo, pp;
            if (sequence[i]*10000000 + sequence[j]*1000000 + sequence[i+1]*100000 + sequence[j-1]*10000 +
                sequence[ip]*1000 + sequence[jp]*100 + sequence[ip-1]*10 + sequence[jp+1] >
                sequence[jp]*10000000 + sequence[ip]*1000000 + sequence[jp+1]*100000 + sequence[ip-1]*10000 +
                sequence[j]*1000 + sequence[i]*100 + sequence[j-1]*10 + sequence[i+1])
            {                
                ii = sequence[jp];
                jj = sequence[ip];
                kk = sequence[jp+1];
                ll = sequence[ip-1];
                mm = sequence[j];
                nn = sequence[i];
                oo = sequence[j-1];            
                pp = sequence[i+1];                
            }
            else             // the symmetric case
            {
                ii = sequence[i];
                jj = sequence[j];
                kk = sequence[i+1];
                ll = sequence[j-1];
                mm = sequence[ip];
                nn = sequence[jp];
                oo = sequence[ip-1];            
                pp = sequence[jp+1];
            }                
             
            if (parsi_int22 == T99)
            {
                if (nn==ii && mm==jj && pp==kk && oo==ll & watson_crick(ii,jj) && !watson_crick(kk,ll))
                {
                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
                    index = structure_type_index (type);
                    counter[index]++;
                    energy += int22[ii][jj][kk][ll][mm][nn][oo][pp];
                }
    
                int ii2, jj2, mm2, nn2;
                if (ii==G && jj==U)   ii2 = A;     else ii2 = ii;
                if (ii==U && jj==G)   jj2 = A;     else jj2 = jj;
                if (mm==G && nn==U)   mm2 = A;     else mm2 = mm;
                if (mm==U && nn==G)   nn2 = A;     else nn2 = nn;
    
                if (watson_crick(kk,ll) || watson_crick(oo,pp))
                {
                    sprintf (type, "misc.internal22_match");
                    index = structure_type_index (type);
                    counter[index]++;
                    energy += misc.internal22_match;
                }
                else if ( ((ii==G && jj==U) || (ii==U && jj==G) || (mm==G && nn==U) || (mm==U && nn==G)) &&
                            (nn2==ii2 && mm2==jj2 && pp==kk && oo==ll))  // the UG closing pairs are the same as UA
                {
                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii2, jj2, kk, ll, mm2, nn2, oo, pp);
                    index = structure_type_index (type);
                    counter[index]++;
                    energy += int22[ii2][jj2][kk][ll][mm2][nn2][oo][pp];
                }
                else if (!(nn==ii && mm==jj && pp==kk && oo==ll))   // was already filled above
                {
                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii2, jj2, kk, ll, jj2, ii2, ll, kk);
                    index = structure_type_index (type);
                    counter[index] += 0.5;
                    energy += (PARAMTYPE) round(0.5 * int22[ii2][jj2][kk][ll][jj2][ii2][ll][kk]);
                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nn2, mm2, pp, oo, mm2, nn2, oo, pp);
                    index = structure_type_index (type);
                    counter[index] += 0.5;
                    energy += (PARAMTYPE) round(0.5 * int22[nn2][mm2][pp][oo][mm2][nn2][oo][pp]);
    
                    int result = check_stability_and_size (kk, ll, oo, pp);
                    switch (result)
                    {
                        case 1:
                            sprintf (type, "misc.internal22_delta_same_size");
                            energy += misc.internal22_delta_same_size;
                            break;
                        case 2:
                            sprintf (type, "misc.internal22_delta_different_size");
                            energy += misc.internal22_delta_different_size;
                            break;
                        case 3:
                            sprintf (type, "misc.internal22_delta_1stable_1unstable");
                            energy += misc.internal22_delta_1stable_1unstable;
                            break;
                        case 4:
                            sprintf (type, "misc.internal22_delta_AC");
                            energy += misc.internal22_delta_AC;
                            break;
                        default:
                            printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kk, ll, oo, pp);
                            exit(1);
                    }
                    index = structure_type_index (type);
                    counter[index]++;
                }
            }
            else
            {
                energy += count_int22_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn, oo, pp);
            }
        }       // end if 2x2
        else
        {
            // this case is not int11, int21, int22                                                                   
            // check if it is a bulge
            if (branch1 == 0 || branch2 == 0)
            {
                l = branch1+branch2;                            
                if (l == 1)
                {
                    if (parsi_bulge1 == T99)
                    {
                        // bulge of size 1
                        // stack[i][j][i+1][j-1]
                        energy += penalty_by_size (l, 'B');
                        count_penalty_by_size (l, 'B', counter);
                        if (1000*sequence[i] + 100*sequence[j] + 10*sequence[ip] + sequence[jp] >
                            1000*sequence[jp] + 100*sequence[ip] + 10*sequence[j] + sequence[i])
                        {
                            sprintf (type, "stack[%d][%d][%d][%d]", sequence[jp], sequence[ip], sequence[j], sequence[i]);
                            energy += stack[sequence[jp]][sequence[ip]][sequence[j]][sequence[i]];
                        }
                        else
                        {
                            sprintf (type, "stack[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[ip], sequence[jp]);
                            energy += stack[sequence[i]][sequence[j]][sequence[ip]][sequence[jp]];
                        }
                        index = structure_type_index (type);
                        counter[index]++;
                    }
                    else
                    {
                        int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
                        if (branch1 == 1)
                        {
                            i2 = sequence[i];
                            j2 = sequence[j];
                            k2 = sequence[i+1];
                            ip2 = sequence[ip];
                            jp2 = sequence[jp];
                        }
                        else        // it's upside down
                        {
                            i2 = sequence[jp];
                            j2 = sequence[ip];
                            k2 = sequence[j-1];
                            ip2 = sequence[j];
                            jp2 = sequence[i];
                        }
                        if (parsi_bulge1 == PARSI)
                        {
                            // energy is size + stack + bulge(k2)
                            // NO! I should NOT add size, it's already in bulge(k2)
                            //penalty_size = penalty_by_size (l, 'B');    
                            //count_penalty_by_size (l, 'B', counter);
                            //energy = IGINF(stack [sequence[i]][sequence[j]][sequence[ip]][sequence[jp]]) + penalty_size;
                            if (1000*sequence[i] + 100*sequence[j] + 10*sequence[ip] + sequence[jp] >
                                1000*sequence[jp] + 100*sequence[ip] + 10*sequence[j] + sequence[i])
                            {
                                sprintf (type, "stack[%d][%d][%d][%d]", sequence[jp], sequence[ip], sequence[j], sequence[i]);
                                energy += stack[sequence[jp]][sequence[ip]][sequence[j]][sequence[i]];
                            }
                            else
                            {
                                sprintf (type, "stack[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[ip], sequence[jp]);
                                energy += stack[sequence[i]][sequence[j]][sequence[ip]][sequence[jp]];
                            }
                            index = structure_type_index (type);
                            counter[index]++;
                            sprintf (type, "bulge%c", int_to_nuc(k2));
                            index = structure_type_index (type);
                            counter[index]++;
                            if (k2 == A)        energy += bulgeA;
                            else if (k2 == C)   energy += bulgeC;
                            else if (k2 == G)   energy += bulgeG;
                            else if (k2 == U)   energy += bulgeU;
                        }
                        else if (parsi_bulge1 == LAVISH)
                        {
                            sprintf (type, "bulge1[%d][%d][%d][%d][%d]", i2, j2, k2, ip2, jp2);
                            index = structure_type_index (type);
                            counter[index]++;
                            energy += bulge1[i2][j2][k2][ip2][jp2];
                        }
                    }
                }
                else
                {
                    // bulge of size bigger than 1
                    // check if (i,j) and (ip,jp) can pair
                    energy += penalty_by_size (l, 'B');
                    count_penalty_by_size (l, 'B', counter);
                    energy +=  AU_penalty (sequence[i],sequence[j]) +
                               AU_penalty (sequence[ip], sequence[jp]);
                    count_AU_penalty (sequence[i],sequence[j], counter);
                    count_AU_penalty (sequence[ip], sequence[jp], counter);
                }
            }
            // it is an internal loop (not a bulge)
            else
            {
                l = branch1+branch2;
                energy += penalty_by_size (l, 'I');
                count_penalty_by_size (l, 'I', counter);
                //printf ("In COUNT, add size %lf\n", penalty_by_size(l,'I'));
                energy += asymmetry_penalty (branch1, branch2);
                count_asymmetry_penalty (branch1, branch2, counter);
                //printf ("In COUNT, add asymmetry %lf\n",asymmetry_penalty (branch1, branch2));
                                                                                                                                                                                                                                                                                                                                                                                                                                    
                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                // If gail_rule is set to 1 in miscloop file,
                // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                {
                    
//#if (MODEL == SIMPLE)   
                    // In the simple model I only use 3 parameters for tstacki, 
                    //  So tstacki[i][j][[0][0] never comes up, I just ignore it
                    /*         
                    i_j_energy  =  IGINF(tstacki[sequence[i]][sequence[j]][0][0]);
                    sprintf (type, "tstacki[%d][%d][0][0]", sequence[i], sequence[j]);
                    index = structure_type_index (type);
                    counter[index]++;

                    ip_jp_energy = IGINF(tstacki[sequence[jp]][sequence[ip]][0][0]);
                    sprintf (type, "tstacki[%d][%d][0][0]", sequence[jp], sequence[ip]);
                    index = structure_type_index (type);
                    counter[index]++;
                    */
//#elif (MODEL == EXTENDED)                                        
                    // actually, just use 3 parameters instead of the tstacki table
                    if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                        ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                    {
                        //internal_AU_closure includes terminal_AU_penalty
                        sprintf (type, "misc.internal_AU_closure");
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += misc.internal_AU_closure;
                        //printf ("energy is %lf \n", misc.internal_AU_closure);
                    }
                    // actually, just use 3 parameters instead of the tstacki table
                    if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                        ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                    {
                        //internal_AU_closure includes terminal_AU_penalty
                        sprintf (type, "misc.internal_AU_closure");
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += misc.internal_AU_closure;
                        //printf ("energy is %lf\n", misc.internal_AU_closure);
                    }                    
//#endif
                }
                else
                {
                    energy += count_special_internal (counter, sequence, i, j, ip, jp);
                    
                    i_j_energy   = IGINF(tstacki[sequence[i]][sequence[j]]
                                          [sequence[i+1]][sequence[j-1]]);

                    if (parsi_tstacki == T99 || parsi_tstacki == PARSI)
                    {
                        // actually, just use 3 parameters instead of the tstacki table
                        if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                            ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                        {
                            //internal_AU_closure includes terminal_AU_penalty
                            //sprintf (type, "misc.terminal_AU_penalty");
                            //index = structure_type_index (type);
                            //counter[index]++;                        
                            sprintf (type, "misc.internal_AU_closure");
                            index = structure_type_index (type);
                            counter[index]++;
                            energy += misc.internal_AU_closure;
                            //printf ("In COUNT, add misc.internal_AU_closure + %lf\n", misc.internal_AU_closure);
                        }
                        if (parsi_tstacki == T99)
                        {
                            if ((sequence[i+1] == A && sequence[j-1] == G) || (sequence[i+1] == G && sequence[j-1] == A))
                            {
                                sprintf (type, "misc.internal_GA_AG_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_GA_AG_mismatch;                         
                            }
                        }
                        else if (parsi_tstacki == PARSI)
                        {
                            if (sequence[i+1] == A && sequence[j-1] == G)
                            {
                                sprintf (type, "misc.internal_AG_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_AG_mismatch;
                                //printf ("In COUNT, add misc.internal_AG_mismatch + %lf\n", misc.internal_AG_mismatch);
                            }
                            if (sequence[i+1] == G && sequence[j-1] == A)
                            {
                                sprintf (type, "misc.internal_GA_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_GA_mismatch;
                                //printf ("In COUNT, add misc.internal_GA_mismatch + %lf\n", misc.internal_GA_mismatch);
                            }
                            if (sequence[i+1] == G && sequence[j-1] == G)
                            {
                                sprintf (type, "misc.internal_GG_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_GG_mismatch;
                                //printf ("In COUNT, add misc.internal_GG_mismatch + %lf\n", misc.internal_GG_mismatch);
                            }
                        }
                        if (sequence[i+1] == U && sequence[j-1] == U)
                        {
                            sprintf (type, "misc.internal_UU_mismatch");
                            index = structure_type_index (type);
                            counter[index]++;
                            energy += misc.internal_UU_mismatch;
                            //printf ("In COUNT, add misc.internal_UU_mismatch + %lf\n", misc.internal_UU_mismatch);
                        }
                    }
                    else if (parsi_tstacki == LAVISH)
                    {
                        // SHOULD NOT apply rule 1 here, they are separate parameters
                        sprintf (type, "tstacki[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[i+1], sequence[j-1]);
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += tstacki[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]];
                    }                   

                    ip_jp_energy = IGINF(tstacki[sequence[jp]][sequence[ip]]
                                          [sequence[jp+1]][sequence[ip-1]]);

                    if (parsi_tstacki == T99 || parsi_tstacki == PARSI)
                    {
                        // actually, just use 3 parameters instead of the tstacki table
                        if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                            ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                        {
                            //internal_AU_closure includes terminal_AU_penalty
                            //sprintf (type, "misc.terminal_AU_penalty");
                            //index = structure_type_index (type);
                            //counter[index]++;                        
                            sprintf (type, "misc.internal_AU_closure");
                            index = structure_type_index (type);
                            counter[index]++;
                            energy += misc.internal_AU_closure;
                            //printf ("In COUNT, add misc.internal_AU_closure + %lf\n", misc.internal_AU_closure);
                        }
                        if (parsi_tstacki == T99)
                        {
                            if ((sequence[jp+1] == A && sequence[ip-1] == G) || (sequence[jp+1] == G && sequence[ip-1] == A))
                            {
                                sprintf (type, "misc.internal_GA_AG_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_GA_AG_mismatch;                                   
                            }
                        }
                        else if (parsi_tstacki == PARSI)
                        {
                            if (sequence[jp+1] == A && sequence[ip-1] == G)
                            {
                                sprintf (type, "misc.internal_AG_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_AG_mismatch;
                                //printf ("In COUNT, add misc.internal_AG_mismatch + %lf\n", misc.internal_AG_mismatch);
                            }
                            if (sequence[jp+1] == G && sequence[ip-1] == A)
                            {
                                sprintf (type, "misc.internal_GA_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_GA_mismatch;
                                //printf ("In COUNT, add misc.internal_GA_mismatch + %lf\n", misc.internal_GA_mismatch);
                            }
                            if (sequence[ip-1] == G && sequence[jp+1] == G)
                            {
                                sprintf (type, "misc.internal_GG_mismatch");
                                index = structure_type_index (type);
                                counter[index]++;
                                energy += misc.internal_GG_mismatch;
                                //printf ("In COUNT, add misc.internal_GG_mismatch + %lf\n", misc.internal_GG_mismatch);
                            }
                        }
                        
                        if (sequence[ip-1] == U && sequence[jp+1] == U)
                        {
                            sprintf (type, "misc.internal_UU_mismatch");
                            index = structure_type_index (type);
                            counter[index]++;
                            energy += misc.internal_UU_mismatch;
                            //printf ("In COUNT, add misc.internal_UU_mismatch + %lf\n", misc.internal_UU_mismatch);
                        }
                    }
                    else if (parsi_tstacki == LAVISH)
                    {
                        // the full tstacki table
                        // SHOULD NOT apply rule 1 here, they are separate parameters
                        sprintf (type, "tstacki[%d][%d][%d][%d]", sequence[jp], sequence[ip], sequence[jp+1], sequence[ip-1]);
                        index = structure_type_index (type);
                        counter[index]++;
                        energy += tstacki[sequence[jp]][sequence[ip]][sequence[jp+1]][sequence[ip-1]];
                    }
                }
            }
        }
    }
    
    
    // now check
    PARAMTYPE energy2 = get_energy (i, j, ip, jp, sequence, NULL);    
    
    if (fabs (energy-energy2) > 0.1)
    {
        printf ("ERROR! The way I compute get_energy and the way I count in s_internal_loop.cpp is different!\n");       
#ifdef DOUBLEPARAMS
        printf ("By counts energy is %.2lf, by get_energy is %.2lf\n", energy, energy2);
#elif LDOUBLEPARAMS
        printf ("By counts energy is %.2Lf, by get_energy is %.2Lf\n", energy, energy2);
#endif
        for (int myi=i; myi <= ip; myi++)    printf ("%c", int_to_nuc(sequence[myi]));
        printf (" ");
        for (int myi=jp; myi <= j; myi++)    printf ("%c", int_to_nuc(sequence[myi]));
        printf (" size1 = %d, size2 = %d, i = %d, j = %d, ip = %d, jp = %d\n", branch1, branch2, i, j, ip, jp);
        exit(1);
    }
//     else
//     {
//         printf ("Counts and energy in s_internal_loop.cpp are equal! %.2lf\n", energy);
//     }
    
}


