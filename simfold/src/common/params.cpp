/***************************************************************************
                          params.cpp  -  description
                             -------------------
    begin                : Sun May 8 2005
    copyright            : (C) 2005 by Mirela Andronescu
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

#include "params.h"
#include "common.h"
#include "externs.h"
#include "constants.h"
#include "string.h"
#include "s_specific_functions.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_internal_loop.h"
#include "s_min_folding.h"
#include "s_sub_folding.h"
#include "s_energy_matrix.h"
#include "simfold.h"

// Ian Wark August 2 2017
// Compilation doesn't like when you directly put "" into a char * argument of a function
// give it this instead
char empty_char[1] = "";

// I played with the following functions to study a bit the parameters

int ignore_AU_penalty = 0;

void print_stacking_energies()
// prints the stacking energies
{
  int i,j,k,l;
  for (i=0; i < NUCL; i++)
      for (j=0; j < NUCL; j++)
          for (k=0; k < NUCL; k++)
          {
              for (l=0; l < NUCL; l++)
                  printf ("stack(%c,%c,%c,%c) = %d\n", int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), int_to_nuc(l), stack[i][j][k][l]);
              printf ("----\n");
          }
}



void print_tstacki_dangling_energies()
// prints the tstacki and dangling energies
{
    int met = 0;
    int i,j,k,l;
    for (k=0; k < NUCL; k++)
        for (l=0; l < NUCL; l++)
        {
            met = 0;
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                {
                    if (tstacki[i][j][k][l] < INF)
                    {
                        met = 1;
                        printf ("tstacki(%c,%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), int_to_nuc(l), tstacki[i][j][k][l]);
                        printf ("d3(%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), dangle_top[i][j][k]);
                        printf ("d5(%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(l), dangle_bot[i][j][l]);
                        printf ("Diff = %5d\n", tstacki[i][j][k][l] - dangle_top[i][j][k] - dangle_bot[i][j][l]);
                    }
                }
            if (met) printf ("---\n");
        }
}


void print_stack_dangling_energies()
// prints the stack and dangling energies
{
    int met = 0;
    int i,j,k,l;
    for (k=0; k < NUCL; k++)
        for (l=0; l < NUCL; l++)
        {
            met = 0;
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                {
                    if (stack[i][j][k][l] < INF)
                    {
                        met = 1;
                        printf ("stack(%c,%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), int_to_nuc(l), stack[i][j][k][l]);
                        printf ("d3(%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), dangle_top[i][j][k]);
                        printf ("d5(%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(l), dangle_bot[i][j][l]);
                        printf ("Diff = %5d\n", stack[i][j][k][l] - dangle_top[i][j][k] - dangle_bot[i][j][l]);
                    }
                }
            if (met) printf ("---\n");
        }
}


void print_stack_equation_dangling()
// prints the stack and dangling energies
{
    int met = 0;
    int i,j,k,l;
    float x,y,z,t;
    x = 0.0592;
    y = -7.0973;
    z = -0.4545;
    t = 7.7252;
    for (k=0; k < NUCL; k++)
        for (l=0; l < NUCL; l++)
        {
            met = 0;
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                {
                    if (stack[i][j][k][l] < INF)
                    {
                        met = 1;
                        printf ("stack(%c,%c,%c,%c) = %5d\t", int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), int_to_nuc(l), stack[i][j][k][l]);
                        printf (" eq = %5.2lf\n", x*dangle_top[i][j][k] + y*dangle_bot[i][j][l] + z*dangle_top[l][k][j] + t*dangle_bot[l][k][i]);
                    }
                }
            if (met) printf ("---\n");
        }
}

void print_int22_tstacki()
// prints the int22 and tstacki energies
{
    int met = 0;
    int i,j,k,l,m,n,o,p;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                            for (o=0; o < NUCL; o++)
                                for (p=0; p < NUCL; p++)
                                {
                                    if (int22[i][j][k][l][m][n][o][p] < INF)
                                    {
                                        printf ("int22(%c,%c,%c,%c,%c,%c,%c,%c) = %5d -- %5d = \t", int_to_nuc(i), int_to_nuc(j),
                                            int_to_nuc(k), int_to_nuc(l), int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(p),
                                            int22[i][j][k][l][m][n][o][p], tstacki[i][j][k][l]+tstacki[n][m][p][o]);
                                        printf ("tstacki(%c,%c,%c,%c) + tstacki(%c,%c,%c,%c) = %5d + %5d\n",
                                            int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), int_to_nuc(l),
                                            int_to_nuc(n), int_to_nuc(m), int_to_nuc(p), int_to_nuc(o), tstacki[i][j][k][l], tstacki[n][m][p][o]);
                                    }
                                }
}


void print_int22()
// prints the int22 and tstacki energies
{
    int met = 0;
    int i,j,k,l,m,n,o,p;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            if (can_pair(i,j))
            {
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                if (can_pair(m,n))
                                {
                                    for (o=0; o < NUCL; o++)
                                        for (p=0; p < NUCL; p++)
                                        {
                                            // exclude duplicates
                                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                                n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                            {
                                                printf ("%.2lf\n",
                                                //int22[%c][%c][%c][%c][%c][%c][%c][%c]=
                                                //int_to_nuc(i), int_to_nuc(j), int_to_nuc(k), int_to_nuc(l),
                                                //int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(p),
                                                int22[i][j][k][l][m][n][o][p]/100.0);
                                            }
                                        }
                                }
            }
}


void test_int22_differences()
// there are some differences here, because in the 53-param model I didn't include all the measured internal loops
{
    PARAMTYPE sint  [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
    int i,j,k,l,m,n,o,p;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                            for (o=0; o < NUCL; o++)
                                for (p=0; p < NUCL; p++)
                                {
                                    sint[i][j][k][l][m][n][o][p] = int22[i][j][k][l][m][n][o][p];
                                    int22[i][j][k][l][m][n][o][p] = INF;
                                }


    simple_internal_energy = 0;
    //save_parameters("turner_parameters_fmnew.txt");
    num_params = create_string_params();

    fill_data_structures_with_new_parameters("turner_parameters_fmnew.txt");

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                            for (o=0; o < NUCL; o++)
                                for (p=0; p < NUCL; p++)
                                {
                                    if (sint[i][j][k][l][m][n][o][p] != int22[i][j][k][l][m][n][o][p])
                                        //printf ("int22[%d][%d][%d][%d][%d][%d][%d][%d]  tur=%d, com=%d, diff=%d\n",
                                        //          i,j,k,l,m,n,o,p, sint[i][j][k][l][m][n][o][p], int22[i][j][k][l][m][n][o][p], int22[i][j][k][l][m][n][o][p]-sint[i][j][k][l][m][n][o][p]);
                                        printf ("int22[5'%c%c%c%c/3'%c%c%c%c]    T99=%3d    MAextra=%3d    T99-MAextra=%d\n",
                                                  int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(j), int_to_nuc(l), int_to_nuc(p), int_to_nuc(n), sint[i][j][k][l][m][n][o][p], int22[i][j][k][l][m][n][o][p], sint[i][j][k][l][m][n][o][p]-int22[i][j][k][l][m][n][o][p]);
                                }
}

void test_tstacki_differences()
// this is fine, no differences here
{
    PARAMTYPE sint  [NUCL] [NUCL] [NUCL] [NUCL];
    int i,j,k,l;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    sint[i][j][k][l] = tstacki[i][j][k][l];
                }


    simple_internal_energy = 0;
    //save_parameters("turner_parameters_fmnew.txt");
    num_params = create_string_params();

    fill_data_structures_with_new_parameters("turner_parameters_fmnew.txt");

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (sint[i][j][k][l] != tstacki[i][j][k][l])
                        printf ("tstacki[%d][%d][%d][%d]  tur=%d, com=%d, diff=%d\n",
                                    i,j,k,l, sint[i][j][k][l], tstacki[i][j][k][l], sint[i][j][k][l]-tstacki[i][j][k][l]);
                }
}


void test_int21_differences()
// they are the same as turner's, within 0.05 kcal/mol
{
    PARAMTYPE sint  [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
    int i,j,k,l,m,n,o;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                            for (o=0; o < NUCL; o++)
                            {
                                sint[i][j][k][l][m][n][o] = int21[i][j][k][l][m][n][o];
                            }

    fill_data_structures_with_new_parameters("turner_parameters_fmnew.txt");

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                            for (o=0; o < NUCL; o++)
                            {
                                if (sint[i][j][k][l][m][n][o] != int21[i][j][k][l][m][n][o])
                                    printf ("int21[%d][%d][%d][%d][%d][%d][%d]  tur=%d, com=%d, diff=%d\n",
                                                i,j,k,l,m,n,o, sint[i][j][k][l][m][n][o], int21[i][j][k][l][m][n][o], int21[i][j][k][l][m][n][o]-sint[i][j][k][l][m][n][o]);
                            }
}


void test_int11_differences()
// they are the same as turner's, within 0.05 kcal/mol
{
    PARAMTYPE sint  [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
    int i,j,k,l,m,n,o;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            sint[i][j][k][l][m][n] = int11[i][j][k][l][m][n];
                            int11[i][j][k][l][m][n] = INF;
                        }

    fill_data_structures_with_new_parameters("turner_parameters_fmnew.txt");

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            if (sint[i][j][k][l][m][n] != int11[i][j][k][l][m][n])
                                printf ("int11[%d][%d][%d][%d][%d][%d]  tur=%d, com=%d, diff=%d\n",
                                            i,j,k,l,m,n, sint[i][j][k][l][m][n], int11[i][j][k][l][m][n], int11[i][j][k][l][m][n]-sint[i][j][k][l][m][n]);
                        }
}

/////////////////////////////////////////////////////////////
// next functions are for playing with parameter learning. //
/////////////////////////////////////////////////////////////


void set_up_index_int_12_34 (char *type, int start0, int start1, int start2, int start3, int &start)
// return newstart
{
    if (type[12] == '0' || (type[6] == '_' && type[34] == '0'))  //int21[1][x][0]
        start = start0;
    if (type[12] == '1' || (type[6] == '_' && type[34] == '1'))  //int21[1][x][1]
        start = start1;
    if (type[12] == '2' || (type[6] == '_' && type[34] == '2'))  //int21[1][x][2]
        start = start2;
    if (type[12] == '3' || (type[6] == '_' && type[34] == '3'))  //int21[1][x][3]
        start = start3;
}

int structure_type_index (const char input[])
  // Mirela: Nov 23, 2003
  // Given the type as a string, return the index in string_params
// TO EXTEND
{
  // Ian Wark July 7 2017
  // A kind of hacky way of supressing the warnings
  // structure_type_index is called in mulitple places with just a string,
  // ex. strucute_type_index("blugeA"), so the function should have const char* as input
  // It needs to be normal char*, so instead just copy the const char* input into a normal char*
  char type[strlen(input)+1];
  strncpy(type, input, strlen(input)+1);

  int i, found;
  found = 0;

  // if type is internal/hairpin/bulge penalty by size, which is bigger than MAXLOOP_X,
  //    then type is actually MAXLOOP_X
  int size;
  // do it later, only if this is penalty, for faster runtime. This actually helps a bit.


  // just a hack to make this faster than linear. NOT GOOD IF THE (NUMBER OF) PARAMETERS CHANGES!!!
    int start = 0;

    /*
    #if (MODEL == SIMPLE)
    if (type[0] == 's')     // stack
        start = 0;
    else if (type[0] == 't')     // tstack or tloop
    {
        if (type[1] == 's')  // tstack
        {
            if (type[8] == '0')   // tstackh[0]
                start = 21;
            else if (type[8] == '1')    //tstackh[1]
                start = 37;
            else if (type[8] == '2')    //tstackh[2]
                start = 53;
            else if (type[8] == '3')    //tstackh[3]
                start = 85;
        }
        else if (type[1] == 'l')    // tloop
            start = 333;
    }
    else if (type[0] == 'm')    // misc....
    {
        if (type[5] == 'i')
        {
            if (type[13] == '_')    // misc.internal_
                start = 117;
            else if (type[13] == '1')   // misc.internal11_
                start = 151;
            else if (type[13] == '2')
            {
                if (type[15] == '1')     // misc.internal21_
                    start = 205;
                else if (type[15] == '2')   // misc.internal22_
                    start = 255;
            }
        }
        else    // al other misc
            start = 324;
    }
    else if (type[0] == 'i')    // int11, int21, int22 or internal_penalty_by_size
    {
        if (type[3] == '1')  // int11
            start = 120;
        else if (type[3] == '2')    // int21 or int22
        {
            if (type[4] == '1')  // int21
            {
                if (type[6] == '1')     // int21[1]...
                    start = 153;
                else if (type[6] == '2')     // int21[2]...
                    start = 179;
            }
            else if (type[4] == '2')  // int22
            {
                if (type[6] == '0')     // int22[0]...
                    start = 207;
                else if (type[6] == '1')     // int22[1]...
                    start = 219;
                else if (type[6] == '2')     // int22[2]...
                    start = 231;
                else if (type[6] == '3')     // int22[3]...
                    start = 243;
            }
        }
        else if (type[3] == 'e')    // internal_penalty_by_size
        {
            start = 308;
            // next internal
            if (sscanf (type, "internal_penalty_by_size[%d]", &size))
            {
                if (size > MAXLOOP_I)
                {
                    //printf ("size = %d\n", size);
                    sprintf (type, "internal_penalty_by_size[%d]", MAXLOOP_I);
                }
            }
        }
    }
    else if (type[0] == 'd')    // dangle_top or dangle_bot
    {
        if (type[7] == 't')  // dangle_top
            start = 260;
        else if (type[7] == 'b')  // dangle_bot
            start = 284;
    }
    else if (type[0] == 'b')    // bulge_penalty_by_size
    {
        start = 311;
        // next bulge
        if (sscanf (type, "bulge_penalty_by_size[%d]", &size))
        {
            if (size > MAXLOOP_B)
            {
                //printf ("size = %d\n", size);
                sprintf (type, "bulge_penalty_by_size[%d]", MAXLOOP_B);
            }
        }
    }
    else if (type[0] == 'h')    // hairpin_penalty_by_size
    {
        start = 317;
        // first hairpin
        if (sscanf (type, "hairpin_penalty_by_size[%d]", &size))
        {
            if (size > MAXLOOP_H)
            {
                //printf ("size = %d\n", size);
                sprintf (type, "hairpin_penalty_by_size[%d]", MAXLOOP_H);
            }
        }
    }
    #elif (MODEL == EXTENDED)
    */
                // next internal
    int end;
    if (parsi_length == T99)               end = MAXLOOP_I_T99;
    else if (parsi_length == PARSI || parsi_length == ZL)        end = MAXLOOP_I_PARSI;
    else if (parsi_length == LAVISH)       end = MAXLOOP_I_LAVISH;
    if (sscanf (type, "internal_penalty_by_size[%d]", &size))
    {
        if (size > end)
        {
            sprintf (type, "internal_penalty_by_size[%d]", end);
        }
    }
    if (parsi_length == T99)               end = MAXLOOP_B_T99;
    else if (parsi_length == PARSI || parsi_length == ZL)        end = MAXLOOP_B_PARSI;
    else if (parsi_length == LAVISH)       end = MAXLOOP_B_LAVISH;
    if (sscanf (type, "bulge_penalty_by_size[%d]", &size))
    {
        if (size > end)
        {
                //printf ("size = %d\n", size);
            sprintf (type, "bulge_penalty_by_size[%d]", end);
        }
    }
    if (parsi_length == T99)               end = MAXLOOP_H_T99;
    else if (parsi_length == PARSI || parsi_length == ZL)        end = MAXLOOP_H_PARSI;
    else if (parsi_length == LAVISH)       end = MAXLOOP_H_LAVISH;
    if (sscanf (type, "hairpin_penalty_by_size[%d]", &size))
    {
        if (size > end)
        {
                //printf ("size = %d\n", size);
            sprintf (type, "hairpin_penalty_by_size[%d]", end);
        }
    }

    // instead of giving the raw numbers, first figure out each category once,
    //        by calling the function set_starters
    if (type[0] == 's')     // stack or special_hl
    {
        if (type[1] == 't')   // stack
            start = 0;
        else if (type[1] == 'p')    // special_hl
        {
            start = start_special_hl;
        }
    }

    else if (type[0] == 'm')    // misc
    {
        if (type[5] == 'h')     // misc.hairpin_AU etc or misc.hairpin_GGG etc
        {
            if (type[13] == 'c' || (type[13]=='G' && type[14]=='G' && type[15]=='G'))   // misc.hairpin_GGG etc
                start = start_misc_last;
            else            start = 21;     // misc.hairpin_AU etc
        }
        else if (type[5] == 'i')    // misc.internal or misc.intermolecular
        {
            if (type[10] == 'n')    // misc.internal
            {
                if (type[14] == 's')    // misc.internal_special
                    start = start_misc_last;
                else if (type[13] == '_')    // misc.internal_AU_closure etc
                    start = start_internal;
                else if (type[13] == '1')   // misc.internal11_AU_closure etc
                    start = start_internal11;
                else if (type[13] == '2')   // misc.internal21 and misc.internal22
                {
                    if (type[14] == '1')    // misc.internal21
                        start = start_internal21;
                    else if (type[14] == '2')   // misc.internal22
                        start = start_internal22;
                }
            }
            else   start = start_misc_last;                 // misc.intermolecular
        }
        else   start = start_misc_last;     // misc.terminal_AU_penalty, misc.multi etc.
    }
    else if (type[0] == 't')    // tstackh or tstacki
    {
        if (type[6] == 'h')         start = 21;     // tstackh
        else if (type[6] == 'i')    start = start_internal;   // tstacki      // same as misc.internal_AU_closure

    }
    else if (type[0] == 'b')    // bulgeA-U, bulge1 or bulge_penalty_by_size
    {
        if (type[5] == '_')  start = start_bulge_size;  // bulge_penalty_by_size
        else                 start = start_bulge;   // bulgeA-U or bulge1
    }
    else if (type[0] == 'h')   start = start_hairpin_size;   // hairpin_penalty_by_size

    else if (type[0] == 'i')    // internal_penalty_by_size, internal_asymmetry, int11, int21, int22,
    {
        if (type[3] == 'e')     // internal_penalty_by_size or internal_asymmetry
            start = start_internal_size;
        else if (type[3] == '1')     // int11
        {
            if (type[6] == '0' || (type[6] == '_' && type[28] == '0'))  //int11[0]
                start = start_internal11;
            else if (type[6] == '1' || (type[6] == '_' && type[28] == '1'))  //int11[1]
                start = start_internal11_C;
            else if (type[6] == '2' || (type[6] == '_' && type[28] == '2'))  //int11[2]
                start = start_internal11_G;
            else if (type[6] == '3' || (type[6] == '_' && type[28] == '3'))  //int11[3]
                start = start_internal11_U;
        }
        else if (type[3] == '2')    // int21 or int22
        {
            if (type[4] == '1')     // int21
            {
                if (type[6] == '0' || (type[6] == '_' && type[28] == '0'))  //int21[0]
                {
                    set_up_index_int_12_34 (type, start_internal21_AUA, start_internal21_AUC,
                            start_internal21_AUG, start_internal21_AUU, start);
                }
                if (type[6] == '1' || (type[6] == '_' && type[28] == '1'))  //int21[1]
                {
                    set_up_index_int_12_34 (type, start_internal21_CGA, start_internal21_CGC,
                            start_internal21_CGG, start_internal21_CGU, start);
                }
                if (type[6] == '2' || (type[6] == '_' && type[28] == '2'))  //int21[2]
                {
                    if (type[9] == '1' || (type[6] == '_' && type[31] == '1'))  //int21[2][1]
                    {
                        set_up_index_int_12_34 (type, start_internal21_GCA, start_internal21_GCC,
                                start_internal21_GCG, start_internal21_GCU, start);
                    }
                    else if (type[9] == '3' || (type[6] == '_' && type[31] == '3'))  //int21[2][3]
                    {
                        set_up_index_int_12_34 (type, start_internal21_GUA, start_internal21_GUC,
                                start_internal21_GUG, start_internal21_GUU, start);
                    }
                }
                if (type[6] == '3' || (type[6] == '_' && type[28] == '3'))  //int21[3]
                {
                    if (type[9] == '0' || (type[6] == '_' && type[31] == '0'))  //int21[3][0]
                    {
                        set_up_index_int_12_34 (type, start_internal21_UAA, start_internal21_UAC,
                                start_internal21_UAG, start_internal21_UAU, start);
                    }
                    else if (type[9] == '2' || (type[6] == '_' && type[31] == '2'))  //int21[3][2]
                    {
                        set_up_index_int_12_34 (type, start_internal21_GUA, start_internal21_UGC,
                                start_internal21_UGG, start_internal21_UGU, start);
                    }
                }
            }
            else                    // int22
            {
                if (type[6] == '0' || (type[6] == '_' && type[28] == '0'))  //int22[0]
                {
                    set_up_index_int_12_34 (type, start_internal22_AUA, start_internal22_AUC,
                            start_internal22_AUG, start_internal22_AUU, start);
                }
                if (type[6] == '1' || (type[6] == '_' && type[28] == '1'))  //int22[1]
                {
                    set_up_index_int_12_34 (type, start_internal22_CGA, start_internal22_CGC,
                            start_internal22_CGG, start_internal22_CGU, start);
                }
                if (type[6] == '2' || (type[6] == '_' && type[28] == '2'))  //int22[2]
                {
                    if (type[9] == '1' || (type[6] == '_' && type[31] == '1'))  //int22[2][1]
                    {
                        set_up_index_int_12_34 (type, start_internal22_GCA, start_internal22_GCC,
                                start_internal22_GCG, start_internal22_GCU, start);
                    }
                    else if (type[9] == '3' || (type[6] == '_' && type[31] == '3'))  //int22[2][3]
                    {
                        set_up_index_int_12_34 (type, start_internal22_GUA, start_internal22_GUC,
                                start_internal22_GUG, start_internal22_GUU, start);
                    }
                }
                if (type[6] == '3' || (type[6] == '_' && type[28] == '3'))  //int22[3]
                {
                    if (type[9] == '0' || (type[6] == '_' && type[31] == '0'))  //int22[3][0]
                    {
                        set_up_index_int_12_34 (type, start_internal22_UAA, start_internal22_UAC,
                                start_internal22_UAG, start_internal22_UAU, start);
                    }
                    else if (type[9] == '2' || (type[6] == '_' && type[31] == '2'))  //int22[3][2]
                    {
                        set_up_index_int_12_34 (type, start_internal22_GUA, start_internal22_UGC,
                                start_internal22_UGG, start_internal22_UGU, start);
                    }
                }
            }
        }
    }   // end if starts with i
    else if (type[0] == 'd')   start = start_dangle;     // dangle


    //printf ("start = %d, num_params = %d, type = %s\n", start, num_params, type);
  // now traverse all params, to see which is the right index
  // veeery slow!!! // it's faster now, since I added the above
  // for the EXTENDED model I added some sort of cashing system
    for (i=start; i < num_params; i++)
    {
        //printf ("type = %s, string_params[%d] = %s\n", type, i, string_params[i]);
        if (strcmp (type, string_params[i]) == 0)
        // Note: this is very slow
        {
            found = 1;
            //printf ("%s found in %d steps\n", type, i-start+1);
            break;
        }
    }
    if (!found)
    {
        fprintf (stderr, "!!!TYPE NOT found: %s!!!\n", type);
        exit(1);
    }
  return i;
}



void set_starters ()
// set all the starters to the appropriate values, to be used in structure_type_index
{
    if (parsi_bulge1 != T99)
        start_bulge = structure_type_index("bulgeA");
    start_internal = structure_type_index("misc.internal_AU_closure");
    if (parsi_int11 == T99)
        start_internal11 = structure_type_index("int11[0][3][3][3][0][3]");
    else
        start_internal11 = structure_type_index("misc.internal11_AU_closure");
    if (parsi_int11 == LAVISH)
    {
        start_internal11_C = structure_type_index("int11[1][2][0][0][0][3]");
        start_internal11_G = structure_type_index("int11[2][1][0][0][0][3]");
        //start_internal11_U = structure_type_index("int11_experimental_addition[3][0][0][0][0][3]");
        start_internal11_U = structure_type_index("int11[3][0][0][0][0][3]");
    }

    if (parsi_int21 == T99)
        start_internal21 = structure_type_index("int21[1][2][0][0][1][2][0]");
    else
        start_internal21 = structure_type_index("misc.internal21_initiation");
    if (parsi_int21 == LAVISH)
    {
        start_internal21_AUA = structure_type_index("int21[0][3][0][0][0][3][0]");
        start_internal21_AUC = structure_type_index("int21[0][3][1][0][0][3][0]");
        start_internal21_AUG = structure_type_index("int21[0][3][2][0][0][3][0]");
        start_internal21_AUU = structure_type_index("int21[0][3][3][0][0][3][0]");
        //start_internal21_CGA = structure_type_index("int21_experimental_addition[1][2][0][0][0][3][0]");
        start_internal21_CGA = structure_type_index("int21[1][2][0][0][0][3][0]");
        start_internal21_CGC = structure_type_index("int21[1][2][1][0][0][3][0]");
        start_internal21_CGG = structure_type_index("int21[1][2][2][0][0][3][0]");
        start_internal21_CGU = structure_type_index("int21[1][2][3][0][0][3][0]");
        start_internal21_GCA = structure_type_index("int21[2][1][0][0][0][3][0]");
        start_internal21_GCC = structure_type_index("int21[2][1][1][0][0][3][0]");
        start_internal21_GCG = structure_type_index("int21[2][1][2][0][0][3][0]");
        start_internal21_GCU = structure_type_index("int21[2][1][3][0][0][3][0]");
        start_internal21_GUA = structure_type_index("int21[2][3][0][0][0][3][0]");
        start_internal21_GUC = structure_type_index("int21[2][3][1][0][0][3][0]");
        start_internal21_GUG = structure_type_index("int21[2][3][2][0][0][3][0]");
        start_internal21_GUU = structure_type_index("int21[2][3][3][0][0][3][0]");
        //start_internal21_UAA = structure_type_index("int21_experimental_addition[3][0][0][0][0][3][0]");
        start_internal21_UAA = structure_type_index("int21[3][0][0][0][0][3][0]");
        start_internal21_UAC = structure_type_index("int21[3][0][1][0][0][3][0]");
        start_internal21_UAG = structure_type_index("int21[3][0][2][0][0][3][0]");
        start_internal21_UAU = structure_type_index("int21[3][0][3][0][0][3][0]");
        start_internal21_UGA = structure_type_index("int21[3][2][0][0][0][3][0]");
        start_internal21_UGC = structure_type_index("int21[3][2][1][0][0][3][0]");
        start_internal21_UGG = structure_type_index("int21[3][2][2][0][0][3][0]");
        start_internal21_UGU = structure_type_index("int21[3][2][3][0][0][3][0]");
    }

    if (parsi_int22 == T99)
        start_internal22 = structure_type_index("int22[0][3][0][0][3][0][0][0]");
    else
        start_internal22 = structure_type_index("misc.internal22mid_group1");
    if (parsi_int22 == LAVISH)
    {
        start_internal22_AUA = structure_type_index("int22[0][3][0][0][0][3][0][0]");
        start_internal22_AUC = structure_type_index("int22[0][3][1][0][0][3][0][0]");
        start_internal22_AUG = structure_type_index("int22[0][3][2][0][0][3][0][0]");
        start_internal22_AUU = structure_type_index("int22[0][3][3][0][0][3][0][0]");
        start_internal22_CGA = structure_type_index("int22[1][2][0][0][0][3][0][0]");
        start_internal22_CGC = structure_type_index("int22[1][2][1][0][0][3][0][0]");
        start_internal22_CGG = structure_type_index("int22[1][2][2][0][0][3][0][0]");
        start_internal22_CGU = structure_type_index("int22[1][2][3][0][0][3][0][0]");
        start_internal22_GCA = structure_type_index("int22[2][1][0][0][0][3][0][0]");
        start_internal22_GCC = structure_type_index("int22[2][1][1][0][0][3][0][0]");
        start_internal22_GCG = structure_type_index("int22[2][1][2][0][0][3][0][0]");
        start_internal22_GCU = structure_type_index("int22[2][1][3][0][0][3][0][0]");
        start_internal22_GUA = structure_type_index("int22[2][3][0][0][0][3][0][0]");
        start_internal22_GUC = structure_type_index("int22[2][3][1][0][0][3][0][0]");
        start_internal22_GUG = structure_type_index("int22[2][3][2][0][0][3][0][0]");
        start_internal22_GUU = structure_type_index("int22[2][3][3][0][0][3][0][0]");
        //start_internal22_UAA = structure_type_index("int22_experimental_addition[3][0][0][0][0][3][0][0]");
        start_internal22_UAA = structure_type_index("int22[3][0][0][0][0][3][0][0]");
        //start_internal22_UAC = structure_type_index("int22_experimental_addition[3][0][1][0][0][3][0][1]");
        start_internal22_UAC = structure_type_index("int22[3][0][1][0][0][3][0][1]");
        //start_internal22_UAG = structure_type_index("int22_experimental_addition[3][0][2][0][0][3][0][2]");
        start_internal22_UAG = structure_type_index("int22[3][0][2][0][0][3][0][2]");
        start_internal22_UAU = structure_type_index("int22[3][0][3][0][0][3][0][3]");
        //start_internal22_UGA = structure_type_index("int22_experimental_addition[3][2][0][0][2][3][0][0]");
        start_internal22_UGA = structure_type_index("int22[3][2][0][0][2][3][0][0]");
        //start_internal22_UGC = structure_type_index("int22_experimental_addition[3][2][1][0][2][3][0][1]");
        start_internal22_UGC = structure_type_index("int22[3][2][1][0][2][3][0][1]");
        //start_internal22_UGG = structure_type_index("int22_experimental_addition[3][2][2][0][2][3][0][2]");
        start_internal22_UGG = structure_type_index("int22[3][2][2][0][2][3][0][2]");
        start_internal22_UGU = structure_type_index("int22[3][2][3][0][2][3][0][3]");
    }

    if (parsi_dangles == T99 | parsi_dangles == LAVISH)
        start_dangle = structure_type_index("dangle_top[0][3][0]");
    start_internal_size = structure_type_index("internal_penalty_by_size[4]");
    if (parsi_bulge1 == T99)    // there's no bulge1, so bulge_penalty_by_size starts from 1
        start_bulge_size = structure_type_index("bulge_penalty_by_size[1]");
    else
        start_bulge_size = structure_type_index("bulge_penalty_by_size[2]");
    start_hairpin_size = structure_type_index("hairpin_penalty_by_size[3]");
    start_misc_last = structure_type_index("misc.terminal_AU_penalty");
    if (parsi_special == LAVISH || parsi_special == T99_LAVISH)
        start_special_hl = structure_type_index("special_hl[0].energy");
}

void count_AU_penalty (int base_i, int base_j, double *counter)
// PRE:  base_i and base_j make a pair
// POST: increment counter at the appropriate position
// Mirela: Nov 23, 2003
{
    int index;
    if ((base_i != C || base_j != G) &&
        (base_i != G || base_j != C))
      {
        index = structure_type_index ("misc.terminal_AU_penalty");
        counter[index]++;
      }
}


void count_dangling_energy (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter)
//      (   )...(   )
//      i1..i2..i3..i4
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling energy between i2 and i3
// Mirela: Nov 23, 2003
// Feb 28, 2008. We might have a situation like this: <   >...(   ) or like this: (   )...<   >.
//  In that case, only add the parameter dangling onto the () pair, if it's at least 2 unpaired bases away
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    char type[100];
    int index_top, index_bot;
    int first_index;    // first index should be dangle_top[0][3][0]
    first_index = structure_type_index ("dangle_top[0][3][0]");

    d_top = 0;
    d_bot = 0;

    if (i2 != link && structure[i2] != '>')
      {
        //d_top = MIN (0, IGINF(dangle_top[sequence[i2]][sequence[i1]][sequence[i2+1]]));
        d_top = dangle_top[sequence[i2]] [sequence[i1]] [sequence[i2+1]];
        sprintf (type, "dangle_top[%d][%d][%d]",sequence[i2], sequence[i1], sequence[i2+1]);
        index_top = structure_type_index (type);
      }
    if (i3-1 != link && structure[i3] != '<')
      {
        //d_bot = MIN (0, IGINF(dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]]));
        d_bot = dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]];
        sprintf (type, "dangle_bot[%d][%d][%d]",sequence[i4], sequence[i3], sequence[i3-1]);
        index_bot = structure_type_index (type);
      }

    if (structure[i2] == '>' && structure[i3] == '(')   // pseudoknot, ignore dangling end dangling on it
    {
        if (i3 <= i2+2)     // >.( or >(   ignore completely
            energy = 0;
        else                // >...(
        {
            energy = d_bot;
            counter[index_bot]++;
        }
    }
    else if (structure[i2] == ')' && structure[i3] == '<')   // pseudoknot, ignore dangling end dangling on it
    {
        if (i3 <= i2+2)     // ).< or )<   ignore completely
            energy = 0;
        else                // )...<
        {
            energy = d_top;
            counter[index_top]++;
        }
    }
    else if (structure[i2] == '>' && structure[i3] == '<')  // case >..<  ignore completely
    {
        energy = 0;
    }
    else if (i2+1 == i3-1 && i2 == link)
      {
        energy = d_bot;
        counter[index_bot]++;
      }
    else if (i2+1 == i3-1 && i3-1 == link)
      {
        energy = d_top;
        counter[index_top]++;
      }
    else if (i2+1 == i3-1)     // see which is smaller
    {
        //energy = d_top < d_bot ? d_top : d_bot;
        // NOTE: the comparison of d_top with d_bot is not right!
        // NO! This is not right if we don't know which of d_top and d_bot is smaller

        // if we restrict the 3' dangling ends to be less than the 5' ones, then it's ok to do what follows

        // if we fix the dangling ends to the Turner parameters, we have to count this
        //if (d_top < d_bot) counter[index_top]++;
        //else counter[index_bot]++;

        if (simple_dangling_ends)
        {
            energy = d_top;
            counter[index_top]++;
        }
        else
        {
            if (d_top < d_bot)
            {
                energy = d_top;
                counter[index_top]++;
            }
            else
            {
                energy = d_bot;
                counter[index_bot]++;
            }
        }

        // if we introduce another variable as min, we need to do this
        //counter_min_dangle[index_top-first_index][index_bot-first_index]++;
    }
    else if (i2+1 < i3-1)
    {
        energy = d_top + d_bot;
        counter[index_top]++;
        counter[index_bot]++;
    }

    else // if there is no free base between the two branches, return 0
        energy = 0;
    //    return energy;
}

void count_dangling_energy_left (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter)
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i1 and i3
// Mirela: Nov 23, 2003
// Feb 28, 2008. We might have a situation like this:  (....<    >   ). In that case, only add the 3' dangling end.
// If it's (.<    > ), don't add any
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    char type[100];
    int index_top, index_bot;
    d_top = 0;
    d_bot = 0;
    int first_index;    // first index should be dangle_top[0][3][0]
    first_index = structure_type_index ("dangle_top[0][3][0]");

    // this will be used in multi-loops.
    // add the dangle_top, even if it is positive
    if (i1 != link)
      {
        //d_top = MIN (0, IGINF(dangle_top[sequence[i1]] [sequence[i2]] [sequence[i1+1]]));
        d_top = dangle_top[sequence[i1]] [sequence[i2]] [sequence[i1+1]];
        sprintf (type, "dangle_top[%d][%d][%d]",sequence[i1], sequence[i2], sequence[i1+1]);
        index_top = structure_type_index (type);
      }
    // in the other parts of the multi-loop, the dangles are added only if they are negative
    if (i3-1 != link && structure[i3] != '<')
      {
        //d_bot = MIN (0, IGINF(dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]]));
        d_bot = dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]];
        sprintf (type, "dangle_bot[%d][%d][%d]",sequence[i4], sequence[i3], sequence[i3-1]);
        index_bot = structure_type_index (type);
      }

    if (structure[i3] == '<')   // pseudoknot inside, ignore dangling end dangling on it
    {
        if (i3 <= i1+2)     // (< or (.<, ignore completely
            energy = 0;
        else                // (....<
        {
            energy = d_top;
            counter[index_top]++;
        }
    }
    else if (i1+1 == i3-1 && i1 == link)
      {
        energy = d_bot;
        counter[index_bot]++;
      }
    else if (i1+1 == i3-1 && i3-1 == link)
      {
        energy = d_top;
        counter[index_top]++;
      }
    else if (i1+1 == i3-1)     // see which is smaller
    {
        //energy = d_top < d_bot ? d_top : d_bot;
        // NOTE: the comparison of d_top with d_bot is not right!
        //if (d_top < d_bot) counter[index_top]++;
        //else counter[index_bot]++;
        //counter_min_dangle[index_top-first_index][index_bot-first_index]++;

        if (simple_dangling_ends)
        {
            energy = d_top;
            counter[index_top]++;
        }
        else
        {
            if (d_top < d_bot)
            {
                energy = d_top;
                counter[index_top]++;
            }
            else
            {
                energy = d_bot;
                counter[index_bot]++;
            }
        }
    }
    else if (i1+1 < i3-1)
    {
        energy = d_top + d_bot;
        counter[index_top]++;
        counter[index_bot]++;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    //    return energy;
}


void count_dangling_energy_right (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter)
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2
// Mirela: Nov 23, 2003
// Feb 28, 2008. We might have a situation like this:  (    <    >...)
//  In that case, only add the 5' dangling end if it's at least two unpaired bases away
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    char type[100];
    int index_top, index_bot;
    int first_index;    // first index should be dangle_top[0][3][0]
    first_index = structure_type_index ("dangle_top[0][3][0]");

    d_top = 0;
    d_bot = 0;

    if (i4 != link && structure[i3] != '<')
      {
        //d_top = MIN (0, IGINF(dangle_top[sequence[i4]] [sequence[i3]] [sequence[i4+1]]));
        d_top = dangle_top[sequence[i4]] [sequence[i3]] [sequence[i4+1]];
        sprintf (type, "dangle_top[%d][%d][%d]",sequence[i4], sequence[i3], sequence[i4+1]);
        index_top = structure_type_index (type);
      }
    if (i2-1 != link)
      {
        //d_bot = MIN (0, IGINF(dangle_bot[sequence[i1]] [sequence[i2]] [sequence[i2-1]]));
        d_bot = dangle_bot[sequence[i1]] [sequence[i2]] [sequence[i2-1]];
        sprintf (type, "dangle_bot[%d][%d][%d]",sequence[i1], sequence[i2], sequence[i2-1]);
        index_bot = structure_type_index (type);
      }

    if (structure[i4] == '>')   // pseudoknot inside, ignore dangling end dangling on it
    {
        if (i2 <= i4+2)     // >.) or >)   ignore completely
            energy = 0;
        else                // >...)
        {
            energy = d_bot;
            counter[index_bot]++;
        }
    }
    else if (i4+1 == i2-1 && i4 == link)
      {
        energy = d_bot;
        counter[index_bot]++;
      }
    else if (i4+1 == i2-1 && i2-1 == link)
      {
        energy = d_top;
        counter[index_top]++;
      }
    else if (i4+1 == i2-1)     // see which is smaller
    {
        //energy = d_top < d_bot ? d_top : d_bot;
        // NOTE: the comparison of d_top with d_bot is not right!
        //if (d_top < d_bot) counter[index_top]++;
        //else counter[index_bot]++;
        //counter_min_dangle[index_top-first_index][index_bot-first_index]++;

        if (simple_dangling_ends)
        {
            energy = d_top;
            counter[index_top]++;
        }
        else
        {
            if (d_top < d_bot)
            {
                energy = d_top;
                counter[index_top]++;
            }
            else
            {
                energy = d_bot;
                counter[index_bot]++;
            }
        }
    }
    else if (i4+1 < i2-1)
    {
        energy = d_top + d_bot;
        counter[index_top]++;
        counter[index_bot]++;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    //    return energy;
}

void count_penalty_by_size (int size, char type, double *counter)
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop
// Mirela: Nov 23, 2003
{
    char t[100];
    int index;
    PARAMTYPE penaltyMAX, penalty;
    double logval;
    int end;
    // the penalties for size <= MAXLOOP should be read from the file "loop"
    // actually only up to some size are measured, for the others apply the formula (see constants.h)
    if (parsi_length == T99)
    {
        if (type == 'H')    end = MAXLOOP_H_T99;
        if (type == 'B')    end = MAXLOOP_B_T99;
        if (type == 'I')    end = MAXLOOP_I_T99;
    }
    else if (parsi_length == PARSI || parsi_length == ZL)
    {
        if (type == 'H')    end = MAXLOOP_H_PARSI;
        if (type == 'B')    end = MAXLOOP_B_PARSI;
        if (type == 'I')    end = MAXLOOP_I_PARSI;
    }
    else if (parsi_length == LAVISH)
    {
        if (type == 'H')    end = MAXLOOP_H_LAVISH;
        if (type == 'B')    end = MAXLOOP_B_LAVISH;
        if (type == 'I')    end = MAXLOOP_I_LAVISH;
    }

    if (type == 'H' && size <= end)
    {
        //return hairpin_penalty_by_size[size];
        sprintf (t, "hairpin_penalty_by_size[%d]", size);
        index = structure_type_index(t);
        counter[index]++;
        return;
    }
    if (type == 'I' && size <= end)
    {
        //return internal_penalty_by_size[size];
        sprintf (t, "internal_penalty_by_size[%d]", size);
        index = structure_type_index(t);
        counter[index]++;
        return;
    }
    if (type == 'B' && size <= end)
    {
        //return bulge_penalty_by_size[size];
        sprintf (t, "bulge_penalty_by_size[%d]", size);
        index = structure_type_index(t);
        counter[index]++;
        return;
    }

    // Mirela: Nov 23, 2003
    // We keep the formula
    // Extrapolation for large loops based on polymer theory
    // internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30)
    // and imcrement for x_penalty_by_size[MAXLOOP] and for misc.param_greater30

    // size > MAXLOOP _B, _H or _I
    if (type == 'H')
      {
        penaltyMAX = hairpin_penalty_by_size[end];
        logval = log (1.0*size/end);
        sprintf (t, "hairpin_penalty_by_size[%d]", end);
        index = structure_type_index (t);
        counter[index]++;
        // added next line on Dec 23, 2006
        // TODO Mar 18, 2008: should that be round???
        //counter[num_params] += round(1.079*logval*100.0);
        counter[num_params] += misc.param_greater30 * logval * 100.0;
      }
    else if (type == 'I')
      {
        penaltyMAX = internal_penalty_by_size[end];
        logval = log (1.0*size/end);
        sprintf (t, "internal_penalty_by_size[%d]", end);
        index = structure_type_index (t);
        counter[index]++;
        // added next line on Dec 23, 2006
        //printf ("Adding %d to total\n", (int)(1.079*logval*100.0));
        //counter[num_params] += round(1.079*logval*100.0);
        counter[num_params] += misc.param_greater30 * logval * 100.0;
      }
    else
      {
        penaltyMAX = bulge_penalty_by_size[end];
        logval = log (1.0*size/end);
        sprintf (t, "bulge_penalty_by_size[%d]", end);
        index = structure_type_index (t);
        counter[index]++;
        // added next line on Dec 23, 2006
        //counter[num_params] += round(1.079*logval*100.0);
        counter[num_params] += misc.param_greater30 * logval * 100.0;
      }

    // Dec 4, 2005: let's keep misc.param_greater30 fixed, i.e. not consider it a parameter
    //penalty = (PARAMTYPE) (penaltyMAX + misc.param_greater30 * logval);
    //index = structure_type_index ("misc.param_greater30");
    // TODO: to store double instead of int
    //counter[index]+= (int) LogValue[size];

    //    return penalty;
    return;
}

void count_asymmetry_penalty (int size1, int size2, double *counter)
// PRE:  size1 and size2 are the sizes of the two free base groups
// POST: Calculate the asymmetry penalty for internal loops
//       Note that if size1 == size2, pen is 0
// Mirela: Nov 23, 2003
{
    if (parsi_asymmetry == T99)
    {
        // count a fixed value for now
        PARAMTYPE ass, pen;
        ass = misc.asymmetry_penalty_array [MIN (2, MIN (size1, size2))-1];
        pen = MIN (misc.asymmetry_penalty_max_correction, abs (size1-size2) * ass);

        //printf ("Adding asym penalty %d to total\n", pen);
        counter[num_params] += pen;

        // don't count anything for now, considered them fixed.
        // TODO: to add new parameter t, for the min(a,b), and then add the restriction t <= a, t <= b and minimize (-t)
        /*
        int ass, pen, index;
        ass = misc.asymmetry_penalty_array [MIN (2, MIN (size1, size2))-1];
        pen = MIN (misc.asymmetry_penalty_max_correction, abs (size1-size2) * ass);

        // store misc.asymmetry_penalty_array[0] and [1], and misc.asymmetry_penalty_max_correction
        index = structure_type_index("misc.asymmetry_penalty_max_correction"); counter[index]++;
        if (size1-size2 == 1)
        {
            index = structure_type_index("misc.asymmetry_penalty_array[0]"); counter[index]++;
        }
        else
        {
            index = structure_type_index("misc.asymmetry_penalty_array[1]"); counter[index]++;
            // [2] and [3] are never used
        }
        */
        // pen can only be 0.5, 1, 1.5, 2, 2.5 or 3 with the current model;
        /*switch (size1-size2)
        {
            case 1: index = structure_type_index("misc.asymmetry_penalty[1]"); counter[index]++; break;
            case 2: index = structure_type_index("misc.asymmetry_penalty[2]"); counter[index]++; break;
            case 3: index = structure_type_index("misc.asymmetry_penalty[3]"); counter[index]++; break;
            case 4: index = structure_type_index("misc.asymmetry_penalty[4]"); counter[index]++; break;
            case 5: index = structure_type_index("misc.asymmetry_penalty[5]"); counter[index]++; break;
            default:index = structure_type_index("misc.asymmetry_penalty[6]"); counter[index]++; break;
        }*/

        //    return pen;
    }
    else if (parsi_asymmetry == PARSI || parsi_asymmetry == LAVISH)
    {
        int index;
        char type[50];
        // assume the size1 + size2 <= MAXLOOP_I. If it's greater, just use the value of the last parameter
        // first symmetric

        if (size1 == size2)     return;
        if (parsi_asymmetry == PARSI)
        {
            sprintf (type, "internal_asymmetry_initiation");
            index = structure_type_index(type);
            counter[index]++;

            sprintf (type, "internal_asymmetry_slope");
            index = structure_type_index(type);
            counter[index] += log (abs (size1-size2));
        }
        else if (parsi_asymmetry == LAVISH)
        {
            if (abs (size1-size2) < MAXLOOP_ASYM)
            {
                // we assume the following model: from asymmetry 1 to 4, we use initiation, slope and int_asym (like an addition)
                if (abs (size1-size2) <= MAX_EXP_ASYM)
                {
                    sprintf (type, "internal_asymmetry_initiation");
                    index = structure_type_index(type);
                    counter[index]++;

                    sprintf (type, "internal_asymmetry_slope");
                    index = structure_type_index(type);
                    counter[index] += log (abs (size1-size2));

                    sprintf (type, "internal_asymmetry[%d]", abs (size1-size2));
                    index = structure_type_index(type);
                    counter[index]++;
                }
                else
                {
                    sprintf (type, "internal_asymmetry[%d]", abs (size1-size2));
                    index = structure_type_index(type);
                    counter[index]++;
                }
            }
            else
            {
                sprintf (type, "internal_asymmetry_initiation");
                index = structure_type_index(type);
                counter[index]++;

                sprintf (type, "internal_asymmetry_slope");
                index = structure_type_index(type);
                counter[index] += log (abs (size1-size2));
            }
        }
    }
}


double count_types (int link, int *sequence, char *csequence, char *structure, char *restricted, str_features *f, double *counter)
// Mirela: Nov 23, 2003
// PRE: string_params have been filled, i.e. by num_params = calling create_string_params() or num_params = create_building_block_strings()
//        where num_params is the global variable
// instead of summing up the energies of each type, increment the corresponding counter
// the free value is counter[num_params]
// Feb 28, 2008: added the situation when the structure can have <xxxx>, meaning pseudoknot and ignore from the energy model
// return the energy
// IF counter is NULL, then don't compute the counts, just the free energy
{
    int i;
    PARAMTYPE energy, en, AUpen;
    PARAMTYPE dang;
    PARAMTYPE misc_energy;
    int h,l;
    static int cannot_add_dangling[MAXSLEN];
    char type[100];
    int index;

    int nb_nucleotides = strlen(csequence);
    for (i=0; i < nb_nucleotides; i++) cannot_add_dangling[i] = 0;

    int *p_table = NULL;
    if (restricted[0] != '\0')
    {
        p_table = new int[nb_nucleotides];
        detect_original_pairs (restricted, p_table);
    }

    energy = 0;
    AUpen = 0;

    for (i=0; i < nb_nucleotides; i++)
    {
        if (debug)
            printf ("i=%d, f[i].pair=%d, f[i].type=%c\n", i, f[i].pair, f[i].type);
        // add some AU_penalties
        if ((i==0 || (i-1 == link && !cannot_add_dangling[i-1] && f[i-1].pair == -1))
             && f[i].pair > i && structure[i] != '<')
        {
            if (!ignore_AU_penalty)
            {
                AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
                if (debug)
                    printf ("%d - AUpen1 \t- add energy %6g\n", i, AUpen);
                energy += AUpen;
                if (counter != NULL)    count_AU_penalty (sequence[i], sequence[f[i].pair], counter);
            }
        }
        else if ( i > 0 && f[i].pair > i  && structure[i] != '<' && f[i-1].pair < i-1 &&
             f[i-1].pair != -1 && !cannot_add_dangling[i])
            //  )(
        {
            AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
            if (debug)
                printf ("%d - AUpen2 \t- add energy %6g\n", i, AUpen);
            energy += AUpen;
            if (counter != NULL)    count_AU_penalty (sequence[i], sequence[f[i].pair], counter);
        }

        // add dangling energies and AU_penalties
        if (f[i].pair == -1 && !cannot_add_dangling[i])
        {
            if ((i == 0 || i-1 == link || (i > 0 && f[i-1].pair == -1 && i != link)) &&
                 i < nb_nucleotides-1 && f[i+1].pair > i+1 && structure[i+1] != '<')
                // .( or ..(
            {
                if (no_dangling_ends)
                    dang = 0;
                else
                {
                    //dang = MIN (0, IGINF(dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]]));
                    dang = IGINF(dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]]);
                    sprintf (type, "dangle_bot[%d][%d][%d]",sequence[f[i+1].pair], sequence[i+1], sequence[i]);
                    index = structure_type_index (type);
                    if (counter != NULL)    counter[index]++;
                }
                AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
                if (debug)
                {
                    printf ("%d - dangle1 \t- add energy %6g\n", i, dang);
                    printf ("%d - AUpen3 \t- add energy %6g\n", i, AUpen);
                }
                energy += dang + AUpen;
                if (counter != NULL)    count_AU_penalty (sequence[i+1], sequence[f[i+1].pair], counter);
            }
            else if ((i == nb_nucleotides-1 || i == link ||
                     (i < nb_nucleotides-1 && f[i+1].pair == -1 && i-1 != link)) &&
                     i > 0 && f[i-1].pair > -1 && f[i-1].pair < i-1 && structure[i-1] != '>')
                // ). or )..
            {
                if (no_dangling_ends)
                    dang = 0;
                else
                {
                    //dang = MIN (0, IGINF(dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]]));
                    dang = IGINF(dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]]);
                    if (debug)
                        printf ("%d - dangle2 \t- add energy %6g\n", i, dang);
                    energy += dang;
                    sprintf (type, "dangle_top[%d][%d][%d]",sequence[i-1], sequence[f[i-1].pair], sequence[i]);
                    index = structure_type_index (type);
                    if (counter != NULL)    counter[index]++;
                }
            }
            else if (i < nb_nucleotides-1 && f[i+1].pair > i+1 && f[i-1].pair < i-1 && f[i-1].pair != -1
                    && structure[i+1] != '<' && structure[i-1] != '>')
               // ).(
            {
                if (no_dangling_ends)
                    dang = 0;
                else
                {
                    //dang = MIN (0, IGINF(s_dangling_energy (sequence, f[i-1].pair, i-1, i+1, f[i+1].pair)));
                    dang = IGINF(s_dangling_energy (sequence, structure, f[i-1].pair, i-1, i+1, f[i+1].pair));
                    if (counter != NULL)    count_dangling_energy (sequence, structure, link, f[i-1].pair, i-1, i+1, f[i+1].pair, counter);
                }
                AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
                if (debug)
                {
                    printf ("%d - dangle1 \t- add energy %6g\n", i, dang);
                    printf ("%d - AUpen4 \t- add energy %6g\n", i, AUpen);
                }
                energy += dang + AUpen;
                if (counter != NULL)    count_AU_penalty (sequence[i+1], sequence[f[i+1].pair], counter);
            }
            else
            {
                continue;
            }
        }

        if (f[i].pair < i || f[i].type == NONE)
        {
            continue;
        }

        if (f[i].type == STACK)
        {
            en = s_stacked_pair::get_energy (i, f[i].pair, sequence);
            if (debug)
                printf ("%d stack \t- add energy %6g\n", i, en);
            energy += en;
            if (counter != NULL)    s_stacked_pair::count_get_energy (i, f[i].pair, sequence, counter);
        }
        else if (f[i].type == HAIRP)    // means we don't have x's inside
        {
            if (link > -1 && i <= link && link < f[i].pair)
            {
                // add intermolecular initiation
                misc_energy = misc.intermolecular_initiation;
                index = structure_type_index ("misc.intermolecular_initiation");
                if (debug)
                    printf ("%d intermol \t- add energy %6g\n", i, misc.intermolecular_initiation);
                if (counter != NULL)    counter[index]++;
                energy += misc_energy;
                // add AU penalty
                if (counter != NULL)    count_AU_penalty (sequence[i], sequence[f[i].pair], counter);
                AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
                energy += AUpen;
                // add dangling ends
                if (link > i)    // there is a dangle_top    (. )
                {
                    if (no_dangling_ends)   dang = 0;
                    else
                    {
                        //dang = MIN (0, IGINF(dangle_top [sequence[i]] [sequence[f[i].pair]] [sequence[i+1]]));
                        dang = IGINF(dangle_top [sequence[i]] [sequence[f[i].pair]] [sequence[i+1]]);
                        if (debug)
                        {
                            printf ("%d - dangle-spec-hairp \t- add energy %6g\n", i, dang);
                        }
                        energy += dang;
                        sprintf (type, "dangle_top[%d][%d][%d]", sequence[i], sequence[f[i].pair], sequence[i+1]);
                        index = structure_type_index (type);
                        if (counter != NULL)    counter[index]++;
                    }
                }
                if (link < f[i].pair-1)    // there is a dangle_bot    ( .)
                {
                    if (no_dangling_ends)   dang = 0;
                    else
                    {
                        //dang = MIN (0, IGINF(dangle_bot [sequence[i]] [sequence[f[i].pair]] [sequence[f[i].pair-1]]));
                        dang = IGINF(dangle_bot [sequence[i]] [sequence[f[i].pair]] [sequence[f[i].pair-1]]);
                        if (debug)
                        {
                            printf ("%d - dangle-spec-hairp \t- add energy %6g\n", i, dang);
                        }
                        energy += dang;
                        sprintf (type, "dangle_bot[%d][%d][%d]", sequence[i], sequence[f[i].pair], sequence[f[i].pair-1]);
                        index = structure_type_index (type);
                        if (counter != NULL)    counter[index]++;
                    }
                }

            }
            else
            {
                en = s_hairpin_loop::get_energy (i, f[i].pair, sequence, csequence, p_table);
                if (debug)
                    printf ("%d hairpin \t- add energy %6g\n", i, en);
                energy += en;
                if (counter != NULL)    s_hairpin_loop::count_get_energy (i, f[i].pair, sequence, csequence, counter);
            }
        }
        else if (f[i].type == INTER)
        {
            int ip, jp;
            ip = f[i].bri[0];
            jp = f[f[i].bri[0]].pair;
            cannot_add_dangling[ip-1] = 1;
            cannot_add_dangling[jp+1] = 1;
            en = s_internal_loop::get_energy (i, f[i].pair, ip, jp, sequence, p_table);
            if (debug)
                printf ("%d internal \t- add energy %6g\n", i, en);
            energy += en;
            if (counter != NULL)    s_internal_loop::count_get_energy (i, f[i].pair, ip, jp, sequence, counter);
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

            if (f[i].num_branches > MAX_BRANCHES)
            {
                fprintf (stderr, "ERROR! f[i].num_branches=%d, but MAX_BRANCHES=%d\n", f[i].num_branches, MAX_BRANCHES);
                exit(1);
            }
            while (l < f[i].bri[0] && !special)
                if (l++ == link) special = 1;
            h = 0;
            while (h < f[i].num_branches-1 && !special)
            {
                l = f[f[i].bri[h]].pair;
                while (l < f[i].bri[h+1] && !special)
                    if (l++ == link) special = 1;
                h++;
            }
            l = f[f[i].bri[f[i].num_branches-1]].pair;
            while (l < f[i].pair && !special)
                if (l++ == link) special = 1;
            // now we now if it a special multi-loop or not

            if (!special)
            {
//                printf ("Regular ML\n");
              // consider the contribution of unpaired bases
              index = structure_type_index ("misc.multi_free_base_penalty");
                for (l=i+1; l < f[i].bri[0]; l++)
                  {
                    misc_energy += misc.multi_free_base_penalty;
                    if (counter != NULL)    counter[index]++;
                  }
                for (h=0; h < f[i].num_branches-1; h++)
                {
                    for (l = f[f[i].bri[h]].pair + 1; l < f[i].bri[h+1]; l++)
                      {
                        misc_energy += misc.multi_free_base_penalty;
                        if (counter != NULL)    counter[index]++;
                      }
                }
                for (l = f[f[i].bri[f[i].num_branches-1]].pair + 1; l < f[i].pair; l++)
                  {
                    misc_energy += misc.multi_free_base_penalty;
                    if (counter != NULL)    counter[index]++;
                  }
                // done considering the contribution of unpaired bases
                misc_energy += misc.multi_offset;
                misc_energy += misc.multi_helix_penalty * (f[i].num_branches + 1);
                index = structure_type_index ("misc.multi_offset");
                if (counter != NULL)    counter[index]++;
                index = structure_type_index ("misc.multi_helix_penalty");
                if (counter != NULL)    counter[index]+= f[i].num_branches + 1;
            }
            /*
            else
            {
//                printf ("Special ML\n");
                misc_energy = misc.intermolecular_initiation;
                index = structure_type_index ("misc.intermolecular_initiation");
                if (counter != NULL)    counter[index]++;
            }
            */
            // add AU_penalties for multi-loop
            // the closing base pair can't be <>
            AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);
            if (counter != NULL)    count_AU_penalty (sequence[i], sequence[f[i].pair], counter);
            for (h=0; h < f[i].num_branches; h++)
            {
                // ignore if the base pair is <>
                if (structure[f[i].bri[h]] != '<')
                {
                    AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);
                    if (counter != NULL)    count_AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair], counter);
                }
            }

            // add dangling energies for multi-loop
            if (no_dangling_ends)
                dang = 0;
            else
            {
                dang += s_dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
                if (counter != NULL)    count_dangling_energy_left (sequence, structure, link, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair, counter);
                for (l=0; l < f[i].num_branches - 1; l++)
                {
                    dang += s_dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
                    if (counter != NULL)    count_dangling_energy (sequence, structure, link, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair, counter);
                }
                dang += s_dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);
                if (counter != NULL)    count_dangling_energy_right (sequence, structure, link, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair, counter);
            }
            // add "no-dangling" restriction
            for (l=0; l < f[i].num_branches; l++)
            {
                cannot_add_dangling [f[i].bri[l] -1] = 1;
                cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
            }
            if (debug)
            {
                printf ("%d - multi m\t- add energy %6g\n", i, misc_energy);
                printf ("%d - multi d\t- add energy %6g\n", i, dang);
                printf ("%d - multi AU\t- add energy %6g\n", i, AUpen);
            }
            energy += misc_energy + dang + AUpen;
        }
    }
    if (restricted[0] != '\0')  delete [] p_table;
    //printf ("Energy: %d\n", energy);
    //printf ("Exiting count_types\n");
    return energy/100.0;
}


double get_feature_counts (char *sequence, char *structure, char *restricted, double *c, double &f)
// a wrapper around count_each_structure_type
// returns the energy function
{
    return count_each_structure_type (sequence, structure, restricted, c, f, 1);
}


double get_feature_counts_quadratic (char *sequence, char *structure, char *restricted, double **quadratic, double *linear, double &f)
// just a dummy function, equivalent with get_feature_counts
// Used to test create_structural_constraints_simfold.cpp for the quadratic case
{
    return count_each_structure_type (sequence, structure, restricted, linear, f, 1);
}


double get_feature_counts_restricted (char *sequence, char *structure, double *c, double &f, int reset_c, int ignore_dangles, int ignore_first_AU_penalty)
// This function computes the c vector and the f value of a linear energy function c' x + f
// structure contains the following characters:
//      - left and right parentheses mean base pairs, to be considered in the energy model
//      - dots mean unpaired bases, to be considered in the energy model
//      - a sequence of x's means ignore that region from the energy model
//      - left and right angles < > mean they enclose a pseudoknot, and ignore that region from the energy model.
// I'm assuming the following are true about structure:
//  - Regions with x are always inside parentheses or angle brackets: (xxxx) or <xxxx>.
//  - (..xx) or <...> or <...((..))> is fine too, I ignore everything that's inside <>
//      and I ignore everything that's inside a hairpin loop if it has at least one x.
//  - For a case like: (..(...)[[)]], this function will be called 3 times,
//      so structure can never be: (..(...)xx)
//  - An angle pair is always nested within parentheses: you can have this: (<xxx>) but not this (<)>
//  - If an angle pair is inside a multi-loop, consider it is a multi-loop branch
//      Example: (((...(((...)))..<xxxxxxx>.)))
//  - If an angle pair is nested inside another base pair, consider it is a multi-loop with 2 branches.
//      Example: ..(((...<xxxxxxx>..)))
//  - A multi-loop branch of type <xxxx> contributes 1 branch penalty
//
// f is the free value of the energy function c' x + f
// If reset_counter is 1, then counter[i] is set to 0 for all i, 0 <= i < num_params.
//      Otherwise, the new values get added to whatever was in counter.
// If ignore_dangles is 1, then dangling ends are not included anywhere
//  Otherwise, dangling ends are included whereever needed only in the parts containing parentheses and dots.
//      For example, in the case (..<xxx>..), dangling ends are included near the parentheses, but not near the angles.
// If ignore_AU is 1, then AU_penalty is not added at the end of the stem.
//      Otherwise, they are.
// Returns the value of the energy function
// Added on Feb 28, 2008
// Modification on May 26, 2008: if c is NULL, then it doesn't compute the c and f, it only returns the energy
{
    if (ignore_dangles)     no_dangling_ends = 1;
    else                    no_dangling_ends = 0;
    ignore_AU_penalty = ignore_first_AU_penalty;
    double energy = count_each_structure_type (sequence, structure, empty_char, c, f, reset_c);
    return energy;

    // I tried to some check below, but if reset_c is 0, then it won't work.
    /*
    PARAMTYPE params_array[num_params];
    double params_array_double[num_params];
    save_parameters_in_array (params_array);
    for (int i=0; i < num_params; i++)  params_array_double[i] = params_array[i]/100.0;
    if (!check_counts_linear (num_params, params_array_double, c, f, energy))
    {
        fprintf (stderr, "ERROR in simfold's get_feature_counts_restricted!");
        exit(1);
    }
    return energy;
    */
}


int check_counts_linear (int numpars, double *params, double *c, double f, double energy)
// Return 1 if energy ~= c'x + f    (where x is params)
// Return 0 otherwise
// Added on Feb 29, 2008
{
    double energy_c = 0;
    int i;
    for (i=0; i < numpars; i++)
    {
        energy_c +=  c[i] * params[i];
    }
    energy_c += f;

    if (fabs(energy_c - energy) > 0.2)  // not really sure what's a good threshold. Maybe 0 is fine
    {
        printf ("ERROR! Something is wrong with the counts or the free energy: c'x+f = %.2lf, energy = %.2lf diff=%.2lf\n", energy_c, energy, fabs(energy_c-energy));
        return 0;
    }
    return 1;
}


int check_counts_quadratic (int numpars, double *params, double **P, double *c, double f, double energy)
// Return 1 if energy ~= x'Px + c'x + f    (where x is params)
// Return 0 otherwise
// Added on Feb 29, 2008
{
    double energy_c = 0;
    int i, j;
    for (i=0; i < numpars; i++)
    {
        energy_c +=  c[i] * params[i];
    }
    energy_c += f;

    for (i=0; i < numpars; i++)
    {
        for (j=i; j < numpars; j++)
        {
            energy_c += P[i][j]*params[i]*params[j];
        }
    }
    if (fabs(energy_c - energy) > 0.2)  // not really sure what's a good threshold. Maybe 0 is fine
    {
        printf ("ERROR! Something is wrong with the counts or the free energy: x'Px + c'x + f = %g, energy = %g diff=%g\n", energy_c, energy, fabs(energy_c-energy));
        return 0;
    }
    return 1;
}



void print_counter (double *counter, double free_value)
// prints to the screen the counts of each parameter
// Added on Feb 29, 2008
{
    //num_params = create_string_params();
    int i;
    PARAMTYPE params_array[MAXNUMPARAMS];

    save_parameters_in_array (params_array);
    printf ("==================================================================\n");
    printf ("%4s\t%30s\t%10s\t%4s\n", "Par#", "Parameter name", "Value", "Counter");
    printf ("------------------------------------------------------------------\n");
    for (i=0; i < num_params; i++)
    {
        if (counter[i] != 0)
        {
            printf ("%4d\t%30s\t% 10.2lf\t%4.2lf\n", i, string_params_human_readable[i], (double)params_array[i], counter[i]);
        }
    }
    if (free_value != 0)
    {
        printf ("%4d\t%30s\t% 10s\t%4.2lf\n", num_params, "FREE PARAMETER", "", (double)free_value);
    }
    printf ("==================================================================\n");
}


double count_each_structure_type (char *sequence, char *structure, char *restricted, double *counter, double &free_value, int reset)
// Mirela: Nov 23, 2003
// PRE: string_params have been filled, i.e. by calling create_string_params() or create_building_block_strings()
// GIven the sequence and the structure, returns the vector counter, with the #of elementary structures
// The sequence and structure may have a space inside, meaning duplex
// The free_value is f, the free value of the energy function c'x + f
// If reset is 1, then counter will be reset to 0 first, otherwise, the values will be added to the current value of counter
// returns the energy value
{
    int i, j;
    str_features *f;
    int *int_sequence;
    char *actual_seq;
    char *actual_str;
    int len, nb_nucleotides;
    char *space;
    int link;

    if (reset && counter != NULL)
    {
        //first reset counter
        // counter[num_params] used to be the free value, but now free_value is the free value
        for (i=0; i < num_params; i++)
        {
            counter[i] = 0;
        }
        free_value = 0;
    }

    //num_params = create_string_params ();
    len = strlen(sequence);
    // look for the space
    space = strstr (sequence, " ");
    if (space == NULL)    // single sequence
    {
        // make sure the structure doesn't have space either
        if (strstr (structure, " ") != NULL)
        {
            fprintf (stderr, "Structure has space, and sequence not!!\n%s\n%s\n", sequence, structure);
            exit (1);
        }
        actual_seq = sequence;
        actual_str = structure;
        link = -1;
    }
    else
    {
        // make sure structure has space at the same place
        int space_location;
        space_location = space-sequence;
        link = space_location - 1;
        if (structure[space_location] != ' ')
        {
            fprintf (stderr, "Structure doesn't have space at the same place as sequence!!\n%s\n%s\n", sequence, structure);
            exit (1);
        }
        actual_seq = new char[len+1];
        actual_str = new char[len+1];
        for (i=0; i < space_location; i++)
        {
            actual_seq[i] = sequence[i];
            actual_str[i] = structure[i];
        }
        for (i=space_location+1; i < len; i++)
        {
            actual_seq[i-1] = sequence[i];
            actual_str[i-1] = structure[i];
        }
        actual_seq[len-1] = '\0';
        actual_str[len-1] = '\0';
    }

    nb_nucleotides = strlen(actual_seq);
    if ((f = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    // this function call should be the same for single sequence of duplex
    //printf ("Structure:\n%s\n", structure);
    detect_structure_features (actual_str, f);

    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(actual_seq[i]);
    // count_types adds the free value at the end of counter. I don't want to look through all those functions, so I'm just hacking it here

    double energy;
    if (counter == NULL)
    {
        energy = count_types (link, int_sequence, actual_seq, structure, restricted, f, NULL);
    }
    else
    {
        double *counter_and_free_value = new double[num_params+1];
        for (i=0; i <= num_params; i++)  counter_and_free_value[i] = 0;
        energy = count_types (link, int_sequence, actual_seq, structure, restricted, f, counter_and_free_value);
        for (i=0; i < num_params; i++)  counter[i] += counter_and_free_value[i];
        free_value += counter_and_free_value[num_params]/100.0;
        delete [] counter_and_free_value;
    }

    if (link > -1)
    {
        delete [] actual_seq;
        delete [] actual_str;
    }
    delete [] int_sequence;
    delete [] f;
    return energy;
}


int is_int22_group_1 (int i, int j, int k, int l)
// returns 1 if it is in group 1 according to Christiansen_Znosko_2008 and Shankar_Turner_2007
{
    // a U · U pair adjacent to an R · R pair,
    if (i==U && j==U && (k==A || k==G) && (l==A || l==G))   return 1;
    if (k==U && l==U && (i==A || i==G) && (j==A || j==G))   return 1;
    //a G · A  or A · G pair adjacent to a Y · Y pair,
    if (((i==A && j==G) || (i==G && j==A)) && (k==C || k==U) && (l==C || l==U))     return 1;
    if (((k==A && l==G) || (k==G && l==A)) && (i==C || i==U) && (j==C || j==U))     return 1;
    //or  any combination of A · C, U · C, C · U,  C · C, C · A, or A · A pairs
    if (((i==A && j==C) || (i==U && j==C) || (i==C && j==U) || (i==C && j==C) || (i==C && j==A) || (i==A && j==A)) &&
        ((k==A && l==C) || (k==U && l==C) || (k==C && l==U) || (k==C && l==C) || (k==C && l==A) || (k==A && l==A)))
        return 1;
    return 0;
}

int is_int22_group_2 (int i, int j, int k, int l)
// returns 1 if it is in group 2 according to Christiansen_Znosko_2008 and Shankar_Turner_2007
{
    //any combination of adjacent G · A and A · G pairs
    if (((i==A && j==G) || (i==G && j==A)) && ((k==A && l==G) || (k==G && l==A)))   return 1;
    //or two U · U pairs
    if (i==U && j==U && k==U && l==U)   return 1;
    return 0;
}


int is_int22_group_3 (int i, int j, int k, int l)
// According to Shankar_Turner_2007, group 3 is: (Christiansen_Znosko_2008 missed one case!!)
//  GA or AG pair adjacent to a CA, AC, or AA pair, or a UU pair adjacent to a YY or CA or AC pair.
{
    //GA or AG pair adjacent to a CA, AC, or AA pair
    //  - this is from Shankar_Turner_2007, Christiansen_Znosko_2008 missed it
    if (((i==G && j==A) || (i==A && j==G)) && ((k==A && l==C) || (k==C && l==A) || (k==A && l==A))) return 1;
    // now the mirrored
    if (((k==G && l==A) || (k==A && l==G)) && ((i==A && j==C) || (i==C && j==A) || (i==A && j==A))) return 1;
    //a U · U pair adjacent to a Y · Y (not U · U),  C · A, or A · C pair
    if (i==U && j==U && k==U && l==U)   return 0;
    if (i==U && j==U && ((k==C || k==U) && (l==C || l==U)))     return 1;
    if (i==U && j==U && ((k==C && l==A) || (k==A && l==C)))     return 1;
    // now the mirrored
    if (k==U && l==U && ((i==C || i==U) && (k==C || k==U)))     return 1;
    if (k==U && l==U && ((i==C && j==A) || (i==A && j==C)))     return 1;
    return 0;
}

int is_int22_group_4 (int i, int j, int k, int l)
// returns 1 if it is in group 4 according to Christiansen_Znosko_2008 and Shankar_Turner_2007
{
    //a G · G pair not adjacent to a U · U pair
    if (i==G && j==G)
    {
        if (k==U && l==U)   return 0;
        return 1;
    }
    if (k==G && l==G)
    {
        if (i==U && j==U)   return 0;
        return 1;
    }
    return 0;
}


int apply_rule_1 (int i, int j, int &i_rule1, int &j_rule1)
// I Don't use this when I initialize the variable any more, I just used it for the feature similarity rules
// check if i-j form a base pair: A-U, C-G or G-U.
// If they do, return 1, and write the replacement in i_rule1 and j_rule1
// If they don't, return 0, and i_rule1=i, j_rule1=j.
        // There's a problem in the case of int21: If the free bases are UG/G, after applying this rule we get CG/A, but G/G "pairs" have some contribution in int21
{
    i_rule1 = i;
    j_rule1 = j;
    if (i==A && j==U)
    {
        j_rule1=C;  return 1;
    }
    else if (i==U && j==A)
    {
        i_rule1=C;  return 1;
    }
    else if (i==G && j==C)
    {
        i_rule1=A;  return 1;
    }
    else if (i==C && j==G)
    {
        j_rule1=A;  return 1;
    }
    else if (i==G && j==U)
    {
        i_rule1=A; j_rule1=C;  return 1;
    }
    else if (i==U && j==G)
    {
        i_rule1=C; j_rule1=A;  return 1;
    }
    return 0;
}

void check_int11_parameters (int i, int j, int k, int l, int m, int n)
        // check if int11 is the sum up of experimental addition etc.
{
    if (parsi_int11 != LAVISH && parsi_int11 != HLI)  return;
    PARAMTYPE int11_shouldbe;
    if (int11_experimental_addition[i][j][k][l][m][n] < INF)
    {
        int11_shouldbe = int11_experimental_addition[i][j][k][l][m][n];
        if ((i==A && j==U) || (i==U && j==A))
            int11_shouldbe += misc.internal11_AU_closure;
        if ((m==A && n==U) || (m==U && n==A))
            int11_shouldbe += misc.internal11_AU_closure;
        // look for GU closure
        if ((i==G && j==U) || (i==U && j==G))
            int11_shouldbe += misc.internal11_GU_closure;
        if ((m==G && n==U) || (m==U && n==G))
            int11_shouldbe += misc.internal11_GU_closure;
        // look for AG mismatch
        if ((k==A && l==G) || (k==G && l==A))
            int11_shouldbe += misc.internal11_AG_mismatch;
        // look for GG mismatch
        if (k==G && l==G)
            int11_shouldbe += misc.internal11_GG_mismatch;
        // look for UU mismatch
        if (k==U && l==U)
            int11_shouldbe += misc.internal11_UU_mismatch;
        // check if it is internal11_5YRR_5YRR
        if (isY(i) && isR(j) && isR(k) && isR(l) && isR(m) && isY(n))
            int11_shouldbe += misc.internal11_5YRR_5YRR;
        if ( isR(i) && isY(j) && isY(k) && isY(l) && isY(m) && isR(n) )
            int11_shouldbe += misc.internal11_5RYY_5RYY;
        if ( isY(i) && isR(j) && isY(k) && isY(l) && isR(m) && isY(n) )
            int11_shouldbe += misc.internal11_5YYR_5YYR;
        if ( (isY(i) && isR(j) && isR(k) && isY(l) && isY(m) && isR(n)) ||
              (isR(i) && isY(j) && isY(k) && isR(l) && isR(m) && isY(n)) )
            int11_shouldbe += misc.internal11_5YRY_5RYR;
        if ( (isR(i) && isY(j) && isR(k) && isY(l) && isY(m) && isR(n)) ||
              (isR(i) && isY(j) && isY(k) && isR(l) && isY(m) && isR(n)) )
            int11_shouldbe += misc.internal11_5RRY_5RYY;

        if (fabs(int11_shouldbe - int11[i][j][k][l][m][n]) > 0.01)
        {
            fprintf (stderr, "DIFFERENCE between what int11 should be (%Lg) and what it is (%Lg) at 5'-%c%c%c/%c%c%c-3\n",
                    int11_shouldbe,  int11[i][j][k][l][m][n], int_to_nuc(i), int_to_nuc(k), int_to_nuc(m),
                            int_to_nuc(n), int_to_nuc(l), int_to_nuc(j));
            exit(1);
        }
    }
}

void similarity_bulge_type (int bulged_base, char *similarity)
//  computes the similarity tule
// TODO
{
    int i, j, ip, jp, tindex;
    char type[100];
    int num = 0;
    int num2 = 0;
    // first count how many they are
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        {
            if (! can_pair (i,j))   continue;
            for (ip=0; ip < NUCL; ip++)
                for (jp=0; jp < NUCL; jp++)
                {
                    if (! can_pair(ip, jp)) continue;
                    sprintf (type, "bulge1[%d][%d][%d][%d][%d]", i, j, bulged_base, ip, jp);
                    tindex = structure_type_index (type);
                    if (similarity_rule[tindex][0] != '\0')
                    {
                        num++;
                    }
                }
        }
    // now traverse again
    // TODO: this is veeery slow
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        {
            if (! can_pair (i,j))   continue;
            for (ip=0; ip < NUCL; ip++)
                for (jp=0; jp < NUCL; jp++)
                {
                    if (! can_pair(ip, jp)) continue;
                    sprintf (type, "bulge1[%d][%d][%d][%d][%d]", i, j, bulged_base, ip, jp);
                    tindex = structure_type_index (type);
                    if (similarity_rule[tindex][0] != '\0')
                    {
                        sprintf (similarity, "%s%.4lf * bulge1[5'-%c%c%c/%c%c-3'] - %.4lf * stack[5'-%c%c/%c%c-3']", similarity, 1.0/num,
                            int_to_nuc(i), int_to_nuc(bulged_base), int_to_nuc(ip), int_to_nuc(jp), int_to_nuc(j),
                            1.0/num,  int_to_nuc(i), int_to_nuc(ip), int_to_nuc(jp), int_to_nuc(j));
                        num2++;
                        if (num2 < num)     sprintf (similarity, "%s + ", similarity);
                    }
                }
        }
}



void extrapolate_parameters ()
// Start from the basic set of parameters, according to the model complexity,
//      and fill up the remaining structures.
{
    int index;
    int start;  // used in penalty_by_size
    int i, j, k, l, m, n, o, p, ip, jp;
    index = 0;
    char type[100];
    int tindex;
    int kalt;   // alternative to k
    double param;

    //#if (MODEL == SIMPLE)
    // set a fixed value to param_greater_30 for now
    misc.param_greater30 = 1.079;
    // keep misc.gail_rule this fixed
    misc.gail_rule = 1;

    if (parsi_int11 == T99)
    {
        // fill the int11 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        if ( ((i==C && j==G) || (i==G && j==C)) &&
                                             ((m==C && n==G) || (m==G && n==C)))
                                        {
                                            if (can_pair(k,l))
                                            {
                                                int11[i][j][k][l][m][n] = misc.internal11_basic_mismatch;
                                                // add the duplicate too
                                                //int11[n][m][l][k][j][i] = misc.internal11_basic_mismatch;
                                            }
                                        }

                                        else
                                        {
                                            if (!(watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                            {
                                                if (k==G && l==G)
                                                {
                                                    int11[i][j][k][l][m][n] = misc.internal11_GG_mismatch;
                                                    // add the duplicate too
                                                    //int11[n][m][l][k][j][i] = misc.internal11_GG_mismatch;
                                                }
                                                else
                                                    int11[i][j][k][l][m][n] = misc.internal11_basic_mismatch;
                                                if (has_AU_penalty(i,j))
                                                    int11[i][j][k][l][m][n] += misc.internal_AU_closure;
                                                if (has_AU_penalty(m,n))
                                                    int11[i][j][k][l][m][n] += misc.internal_AU_closure;
                                            }
                                        }

                                        // round it to match Turner parameters
                                        //if (int11[i][j][k][l][m][n] % 10 == 5) int11[i][j][k][l][m][n] += 5;
                                        //if (int11[i][j][k][l][m][n] % 10 == -5) int11[i][j][k][l][m][n] += 5;
                                    }
                }
    }
    if (parsi_int21 == T99)
    {
        // fill the int21 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        for(o=0; o < NUCL; o++)
                                        {
                                            if ((i==C && j==G && m==C && n==G) ||  // these are already filled above, except what can pair inside
                                                (i==G && j==C && m==G && n==C))
                                            {
                                                if (can_pair(k,l) || can_pair(k,o))
                                                    int21[i][j][k][l][m][n][o] = misc.internal21_match;
                                            }
                                            else
                                            {
                                                if (can_pair(k,l) || can_pair(k,o))
                                                    int21[i][j][k][l][m][n][o] = misc.internal21_match;
                                                // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                                                else
                                                    int21[i][j][k][l][m][n][o] = (PARAMTYPE)(int21[C][G][k][l][C][G][o]/2.0) +
                                                    (PARAMTYPE)(int21[G][C][k][l][G][C][o]/2.0);
                                                if (has_AU_penalty(i,j))
                                                    int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                                if (has_AU_penalty(m,n))
                                                    int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                            }
                                            // round it to match Turner parameters - seems to be inconsistent
                                            //if (int21[i][j][k][l][m][n][o] % 10 == 5) int21[i][j][k][l][m][n][o] += 5;
                                        }
                                    }
                }
    }
    if (parsi_int22 == T99)
    {
        int ii, jj, mm, nn;
        // fill the int22 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        for(o=0; o < NUCL; o++)
                                            for (p=0; p < NUCL; p++)
                                            {
                                                /*
                                                if (i==C && j==G && m==C && n==G)
                                                {
                                                    if(watson_crick(k,l) || watson_crick(o,p))
                                                    {
                                                        int22[i][j][k][l][m][n][o][p] = misc.internal22_match;
                                                    }
                                                    // else do nothing, it's parameter
                                                } */
                                                // if a closing pair is wobble, it's the same as if G would be A
                                                if (i==G && j==U)   ii = A;     else ii = i;
                                                if (i==U && j==G)   jj = A;     else jj = j;
                                                if (m==G && n==U)   mm = A;     else mm = m;
                                                if (m==U && n==G)   nn = A;     else nn = n;


                                                if (watson_crick(k,l) || watson_crick(o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] = misc.internal22_match;
                                                }
                                                else if (nn==ii && mm==jj && p==k && o==l)  // the UG closing pairs are the same as UA
                                                {
                                                    int22[i][j][k][l][m][n][o][p] = int22[ii][jj][k][l][mm][nn][o][p];
                                                }
                                                else //if (!(n==i && m==j && p==k && o==l))   // was already filled above
                                                {
                                                    int result = check_stability_and_size (k, l, o, p);
                                                    // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                                                    PARAMTYPE temp = (PARAMTYPE)(int22[ii][jj][k][l][jj][ii][l][k]/2.0) +
                                                               (PARAMTYPE)(int22[nn][mm][p][o][mm][nn][o][p]/2.0);
                                                    // rounf it to match Turner parameters
                                                    //if (temp%10 == 5) temp -= 5; if (temp%10 == -5) temp += 5;
                                                    int22[i][j][k][l][m][n][o][p] =  temp;
                                                    switch (result)
                                                    {
                                                        case 1: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_same_size; break;
                                                        case 2: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_different_size; break;
                                                        case 3: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_1stable_1unstable; break;
                                                        case 4: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_AC; break;
                                                        default: fprintf (stderr, "ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, k, l, o, p); exit(1);
                                                    }
                                                }
                                            }
                                    }
                }
    }

    if (parsi_tstacki == T99)
    {
        // fill the tstacki data structure, now that we also have terminal_AU_penalty - actually I removed it from here
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        if (!can_pair (i, j))
                            tstacki[i][j][k][l] = INF;
                        else
                        {
                            tstacki[i][j][k][l] = 0;
                            if (((i == A || i == G) && j == U) ||
                                ((j == A || j == G) && i == U))
                            {
                                tstacki[i][j][k][l] += misc.internal_AU_closure;
                            }
                            if ((k == A && l == G) ||
                                (l == A && k == G))
                            {
                                tstacki[i][j][k][l] += misc.internal_GA_AG_mismatch;
                            }
                            if (k == U && l == U)
                            {
                                tstacki[i][j][k][l] += misc.internal_UU_mismatch;
                            }
                        }
                    }
    }

    //#elif (MODEL == EXTENDED)
    int rep1, rep2;
    if (parsi_tstackh == PARSI)      // fill up tstackh from 4 parameters
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        if (!can_pair (i,j))    tstackh[i][j][k][l] = INF;
                        else
                        {
                            // exclude duplicates
                            // stack[i][j][k][l] is the same as stack[l][k][j][i]
                            tstackh[i][j][k][l] = 0;
                            if (((i==A || i==G) && j==U) || (i==U && (j==A || j==G)))
                                tstackh[i][j][k][l] += misc.hairpin_AU_closure;
                            if (k == A && l == G)
                                tstackh[i][j][k][l] += misc.hairpin_AG_mismatch;
                            if (k == G && l == A)
                                tstackh[i][j][k][l] += misc.hairpin_GA_mismatch;
                            if (k == U && l == U)
                                tstackh[i][j][k][l] += misc.hairpin_UU_mismatch;
                        }
                    }
        // apply rule 1 -- not so sure I should, it wasn't done in the basic model
//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//                 for (k=0; k < NUCL; k++)
//                     for (l=0; l < NUCL; l++)
//                     {
//                         if (apply_rule_1 (k, l, rep1, rep2))
//                             tstackh[i][j][k][l] = tstackh[i][j][rep1][rep2];
//                     }
    }
    if (parsi_tstacki == PARSI)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        if (!can_pair (i, j))   tstacki[i][j][k][l] = INF;
                        else
                        {
                            tstacki[i][j][k][l] = 0;
                            if (((i == A || i == G) && j == U) || ((j == A || j == G) && i == U))
                                tstacki[i][j][k][l] += misc.internal_AU_closure;
                            if (k == A && l == G)
                                tstacki[i][j][k][l] += misc.internal_AG_mismatch;
                            if (k == G && l == A)
                                tstacki[i][j][k][l] += misc.internal_GA_mismatch;
                            if (k == G && l == G)
                                tstacki[i][j][k][l] += misc.internal_GG_mismatch;
                            if (k == U && l == U)
                                tstacki[i][j][k][l] += misc.internal_UU_mismatch;
                        }
                    }
        // apply rule 1 -- not so sure I should, it wasn't done in the basic model
//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//                 for (k=0; k < NUCL; k++)
//                     for (l=0; l < NUCL; l++)
//                     {
//                         if (apply_rule_1 (k, l, rep1, rep2))
//                             tstacki[i][j][k][l] = tstacki[i][j][rep1][rep2];
//                     }
    }
    ///////////// internal loop asymmetry
    if (parsi_asymmetry == PARSI)
    {
        for (i=1; i < MAXLOOP_ASYM; i++)
        {
            internal_asymmetry[i] = (PARAMTYPE) (internal_asymmetry_initiation + internal_asymmetry_slope * log(i));
        }
    }
    ////////////// int 11
    if (parsi_int11 == PARSI || parsi_int11 == LAVISH || parsi_int11 == HLI)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {
                                if (!can_pair(i,j))         int11[i][j][k][l][m][n] = INF;
                                else if (!can_pair(m,n))    int11[i][j][k][l][m][n] = INF;
                                else
                                {
                                    if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                    {
                                        // do nothing if it's !parsi_int11 and it doesn't have experimental support
                                        if (parsi_int11 == LAVISH && int11_experimental_addition[i][j][k][l][m][n] >= INF)   continue;
                                        // do nothing if parsi_int11 is HLI and there is experimental support
                                        if (parsi_int11 == HLI && int11_experimental_addition[i][j][k][l][m][n] < INF)   continue;
                                        int11[i][j][k][l][m][n] = 0;
                                        // look for AU closure
                                        if ((i==A && j==U) || (i==U && j==A))
                                            int11[i][j][k][l][m][n] += misc.internal11_AU_closure;
                                        if ((m==A && n==U) || (m==U && n==A))
                                            int11[i][j][k][l][m][n] += misc.internal11_AU_closure;
                                        // look for GU closure
                                        if ((i==G && j==U) || (i==U && j==G))
                                            int11[i][j][k][l][m][n] += misc.internal11_GU_closure;
                                        if ((m==G && n==U) || (m==U && n==G))
                                            int11[i][j][k][l][m][n] += misc.internal11_GU_closure;
                                        // look for AG mismatch
                                        if ((k==A && l==G) || (k==G && l==A))
                                            int11[i][j][k][l][m][n] += misc.internal11_AG_mismatch;
                                        // look for GG mismatch
                                        if (k==G && l==G)
                                            int11[i][j][k][l][m][n] += misc.internal11_GG_mismatch;
                                        // look for UU mismatch
                                        if (k==U && l==U)
                                            int11[i][j][k][l][m][n] += misc.internal11_UU_mismatch;
                                        // check the nearest neighbours
                                        // this one is symmetric
                                        if (isY(i) && isR(j) && isR(k) && isR(l) && isR(m) && isY(n))
                                            int11[i][j][k][l][m][n] += misc.internal11_5YRR_5YRR;
                                        // this one is symmetric
                                        if ( isR(i) && isY(j) && isY(k) && isY(l) && isY(m) && isR(n) )
                                            int11[i][j][k][l][m][n] += misc.internal11_5RYY_5RYY;
                                        // this one is symmetric
                                        if ( isY(i) && isR(j) && isY(k) && isY(l) && isR(m) && isY(n) )
                                            int11[i][j][k][l][m][n] += misc.internal11_5YYR_5YYR;
                                        // this one is NOT symmetric, so we have to consider it in both directions
                                        if ( (isY(i) && isR(j) && isR(k) && isY(l) && isY(m) && isR(n)) ||
                                                (isR(i) && isY(j) && isY(k) && isR(l) && isR(m) && isY(n)) )
                                            int11[i][j][k][l][m][n] += misc.internal11_5YRY_5RYR;
                                        // this one is NOT symmetric, so we have to consider it in both directions
                                        if ( (isR(i) && isY(j) && isR(k) && isY(l) && isY(m) && isR(n)) ||
                                                (isR(i) && isY(j) && isY(k) && isR(l) && isY(m) && isR(n)) )
                                            int11[i][j][k][l][m][n] += misc.internal11_5RRY_5RYY;
                                        // Now, if !parsi_int11, add the experimental addition
                                        if (parsi_int11 == LAVISH || parsi_int11 == HLI)
                                            int11[i][j][k][l][m][n] += int11_experimental_addition[i][j][k][l][m][n];
                                    }
                                }
                            }
            // now fill up the mirrored values
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {
                                if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n > n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                {
                                    int11[i][j][k][l][m][n] = int11[n][m][l][k][j][i];
                                }
                            }
        // apply rule 1 -- not so sure I should, it wasn't done in the basic model
    //     for (i=0; i < NUCL; i++)
    //         for (j=0; j < NUCL; j++)
    //             for (k=0; k < NUCL; k++)
    //                 for (l=0; l < NUCL; l++)
    //                     for (m=0; m < NUCL; m++)
    //                         for (n=0; n < NUCL; n++)
    //                         {
    //                             if (apply_rule_1 (k, l, rep1, rep2))
    //                             {
    //                                 int11[i][j][k][l][m][n] = int11[i][j][rep1][rep2][m][n];
    //                             }
    //                         }
    }
    ///////////////// int 21
    // it starts the same for parsi_int21 and !parsi_int21
    // I'm also applying rule 1 in this case!
    // this is essentially the same as what we have in the similarity rules
    if (parsi_int21 == PARSI || parsi_int21 == LAVISH)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for (o=0; o < NUCL; o++)
                                {
                                    if (!can_pair(i,j))         int21[i][j][k][l][m][n][o] = INF;
                                    else if (!can_pair(m,n))    int21[i][j][k][l][m][n][o] = INF;
                                    else
                                    {
                                        if (parsi_int21 == LAVISH && int21_experimental_addition[i][j][k][l][m][n][o] >= INF)   continue;
                                        int21[i][j][k][l][m][n][o] = misc.internal21_initiation;

                                        // look for AU closure
                                        if ((i==A && j==U) || (i==U && j==A))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                        if ((m==A && n==U) || (m==U && n==A))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                        // look for GU closure
                                        if ((i==G && j==U) || (i==U && j==G))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_GU_closure;
                                        if ((m==G && n==U) || (m==U && n==G))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_GU_closure;
                                        // look for AG mismatch - but not applied to 5'RA/3'YG loops
                                        if ((k==A && l==G &&   i!=A && i!=G   &&   j!=U && j!=C) ||
                                            (k==G && l==A) ||
                                            (k==G && o==A &&   m!=U && m!=C   &&   n!=A && n!=G) ||
                                            (k==A && o==G))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_AG_mismatch;
                                        // look for GG mismatch
                                        if (k==G && (l==G || o==G))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_GG_mismatch;
                                        // look for UU mismatch
                                        if (k==U && (l==U || o==U))
                                            int21[i][j][k][l][m][n][o] += misc.internal21_UU_mismatch;
                                        if (!parsi_int21)
                                            int21[i][j][k][l][m][n][o] += int21_experimental_addition[i][j][k][l][m][n][o];
                                    }
                                }
        // apply rule 1 -- not so sure I should, it wasn't done in the basic model

    //     for (i=0; i < NUCL; i++)
    //         for (j=0; j < NUCL; j++)
    //             for (k=0; k < NUCL; k++)
    //                 for (l=0; l < NUCL; l++)
    //                     for (m=0; m < NUCL; m++)
    //                         for (n=0; n < NUCL; n++)
    //                             for (o=0; o < NUCL; o++)
    //                             {
    //                                 int rep3, rep4;
    //                                 apply_rule_1 (k, l, rep1, rep2);
    //                                 apply_rule_1 (rep1, o, rep3, rep4);
    //                                 int21[i][j][k][l][m][n][o] = int21[i][j][rep3][rep2][m][n][rep4];
    //                             }
    }

    ////////////// int 22
    // this is for both parsi_int22 and !parsi_int22
    if (parsi_int22 == PARSI || parsi_int22 == LAVISH)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for (o=0; o < NUCL; o++)
                                    for (p=0; p < NUCL; p++)
                                    {
                                        if (!can_pair(i,j))         int22[i][j][k][l][m][n][o][p] = INF;
                                        else if (!can_pair(m,n))    int22[i][j][k][l][m][n][o][p] = INF;
                                        else
                                        {
                                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                                n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                            {
                                                if (parsi_int22 == LAVISH && int22_experimental_addition[i][j][k][l][m][n][o][p] >= INF)   continue;
                                                int22[i][j][k][l][m][n][o][p] = 0;
                                                if (is_int22_group_1 (k,l,o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22mid_group1;
                                                }
                                                else if (is_int22_group_2 (k,l,o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22mid_group2;
                                                }
                                                else if (is_int22_group_3 (k,l,o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22mid_group3;
                                                }
                                                else if (is_int22_group_4 (k,l,o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22mid_group4;
                                                }
                                                if ((i==A && j==U) || (i==U && j==A))
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22_AU_closure;
                                                if ((m==A && n==U) || (m==U && n==A))
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22_AU_closure;
                                                // look for GU closure
                                                if ((i==G && j==U) || (i==U && j==G))
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22_GU_closure;
                                                if ((m==G && n==U) || (m==U && n==G))
                                                    int22[i][j][k][l][m][n][o][p] += misc.internal22_GU_closure;
                                                if (parsi_int22 == LAVISH)
                                                    int22[i][j][k][l][m][n][o][p] += int22_experimental_addition[i][j][k][l][m][n][o][p];
                                            }
                                        }
                                    }
        // now fill up the mirrored values
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for (o=0; o < NUCL; o++)
                                    for (p=0; p < NUCL; p++)
                                    {
                                        if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p >
                                            n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                        {
                                            int22[i][j][k][l][m][n][o][p] = int22[n][m][p][o][j][i][l][k];
                                        }
                                    }
        // apply rule 1 -- not so sure I should, it wasn't done in the basic model

    //     for (i=0; i < NUCL; i++)
    //         for (j=0; j < NUCL; j++)
    //             for (k=0; k < NUCL; k++)
    //                 for (l=0; l < NUCL; l++)
    //                     for (m=0; m < NUCL; m++)
    //                         for (n=0; n < NUCL; n++)
    //                             for (o=0; o < NUCL; o++)
    //                                 for (p=0; p < NUCL; p++)
    //                                 {
    //                                     if (apply_rule_1 (k, l, rep1, rep2))
    //                                         int22[i][j][k][l][m][n][o][p] = int22[i][j][rep1][rep2][m][n][o][p];
    //                                     if (apply_rule_1 (o, p, rep1, rep2))
    //                                         int22[i][j][k][l][m][n][o][p] = int22[i][j][k][l][m][n][rep1][rep2];
    //                                 }
    }

    ////////////// bulge1
    if (parsi_bulge1 == PARSI)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (ip=0; ip < NUCL; ip++)
                        for (jp=0; jp < NUCL; jp++)
                        {
                            if (!can_pair (i,j))            bulge1[i][j][k][ip][jp] = INF;
                            else if (!can_pair (ip,jp))     bulge1[i][j][k][ip][jp] = INF;
                            else
                            {
                                PARAMTYPE bulged;
                                if (k==A)        bulged = bulgeA;
                                else if (k==C)   bulged = bulgeC;
                                else if (k==G)   bulged = bulgeG;
                                else             bulged = bulgeU;
                                bulge1[i][j][k][ip][jp] = stack[i][j][ip][jp] + bulged;
                            }
                        }
    }
    ////////////// dangling ends
    if (parsi_dangles == PARSI)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                {
                    if (!can_pair(i,j))
                    {
                        dangle_top[i][j][k] = INF;
                        dangle_bot[i][j][k] = INF;
                    }
                    else
                    {
                        dangle_top[i][j][k] = 0;
                        dangle_bot[i][j][k] = 0;
                    }
                }
    }
    //if (parsi_others)   // I have nothing in here for now

    ///////////// length of internal, hairpin and bulge loops
    if (parsi_length == PARSI)
    {
        double logval;
        misc.param_greater30 = 1.079;
        // 1.079 comes from: T * 1.75 * R = 310.15 * 1.75 * 1.98717 / 1000
        for (i=MAXLOOP_H_PARSI+1; i <= MAXLOOP_H_LAVISH; i++)
        {
            logval = log (1.0*i/MAXLOOP_H_PARSI);
            hairpin_penalty_by_size[i] = hairpin_penalty_by_size[MAXLOOP_H_PARSI] + (PARAMTYPE)(100.0*misc.param_greater30 * logval);
        }
        for (i=MAXLOOP_I_PARSI+1; i <= MAXLOOP_I_LAVISH; i++)
        {
            logval = log (1.0*i/MAXLOOP_I_PARSI);
            internal_penalty_by_size[i] = internal_penalty_by_size[MAXLOOP_I_PARSI] + (PARAMTYPE)(100.0*misc.param_greater30 * logval);
        }
        for (i=MAXLOOP_B_PARSI+1; i <= MAXLOOP_B_LAVISH; i++)
        {
            logval = log (1.0*i/MAXLOOP_B_PARSI);
            bulge_penalty_by_size[i] = bulge_penalty_by_size[MAXLOOP_B_PARSI] + (PARAMTYPE)(100.0*misc.param_greater30 * logval);
        }
    }
    else if (parsi_length == ZL)
    {
        double logval;
        double Tem = 310.15;
        double R = 1.98717/1000.0;
        double multiplier;      // according to Zhang_Liang_2008
        for (i=MAXLOOP_H_PARSI+1; i <= MAXLOOP_H_LAVISH; i++)
        {
            multiplier = 1.75;
            logval = log (1.0*i/MAXLOOP_H_PARSI);
            hairpin_penalty_by_size[i] = hairpin_penalty_by_size[MAXLOOP_H_PARSI] +
                (PARAMTYPE)(100.0 * R*Tem*multiplier * logval);
        }
        for (i=MAXLOOP_I_PARSI+1; i <= MAXLOOP_I_LAVISH; i++)
        {
            multiplier = 1.55;
            logval = log (1.0*i/MAXLOOP_I_PARSI);
            internal_penalty_by_size[i] = internal_penalty_by_size[MAXLOOP_I_PARSI] +
                (PARAMTYPE)(100.0* R*Tem*multiplier * logval);
        }
        for (i=MAXLOOP_B_PARSI+1; i <= MAXLOOP_B_LAVISH; i++)
        {
            multiplier = 1.85;
            logval = log (1.0*i/MAXLOOP_B_PARSI);
            bulge_penalty_by_size[i] = bulge_penalty_by_size[MAXLOOP_B_PARSI] +
                (PARAMTYPE)(100.0* R*Tem*multiplier * logval);
        }
    }


    /////////// special hairpin loops and internal loops
    if (parsi_special == PARSI)
    {
        misc.hairpin_GGG = 0;
        misc.hairpin_c1 = 0;
        misc.hairpin_c2 = 0;
        misc.hairpin_c3 = 0;
        nb_special_hl = 0;
    }

}

void initialize_correct_int11_expadd (int ii, int jj, int kk, int ll, int mm, int nn)
// assume int11 is already initialized, which it should, since init_data should already be called
// Initialize int11_experimental_addition with the correct value (i.e. such that int11 is the sum of the other ones)
// this is similar to the function count_int11_MODEL_EXTENDED in s_internal_loop.cpp
{
    if (parsi_int11 != LAVISH && parsi_int11 != HLI)  return;
    // DO NOT apply rule 1, otherwise it skrews up the string_params array
    //apply_rule_1 (kk, ll, kk, ll);

    int11_experimental_addition[ii][jj][kk][ll][mm][nn] = int11[ii][jj][kk][ll][mm][nn];
    if ((ii==A && jj==U) || (ii==U && jj==A))
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_AU_closure;
    if ((mm==A && nn==U) || (mm==U && nn==A))
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_AU_closure;
    // look for GU closure
    if ((ii==G && jj==U) || (ii==U && jj==G))
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_GU_closure;
    if ((mm==G && nn==U) || (mm==U && nn==G))
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_GU_closure;
    // look for AG mismatch
    if ((kk==A && ll==G) || (kk==G && ll==A))
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_AG_mismatch;
    // look for GG mismatch
    if (kk==G && ll==G)
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_GG_mismatch;
    // look for UU mismatch
    if (kk==U && ll==U)
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_UU_mismatch;
    // check if it is internal11_5YRR_5YRR
    if (isY(ii) && isR(jj) && isR(kk) && isR(ll) && isR(mm) && isY(nn))
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_5YRR_5YRR;
    if ( isR(ii) && isY(jj) && isY(kk) && isY(ll) && isY(mm) && isR(nn) )
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_5RYY_5RYY;
    if ( isY(ii) && isR(jj) && isY(kk) && isY(ll) && isR(mm) && isY(nn) )
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_5YYR_5YYR;
    if ( (isY(ii) && isR(jj) && isR(kk) && isY(ll) && isY(mm) && isR(nn)) ||
        (isR(ii) && isY(jj) && isY(kk) && isR(ll) && isR(mm) && isY(nn)) )
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_5YRY_5RYR;
    if ( (isR(ii) && isY(jj) && isR(kk) && isY(ll) && isY(mm) && isR(nn)) ||
        (isR(ii) && isY(jj) && isY(kk) && isR(ll) && isY(mm) && isR(nn)) )
        int11_experimental_addition[ii][jj][kk][ll][mm][nn] -= misc.internal11_5RRY_5RYY;
}


void initialize_correct_int21_expadd (int ii, int jj, int kk, int ll, int mm, int nn, int oo)
// assume int21 is already initialized, which it should, since init_data should already be called
// Initialize int21_experimental_addition with the correct value (i.e. such that int21 is the sum of the other ones)
// this is similar to the function count_int21_MODEL_EXTENDED in s_internal_loop.cpp
{
    if (parsi_int21 != LAVISH)  return;
    // DO NOT apply rule 1, otherwise it skrews up the string_params array
    //apply_rule_1 (kk, ll, kk, ll);
    //apply_rule_1 (kk, oo, kk, oo);

    int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] = int21[ii][jj][kk][ll][mm][nn][oo];
    int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_initiation;

    // look for AU closure
    if ((ii==A && jj==U) || (ii==U && jj==A))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_AU_closure;
    if ((mm==A && nn==U) || (mm==U && nn==A))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_AU_closure;
    // look for GU closure
    if ((ii==G && jj==U) || (ii==U && jj==G))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_GU_closure;
    if ((mm==G && nn==U) || (mm==U && nn==G))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_GU_closure;
    // look for AG mismatch - but not applied to 5'RA/3'YG loops
    if ((kk==A && ll==G &&   ii!=A && ii!=G   &&   jj!=U && jj!=C) ||
        (kk==G && ll==A) ||
        (kk==G && oo==A &&   mm!=U && mm!=C   &&   nn!=A && nn!=G) ||
        (kk==A && oo==G))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_AG_mismatch;
    // look for GG mismatch
    if (kk==G && (ll==G || oo==G))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_GG_mismatch;
    // look for UU mismatch
    if (kk==U && (ll==U || oo==U))
        int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] -= misc.internal21_UU_mismatch;

}


void initialize_correct_int22_expadd (int ii, int jj, int kk, int ll, int mm, int nn, int oo, int pp)
// assume int22 is already initialized, which it should, since init_data should already be called
// Initialize int22_experimental_addition with the correct value (i.e. such that int22 is the sum of the other ones)
// this is similar to the function count_int22_MODEL_EXTENDED in s_internal_loop.cpp
{
    if (parsi_int22 != LAVISH)  return;

    // Apply rule 1 twice
    // DO NOT apply rule 1, otherwise it skrews up the string_params array
    //apply_rule_1 (kk, ll, kk, ll);
    //apply_rule_1 (oo, pp, oo, pp);

    // Applying rule 1 might not get the order of ii, jj, kk, ll, mm, nn, oo, pp to be in the first symmetric part, which is part of the feature set

    int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] = int22[ii][jj][kk][ll][mm][nn][oo][pp];

    if (is_int22_group_1 (kk, ll, oo, pp))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22mid_group1;
    else if (is_int22_group_2 (kk, ll, oo, pp))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22mid_group2;
    else if (is_int22_group_3 (kk, ll, oo, pp))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22mid_group3;
    else if (is_int22_group_4 (kk, ll, oo, pp))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22mid_group4;
    if ((ii == A && jj == U) || (ii == U && jj == A))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22_AU_closure;
    else if ((ii == G && jj == U) || (ii == U && jj == G))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22_GU_closure;
    if ((mm == A && nn == U) || (mm == U && nn == A))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22_AU_closure;
    else if ((mm == G && nn == U) || (mm == U && nn == G))
        int22_experimental_addition[ii][jj][kk][ll][mm][nn][oo][pp] -= misc.internal22_GU_closure;

}


int traverse_features_and_do_work (const char *calling_function, PARAMTYPE *array, const char *filename)
// This function should be called by create_string_params and other functions, with the name of the calling function as argument
// The purpose of it is to traverse the model's features in only one function instead of in many fucreate_string_paramsnctions as it was up until now.
// Make sure the calling_function string is properly dealt with at the beginning of this function
// started on Mar 18, 2008
// array is NULL unless the calling function is save_parameters_in_array or fill_data_structures_with_new_parameters_from_array
{
    int index;
    int start;  // used in penalty_by_size
    int i, j, k, l, m, n, o, p;
    index = 0;
    int job;        // used to figure out what kind of work to do; faster to use int
    int ip, jp;
    char type[100];
    int tindex;
    int kalt;   // alternative to k
    char buffer[100];
    double param;
    int line = 0;
    int inc;
    int sim_index = 0;


    // need for save_parameters
    FILE *file;

    if (strcmp (calling_function, "create_string_params") == 0)
        job = 0;
    else if (strcmp (calling_function, "fill_similarity_rules") == 0)
        job = 1;
    else if (strcmp (calling_function, "save_parameters_in_array") == 0)
        job = 2;
    else if (strcmp (calling_function, "save_parameters") == 0)
        job = 3;
    else if (strcmp (calling_function, "fill_data_structures_with_new_parameters") == 0)
        job = 4;
    else if (strcmp (calling_function, "fill_data_structures_with_new_parameters_from_array") == 0)
        job = 5;    // array contains the array with the parameters


    if (job == 3)
    {
        if ((file = fopen (filename, "w")) == NULL)
        {
            giveup ("Cannot open file", filename);
        }
    }
    else if (job == 4)
    {
        if ((file = fopen (filename, "r")) == NULL)
        {
            giveup ("Cannot open file", filename);
        }
    }

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (stack[i][j][k][l] < INF)
                    {
                        // exclude duplicates
                        // stack[i][j][k][l] is the same as stack[l][k][j][i]
                        if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                        {
                            switch (job)
                            {
                                case 0:
                                    sprintf (string_params[index], "stack[%d][%d][%d][%d]", i, j, k, l);
                                    sprintf (string_params_human_readable[index], "stack[5'-%c%c/%c%c-3']",
                                        int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j));
                                    break;
                                case 2:
                                    array[index] = stack[i][j][k][l];
                                    break;
                                case 3:
                                    fprintf (file, "%.2lf\n", (double)stack[i][j][k][l]/100.0);
                                    break;
                                case 4:
                                    fgets (buffer, sizeof(buffer), file);
                                    sscanf (buffer, "%lf\n", &param);
                                    param *= 100;
                                    stack[i][j][k][l] = (PARAMTYPE) param;
                                    // add the duplicate too
                                    stack[l][k][j][i] = (PARAMTYPE) param;
                                    break;
                                case 5:
                                    stack[i][j][k][l] = array[index];
                                    stack[l][k][j][i] = array[index];
                                    break;
                            }
                            index++;    sim_index++;
                        }
                    }
                }

    if (parsi_tstackh == PARSI)     // 1: parsimonious model. Use the same idea as for tstacki
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_AU_closure");
                sprintf (string_params_human_readable[index], "hairpin_AU_GU_closure_penalty");
                break;
            case 2:
                array[index] = misc.hairpin_AU_closure;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_AU_closure/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_AU_closure = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_AU_closure = array[index];
                break;
        }
        index++;    sim_index++;

        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_AG_mismatch");
                sprintf (string_params_human_readable[index], "hairpin_AG_mismatch");
                break;
            case 2:
                array[index] = misc.hairpin_AG_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_AG_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_AG_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_AG_mismatch = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_GA_mismatch");
                sprintf (string_params_human_readable[index], "hairpin_GA_mismatch");
                break;
            case 2:
                array[index] = misc.hairpin_GA_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_GA_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_GA_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_GA_mismatch = array[index];
                break;
        }
        index++;    sim_index++;

        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_UU_mismatch");
                sprintf (string_params_human_readable[index], "hairpin_UU_mismatch");
                break;
            case 2:
                array[index] = misc.hairpin_UU_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_UU_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_UU_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_UU_mismatch = array[index];
                break;
        }
        index++;    sim_index++;
    }
    else if (parsi_tstackh == LAVISH || parsi_tstackh == T99)   // lavish and turner99 model for tstackh
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        if (tstackh[i][j][k][l] < INF)
                        {
                            // no duplicates here
                            switch (job)
                            {
                                case 0:
                                    sprintf (string_params[index], "tstackh[%d][%d][%d][%d]", i, j, k, l);
                                    sprintf (string_params_human_readable[index], "tstackh[5'-%c%c/%c%c-3']",
                                            int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j));
                                    break;
                                case 1:
                                    if (parsi_tstackh == LAVISH)
                                    {
                                        if (similarity_rule[sim_index][0] == '\0')
                                        {
                                            int k_rule1, l_rule1;
                                            if (apply_rule_1 (k, l, k_rule1, l_rule1))
                                                sprintf (similarity_rule[sim_index], "1 * tstackh[5'-%c%c/%c%c-3']",
                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                        }
                                    }
                                    break;
                                case 2:
                                    array[index] = tstackh[i][j][k][l];
                                    break;
                                case 3:
                                    fprintf (file, "%.2lf\n", (double)tstackh[i][j][k][l]/100.0);
                                    break;
                                case 4:
                                    fgets (buffer, sizeof(buffer), file);
                                    sscanf (buffer, "%lf\n", &param);
                                    param *= 100;
                                    tstackh[i][j][k][l] = (PARAMTYPE) param;
                                    break;
                                case 5:
                                    tstackh[i][j][k][l] = array[index];
                                    break;
                            }
                            index++;    sim_index++;
                        }
                    }

    }
    // the new bulge1 parameters

    // for turner99, there's no bulge1
    if (parsi_bulge1 == PARSI || parsi_bulge1 == LAVISH)     // add 4 features or more
    {
        // first take care of bulgeA, bulgeB, bulgeC, bulgeU
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "bulgeA");
                sprintf (string_params_human_readable[index], "bulgeA");
                break;
            case 1:
                if (parsi_bulge1 == LAVISH)
                    similarity_bulge_type (A, similarity_rule[sim_index]);
                break;
            case 2:
                array[index] = bulgeA;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)bulgeA/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                bulgeA = (PARAMTYPE) param;
                break;
            case 5:
                bulgeA = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "bulgeC");
                sprintf (string_params_human_readable[index], "bulgeC");
                break;
            case 1:
                if (parsi_bulge1 == LAVISH)
                    similarity_bulge_type (C, similarity_rule[sim_index]);
                break;
            case 2:
                array[index] = bulgeC;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)bulgeC/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                bulgeC = (PARAMTYPE) param;
                break;
            case 5:
                bulgeC = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "bulgeG");
                sprintf (string_params_human_readable[index], "bulgeG");
                break;
            case 1:
                if (parsi_bulge1 == LAVISH)
                    similarity_bulge_type (G, similarity_rule[sim_index]);
                break;
            case 2:
                array[index] = bulgeG;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)bulgeG/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                bulgeG = (PARAMTYPE) param;
                break;
            case 5:
                bulgeG = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "bulgeU");
                sprintf (string_params_human_readable[index], "bulgeU");
                break;
            case 1:
                if (parsi_bulge1 == LAVISH)
                    similarity_bulge_type (U, similarity_rule[sim_index]);
                break;
            case 2:
                array[index] = bulgeU;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)bulgeU/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                bulgeU = (PARAMTYPE) param;
                break;
            case 5:
                bulgeU = array[index];
                break;
        }
        index++;    sim_index++;
    }

    if (parsi_bulge1 == LAVISH)      // lavish for bulge1
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (ip=0; ip < NUCL; ip++)
                        for (jp=0; jp < NUCL; jp++)
                        {
                            if (bulge1[i][j][k][ip][jp] < INF)
                            {
                                // no duplicates
                                switch (job)
                                {
                                    case 0:
                                        sprintf (string_params[index], "bulge1[%d][%d][%d][%d][%d]", i, j, k, ip, jp);
                                        sprintf (string_params_human_readable[index], "bulge1[5'-%c%c%c/%c%c-3']",
                                            int_to_nuc(i), int_to_nuc(k), int_to_nuc(ip), int_to_nuc(jp), int_to_nuc(j));
                                        break;
                                    case 1:
                                        if (similarity_rule[sim_index][0] == '\0')
                                        {
                                            // make sure I get the right mirror for stack
                                            if (i*1000 + j*100 + ip*10 + jp <= jp*1000 + ip*100 + j*10 + i)
                                                sprintf (similarity_rule[sim_index], "1 * stack[5'-%c%c/%c%c-3'] + 1 * bulge%c",
                                                    int_to_nuc(i), int_to_nuc(ip), int_to_nuc(jp), int_to_nuc(j), int_to_nuc(k));
                                            else
                                                sprintf (similarity_rule[sim_index], "1 * stack[5'-%c%c/%c%c-3'] + 1 * bulge%c",
                                                    int_to_nuc(jp), int_to_nuc(j), int_to_nuc(i), int_to_nuc(ip), int_to_nuc(k));
                                        }
                                        break;
                                    case 2:
                                        array[index] = bulge1[i][j][k][ip][jp];
                                        break;
                                    case 3:
                                        fprintf (file, "%.2lf\n", (double)bulge1[i][j][k][ip][jp]/100.0);
                                        break;
                                    case 4:
                                        fgets (buffer, sizeof(buffer), file);
                                        sscanf (buffer, "%lf\n", &param);
                                        param *= 100;
                                        bulge1[i][j][k][ip][jp] = (PARAMTYPE) param;
                                        break;
                                    case 5:
                                        bulge1[i][j][k][ip][jp] = array[index];
                                        break;
                                }
                                index++;    sim_index++;

                            }
                        }
    }

    // this is included in all cases, at least for internal loops 1xn
    switch (job)
    {
        case 0:
            sprintf (string_params[index], "misc.internal_AU_closure");
            sprintf (string_params_human_readable[index], "internal_AU_GU_closure_penalty");
            break;
        case 2:
            array[index] = misc.internal_AU_closure;
            break;
        case 3:
            fprintf (file, "%.2lf\n", (double)misc.internal_AU_closure/100.0);
            break;
        case 4:
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf\n", &param);
            param *= 100;
            misc.internal_AU_closure = (PARAMTYPE) param;
            break;
        case 5:
            misc.internal_AU_closure = array[index];
            break;
    }
    index++;    sim_index++;

    if (parsi_tstacki == T99)
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_GA_AG_mismatch");
                sprintf (string_params_human_readable[index], "internal_GA_AG_mismatch");
                break;
            case 2:
                array[index] = misc.internal_GA_AG_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_GA_AG_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_GA_AG_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_GA_AG_mismatch = array[index];
                break;
        }
        index++;    sim_index++;
    }
    else if (parsi_tstacki == PARSI)
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_AG_mismatch");
                sprintf (string_params_human_readable[index], "internal_AG_mismatch");
                break;
            case 2:
                array[index] = misc.internal_AG_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_AG_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_AG_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_AG_mismatch = array[index];
                break;
        }
        index++;    sim_index++;

        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_GA_mismatch");
                sprintf (string_params_human_readable[index], "internal_GA_mismatch");
                break;
            case 2:
                array[index] = misc.internal_GA_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_GA_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_GA_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_GA_mismatch = array[index];
                break;
        }
        index++;     sim_index++;

        switch (job)         // see Schroeder_Turner_2000
        {
            case 0:
                sprintf (string_params[index], "misc.internal_GG_mismatch");
                sprintf (string_params_human_readable[index], "internal_GG_mismatch");
                break;
            case 2:
                array[index] = misc.internal_GG_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_GG_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_GG_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_GG_mismatch = array[index];
                break;
        }
        index++;     sim_index++;
    }
    if (parsi_tstacki == PARSI || parsi_tstacki == T99)
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_UU_mismatch");
                sprintf (string_params_human_readable[index], "internal_UU_mismatch");
                break;
            case 2:
                array[index] = misc.internal_UU_mismatch;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_UU_mismatch/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_UU_mismatch = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_UU_mismatch = array[index];
                break;
        }
        index++;    sim_index++;
    }

    if (parsi_tstacki == LAVISH)
    {
        for (i=0; i < NUCL; i++)
        {
            // these depend on tstackh, so parsi_tstacki and parsi_tstackh should be the same
            // took out the following for now
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        if (tstacki[i][j][k][l] < INF)
                        {
                            // no duplicates here
                            switch (job)
                            {
                                case 0:
                                    sprintf (string_params[index], "tstacki[%d][%d][%d][%d]", i, j, k, l);
                                    sprintf (string_params_human_readable[index], "tstacki[5'-%c%c/%c%c-3']",
                                            int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j));
                                    break;
                                case 1:
                                    if (similarity_rule[sim_index][0] == '\0')
                                    {
                                        // follow a dependency on tstacki, not on tstackh!
                                        // follow RULE 1
                                        int k_rule1, l_rule1;
                                        if (apply_rule_1 (k, l, k_rule1, l_rule1))
                                            sprintf (similarity_rule[sim_index], "1 * tstacki[5'-%c%c/%c%c-3']",
                                                int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                        // make them the same as tstackh
                                        else
                                        {
                                            // I used to use tstackh but I changed my mind
                                            //sprintf (similarity_rule[sim_index], "1 * tstackh[5'-%c%c/%c%c-3']",
                                            //    int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j));

                                            // try to make the same as tstacki[j][i][k][l]
                                            char type[100];
                                            int tindex;
                                            sprintf (type, "tstacki[%d][%d][%d][%d]", j, i, k, l);
                                            tindex = structure_type_index(type);
                                            if (similarity_rule[tindex][0] != '\0')
                                                sprintf (similarity_rule[sim_index], "1 * tstacki[5'-%c%c/%c%c-3']",
                                                    int_to_nuc(j), int_to_nuc(k), int_to_nuc(l), int_to_nuc(i));

                                            if (i==G && j==U)       // make it AU.
                                                // It doesn't matter if the one with AU has experimental support
                                            {
                                                sprintf (similarity_rule[sim_index], "1 * tstacki[5'-%c%c/%c%c-3']",
                                                    int_to_nuc(A), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j));
                                            }
                                            else if (i==U && j==G)       // make it UA
                                            {
                                                sprintf (similarity_rule[sim_index], "1 * tstacki[5'-%c%c/%c%c-3']",
                                                    int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(A));
                                            }
                                            if (i==A && j==U)   // replace with GC
                                            {
                                                sprintf (type, "tstacki[%d][%d][%d][%d]", G, C, k, l);
                                                tindex = structure_type_index(type);
                                                if (similarity_rule[tindex][0] != '\0')
                                                    sprintf (similarity_rule[sim_index], "1 * tstacki[5'-%c%c/%c%c-3'] + 1 * internal_AU_GU_closure_penalty",
                                                        int_to_nuc(G), int_to_nuc(k), int_to_nuc(l), int_to_nuc(C));
                                            }
                                            else if (i==U && j==A)   // replace with CG
                                            {
                                                sprintf (type, "tstacki[%d][%d][%d][%d]", C, G, k, l);
                                                tindex = structure_type_index(type);
                                                if (similarity_rule[tindex][0] != '\0')
                                                    sprintf (similarity_rule[sim_index], "1 * tstacki[5'-%c%c/%c%c-3'] + 1 * internal_AU_GU_closure_penalty",
                                                        int_to_nuc(C), int_to_nuc(k), int_to_nuc(l), int_to_nuc(G));
                                            }
                                        }
                                    }
                                    break;
                                case 2:
                                    array[index] = tstacki[i][j][k][l];
                                    break;
                                case 3:
                                    fprintf (file, "%.2lf\n", (double)tstacki[i][j][k][l]/100.0);
                                    break;
                                case 4:
                                    fgets (buffer, sizeof(buffer), file);
                                    sscanf (buffer, "%lf\n", &param);
                                    param *= 100;
                                    tstacki[i][j][k][l] = (PARAMTYPE) param;
                                    break;
                                case 5:
                                    tstacki[i][j][k][l] = array[index];
                                    break;
                            }
                            index++;    sim_index++;
                        }
                    }
        }
    }

    if (!simple_internal_energy)
    {

        if (parsi_int11 == T99)
        {
            // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                {
                                    if ( (((i==C && j==G) || (i==G && j==C)) &&
                                        ((m==C && n==G) || (m==G && n==C)) &&
                                        !can_pair(k,l)) ||
                                        (watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                    {
                                        if (int11[i][j][k][l][m][n] < INF)
                                        {
                                            // exclude duplicates
                                            // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                            if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                            {
                                                switch (job)
                                                {
                                                    case 0:
                                                        sprintf (string_params[index], "int11[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
                                                        sprintf (string_params_human_readable[index], "int11[5'-%c%c%c/%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(l), int_to_nuc(j));
                                                        break;
                                                    case 2:
                                                        array[index] = int11[i][j][k][l][m][n];
                                                        break;
                                                    case 3:
                                                        fprintf (file, "%.2lf\n", (double)int11[i][j][k][l][m][n]/100.0);
                                                        break;
                                                    case 4:
                                                        fgets (buffer, sizeof(buffer), file);
                                                        sscanf (buffer, "%lf\n", &param);
                                                        param *= 100;
                                                        int11[i][j][k][l][m][n] = (PARAMTYPE) param;
                                                        // add the duplicate too
                                                        int11[n][m][l][k][j][i] = (PARAMTYPE) param;
                                                        break;
                                                    case 5:
                                                        int11[i][j][k][l][m][n] = array[index];
                                                        int11[n][m][l][k][j][i] = array[index];
                                                        break;
                                                }
                                                index++;      sim_index++;
                                            }
                                        }
                                    }
                                }
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_basic_mismatch");
                    sprintf (string_params_human_readable[index], "int11_basic_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal11_basic_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_basic_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_basic_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_basic_mismatch = array[index];
                    break;
            }
            index++;     sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_GG_mismatch");
                    sprintf (string_params_human_readable[index], "int11_GG_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal11_GG_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_GG_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_GG_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_GG_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;
        }
        else if (parsi_int11 == PARSI || parsi_int11 == LAVISH || parsi_int11 == HLI)
        {
            // first let's work with the 10 extra-parameters
            // Also, these are the only ones I use if parsi_int11
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_AU_closure");
                    sprintf (string_params_human_readable[index], "internal11_AU_closure");
                    break;
                case 2:
                    array[index] = misc.internal11_AU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_AU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_AU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_AU_closure = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_GU_closure");
                    sprintf (string_params_human_readable[index], "internal11_GU_closure");
                    break;
                case 2:
                    array[index] = misc.internal11_GU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_GU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_GU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_GU_closure = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_AG_mismatch");
                    sprintf (string_params_human_readable[index], "internal11_AG_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal11_AG_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_AG_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_AG_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_AG_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_GG_mismatch");
                    sprintf (string_params_human_readable[index], "internal11_GG_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal11_GG_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_GG_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_GG_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_GG_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_UU_mismatch");
                    sprintf (string_params_human_readable[index], "internal11_UU_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal11_UU_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_UU_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_UU_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_UU_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_5YRR_5YRR");
                    sprintf (string_params_human_readable[index], "internal11_5YRR_5YRR");
                    break;
                case 2:
                    array[index] = misc.internal11_5YRR_5YRR;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_5YRR_5YRR/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_5YRR_5YRR = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_5YRR_5YRR = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_5RYY_5RYY");
                    sprintf (string_params_human_readable[index], "internal11_5RYY_5RYY");
                    break;
                case 2:
                    array[index] = misc.internal11_5RYY_5RYY;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_5RYY_5RYY/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_5RYY_5RYY = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_5RYY_5RYY = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_5YYR_5YYR");
                    sprintf (string_params_human_readable[index], "internal11_5YYR_5YYR");
                    break;
                case 2:
                    array[index] = misc.internal11_5YYR_5YYR;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_5YYR_5YYR/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_5YYR_5YYR = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_5YYR_5YYR = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_5YRY_5RYR");
                    sprintf (string_params_human_readable[index], "internal11_5YRY_5RYR");
                    break;
                case 2:
                    array[index] = misc.internal11_5YRY_5RYR;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_5YRY_5RYR/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_5YRY_5RYR = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_5YRY_5RYR = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal11_5RRY_5RYY");
                    sprintf (string_params_human_readable[index], "internal11_5RRY_5RYY");
                    break;
                case 2:
                    array[index] = misc.internal11_5RRY_5RYY;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal11_5RRY_5RYY/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal11_5RRY_5RYY = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal11_5RRY_5RYY = array[index];
                    break;
            }
            index++;    sim_index++;

            if (parsi_int11 == LAVISH || parsi_int11 == HLI)
            {
                for (i=0; i < NUCL; i++)
                    for (j=0; j < NUCL; j++)
                        for (k=0; k < NUCL; k++)
                            for (l=0; l < NUCL; l++)
                            {
                                //if (watson_crick (k,l)) continue;     // we want all these
                                for (m=0; m < NUCL; m++)
                                    for (n=0; n < NUCL; n++)
                                    {
                                        if (int11[i][j][k][l][m][n] < INF)
                                        {
                                            // exclude duplicates
                                            // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                            if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                            {
                                                inc = 1;
                                                switch (job)
                                                {
                                                    case 0:
                                                        //if (int11_experimental_addition[i][j][k][l][m][n] < INF)
                                                        // it might not get in here, but then it's ok, because we replace the right values in case 2
                                                        if (similarity_rule[sim_index][0] != '\0')
                                                        // this assumes I called the function fill_similarity_rule_with_optical_melting_reference
                                                        {
                                                            sprintf (string_params[index], "int11_experimental_addition[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
                                                            sprintf (string_params_human_readable[index], "int11_expadd[5'-%c%c%c/%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(l), int_to_nuc(j));
                                                            //printf ("IN create_string_params: %s; similarity_rule is %s\n", string_params_human_readable[index], similarity_rule[sim_index]);
                                                            initialize_correct_int11_expadd (i, j, k, l, m, n);
                                                        }
                                                        else if (parsi_int11 == LAVISH)
                                                        {
                                                            sprintf (string_params[index], "int11[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
                                                            sprintf (string_params_human_readable[index], "int11[5'-%c%c%c/%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(l), int_to_nuc(j));
                                                        }
                                                        else if (parsi_int11 == HLI)  inc = 0;
                                                        break;
                                                    case 1:
                                                        if (similarity_rule[sim_index][0] == '\0')  // this should only happen if LAVISH
                                                        {
                                                            int k_rule1, l_rule1;
                                                            int done = 0;
                                                            // we might have the experimental addition instead of just int11
                                                            if (apply_rule_1 (k, l, k_rule1, l_rule1))
                                                            {
                                                                // after the conversion, it might be that we have the wrong "mirror"
                                                                if (i*100000 + j*10000 + k_rule1*1000 + l_rule1*100 + m*10 + n
                                                                    <= n*100000 + m*10000 + l_rule1*1000+ k_rule1*100 + j*10 + i)
                                                                {
                                                                    if (int11_experimental_addition[i][j][k_rule1][l_rule1][m][n] < INF)
                                                                    {
                                                                        sprintf (similarity_rule[sim_index],
                                                                            "1 * int11_expadd[5'-%c%c%c/%c%c%c-3']",
                                                                            int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(m),
                                                                            int_to_nuc(n), int_to_nuc(l_rule1), int_to_nuc(j));
                                                                    }
                                                                    else
                                                                    {
                                                                        sprintf (similarity_rule[sim_index], "1 * int11[5'-%c%c%c/%c%c%c-3']",
                                                                            int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(m),
                                                                            int_to_nuc(n), int_to_nuc(l_rule1), int_to_nuc(j));
                                                                        done = 1;
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                    if (int11_experimental_addition[n][m][l_rule1][k_rule1][j][i] < INF)
                                                                    {
                                                                        sprintf (similarity_rule[sim_index],
                                                                                "1 * int11_expadd[5'-%c%c%c/%c%c%c-3']",
                                                                                    int_to_nuc(n), int_to_nuc(l_rule1), int_to_nuc(j),
                                                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(m));
                                                                    }
                                                                    else
                                                                    {
                                                                        sprintf (similarity_rule[sim_index], "1 * int11[5'-%c%c%c/%c%c%c-3']",
                                                                                    int_to_nuc(n), int_to_nuc(l_rule1), int_to_nuc(j),
                                                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(m));
                                                                        done = 1;
                                                                    }
                                                                }
                                                            }
                                                            if (done) break;
                                                            // look for AU closure
                                                            char plus[5] = "";
                                                            if (similarity_rule[sim_index][0] != '\0')  strcpy (plus, " + ");
                                                            if (((i==A && j==U) || (i==U && j==A)) && ((m==A && n==U) || (m==U && n==A)))
                                                                sprintf (similarity_rule[sim_index], "%s%s2 * internal11_AU_closure", similarity_rule[sim_index], plus);
                                                            else if ((i==A && j==U) || (i==U && j==A))
                                                                sprintf (similarity_rule[sim_index], "%s%s1 * internal11_AU_closure", similarity_rule[sim_index], plus);
                                                            else if ((m==A && n==U) || (m==U && n==A))
                                                                sprintf (similarity_rule[sim_index], "%s%s1 * internal11_AU_closure", similarity_rule[sim_index], plus);
                                                            // look for GU closure
                                                            if (similarity_rule[sim_index][0] != '\0')  strcpy (plus, " + ");
                                                            else    strcpy (plus, "");
                                                            if (((i==G && j==U) || (i==U && j==G)) && ((m==G && n==U) || (m==U && n==G)))
                                                                sprintf (similarity_rule[sim_index], "%s%s2 * internal11_GU_closure", similarity_rule[sim_index], plus);
                                                            else if ((i==G && j==U) || (i==U && j==G))
                                                                sprintf (similarity_rule[sim_index], "%s%s1 * internal11_GU_closure", similarity_rule[sim_index], plus);
                                                            else if ((m==G && n==U) || (m==U && n==G))
                                                                sprintf (similarity_rule[sim_index], "%s%s1 * internal11_GU_closure", similarity_rule[sim_index], plus);
                                                            // look for AG mismatch
                                                            if ((k==A && l==G) || (k==G && l==A))
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_AG_mismatch", similarity_rule[sim_index]);
                                                            }
                                                            // look for GG mismatch
                                                            if (k==G && l==G)
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_GG_mismatch", similarity_rule[sim_index]);
                                                            }
                                                            // look for UU mismatch
                                                            if (k==U && l==U)
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_UU_mismatch", similarity_rule[sim_index]);
                                                            }
                                                            // check the nearest neighbours
                                                            if (isY(i) && isR(j) && isR(k) && isR(l) && isR(m) && isY(n))
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_5YRR_5YRR", similarity_rule[sim_index]);
                                                            }
                                                            if ( isR(i) && isY(j) && isY(k) && isY(l) && isY(m) && isR(n) )
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_5RYY_5RYY", similarity_rule[sim_index]);
                                                            }
                                                            if ( isY(i) && isR(j) && isY(k) && isY(l) && isR(m) && isY(n) )
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_5YYR_5YYR", similarity_rule[sim_index]);
                                                            }
                                                            if ( (isY(i) && isR(j) && isR(k) && isY(l) && isY(m) && isR(n)) ||
                                                                (isR(i) && isY(j) && isY(k) && isR(l) && isR(m) && isY(n)) )
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_5YRY_5RYR", similarity_rule[sim_index]);
                                                            }
                                                            if ( (isR(i) && isY(j) && isR(k) && isY(l) && isY(m) && isR(n)) ||
                                                                (isR(i) && isY(j) && isY(k) && isR(l) && isY(m) && isR(n)) )
                                                            {
                                                                if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s1 * internal11_5RRY_5RYY", similarity_rule[sim_index]);
                                                            }
                                                        }
                                                        //else    // appears in optical melting experiments
                                                        //{
                                                            // Probably I don't need this, because it's done in case job is 1.
                                                        //    int11_experimental_addition[i][j][k][l][m][n] = int11[i][j][k][l][m][n];
                                                            // replace the corresponding int21 in string_params
                                                    //     sprintf (string_params[index], "int11_experimental_addition[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
                                                    //     sprintf (string_params_human_readable[index], "int11_expadd[5'-%c%c%c/%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(l), int_to_nuc(j));
                                                        //}
                                                        break;
                                                    case 2:
                                                        if (int11_experimental_addition[i][j][k][l][m][n] < INF)
                                                            array[index] = int11_experimental_addition[i][j][k][l][m][n];
                                                        else if (parsi_int11 == LAVISH)
                                                            array[index] = int11[i][j][k][l][m][n];
                                                        else if (parsi_int11 == HLI)     inc = 0;
                                                        break;
                                                    case 3:
                                                        if (int11_experimental_addition[i][j][k][l][m][n] < INF)
                                                            fprintf (file, "%.2lf\n", (double)int11_experimental_addition[i][j][k][l][m][n]/100.0);
                                                        else  if (parsi_int11 == LAVISH)
                                                            fprintf (file, "%.2lf\n", (double)int11[i][j][k][l][m][n]/100.0);
                                                        else if (parsi_int11 == HLI)     inc = 0;
                                                        break;
                                                    case 4:
                                                        fgets (buffer, sizeof(buffer), file);
                                                        sscanf (buffer, "%lf\n", &param);
                                                        param *= 100;
                                                        if (similarity_rule[sim_index][0] != '\0')
                                                        // this assumes I called the function fill_similarity_rule_with_optical_melting_reference
                                                        {
                                                            int11_experimental_addition[i][j][k][l][m][n] = (PARAMTYPE) param;
                                                            // also the symmetric one
                                                            int11_experimental_addition[n][m][l][k][j][i] = int11_experimental_addition[i][j][k][l][m][n];
                                                            //printf ("IN fill_data_structures_with_new_params: %s has value %lf\n", string_params_human_readable[index], param);
                                                        }
                                                        else if (parsi_int11 == LAVISH)
                                                        {
                                                            int11[i][j][k][l][m][n] = (PARAMTYPE) param;
                                                            int11[n][m][l][k][j][i] = int11[i][j][k][l][m][n];
                                                        }
                                                        else if (parsi_int11 == HLI)     inc = 0;
                                                        break;
                                                    case 5:
                                                        if (similarity_rule[sim_index][0] != '\0')
                                                        // this assumes I called the function fill_similarity_rule_with_optical_melting_reference
                                                        {
                                                            int11_experimental_addition[i][j][k][l][m][n] = array[index];
                                                            // also the symmetric one
                                                            int11_experimental_addition[n][m][l][k][j][i] = int11_experimental_addition[i][j][k][l][m][n];
                                                            //printf ("IN fill_data_structures_with_new_params: %s has value %lf\n", string_params_human_readable[index], param);
                                                        }
                                                        else if (parsi_int11 == LAVISH)
                                                        {
                                                            int11[i][j][k][l][m][n] = array[index];
                                                            int11[n][m][l][k][j][i] = int11[i][j][k][l][m][n];
                                                        }
                                                        else if (parsi_int11 == HLI)     inc = 0;
                                                        break;
                                                }
                                                if (inc) index++;
                                                //index++;
                                                sim_index++;
                                            }
                                        }
                                    }
                            }
            }
        }

        if (parsi_int21 == T99)
        {
            // NEXT, int21 parameters
            // go with few parameters, as in Mathews et al 1999
            // closed by CG
            i=C; j=G; m=C; n=G;
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (o=0; o < NUCL; o++)
                        if (!can_pair(k,l) && !can_pair(k,o))
                            if (int21[i][j][k][l][m][n][o] < INF)
                            {
                                // no duplicates here
                                switch (job)
                                {
                                    case 0:
                                        sprintf (string_params[index], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                                        sprintf (string_params_human_readable[index], "int21[5'-%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(l), int_to_nuc(j));
                                        break;
                                    case 2:
                                        array[index] = int21[i][j][k][l][m][n][o];
                                        break;
                                    case 3:
                                        fprintf (file, "%.2lf\n", (double)int21[i][j][k][l][m][n][o]/100.0);
                                        break;
                                    case 4:
                                        fgets (buffer, sizeof(buffer), file);
                                        sscanf (buffer, "%lf\n", &param);
                                        param *= 100;
                                        int21[i][j][k][l][m][n][o] = (PARAMTYPE) param;
                                        break;
                                    case 5:
                                        int21[i][j][k][l][m][n][o] = array[index];
                                        break;
                                }
                                index++;    sim_index++;
                            }
            // closed by GC
            i=G; j=C; m=G; n=C;
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (o=0; o < NUCL; o++)
                        if (!can_pair(k,l) && !can_pair(k,o))
                            if (int21[i][j][k][l][m][n][o] < INF)
                            {
                                // no duplicates here
                                switch (job)
                                {
                                    case 0:
                                        sprintf (string_params[index], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                                        sprintf (string_params_human_readable[index], "int21[5'-%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(l), int_to_nuc(j));
                                        break;
                                    case 2:
                                        array[index] = int21[i][j][k][l][m][n][o];
                                        break;
                                    case 3:
                                        fprintf (file, "%.2lf\n", (double)int21[i][j][k][l][m][n][o]/100.0);
                                        break;
                                    case 4:
                                        fgets (buffer, sizeof(buffer), file);
                                        sscanf (buffer, "%lf\n", &param);
                                        param *= 100;
                                        int21[i][j][k][l][m][n][o] = (PARAMTYPE) param;
                                        break;
                                    case 5:
                                        int21[i][j][k][l][m][n][o] = array[index];
                                        break;
                                }
                                index++;    sim_index++;
                            }
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_match");
                    sprintf (string_params_human_readable[index], "int21_match");
                    break;
                case 2:
                    array[index] = misc.internal21_match;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_match/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_match = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_match = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_AU_closure");
                    sprintf (string_params_human_readable[index], "int21_AU_closure");
                    break;
                case 2:
                    array[index] = misc.internal21_AU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_AU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_AU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_AU_closure = array[index];
                    break;
            }
            index++;    sim_index++;
        }
        else if (parsi_int21 == PARSI || parsi_int21 == LAVISH)
        {
            // first let's work with the 6 extra-parameters
            // Also, these are the only ones I need if parsi_int21
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_initiation");
                    sprintf (string_params_human_readable[index], "internal21_initiation");
                    break;
                case 2:
                    array[index] = misc.internal21_initiation;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_initiation/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_initiation = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_initiation = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_AU_closure");
                    sprintf (string_params_human_readable[index], "internal21_AU_closure");
                    break;
                case 2:
                    array[index] = misc.internal21_AU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_AU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_AU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_AU_closure = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_GU_closure");
                    sprintf (string_params_human_readable[index], "internal21_GU_closure");
                    break;
                case 2:
                    array[index] = misc.internal21_GU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_GU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_GU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_GU_closure = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_AG_mismatch");
                    sprintf (string_params_human_readable[index], "internal21_AG_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal21_AG_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_AG_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_AG_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_AG_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_GG_mismatch");
                    sprintf (string_params_human_readable[index], "internal21_GG_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal21_GG_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_GG_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_GG_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_GG_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal21_UU_mismatch");
                    sprintf (string_params_human_readable[index], "internal21_UU_mismatch");
                    break;
                case 2:
                    array[index] = misc.internal21_UU_mismatch;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal21_UU_mismatch/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal21_UU_mismatch = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal21_UU_mismatch = array[index];
                    break;
            }
            index++;    sim_index++;

            if (parsi_int21 == LAVISH)
            {
                for (i=0; i < NUCL; i++)
                    for (j=0; j < NUCL; j++)
                        for (k=0; k < NUCL; k++)
                            for (l=0; l < NUCL; l++)
                                for (m=0; m < NUCL; m++)
                                    for (n=0; n < NUCL; n++)
                                        for(o=0; o < NUCL; o++)
                                        {
                                            if (int21[i][j][k][l][m][n][o] < INF)
                                            {
                                                // no duplicates here
                                                switch (job)
                                                {
                                                    case 0:
                                                        //if (int21_experimental_addition[i][j][k][l][m][n][o] < INF)
                                                        // it might not get in here, but then it's ok, because we replace the right values in case 2
                                                        if (similarity_rule[sim_index][0] != '\0')
                                                        // this assumes I called the function fill_similarity_rule_with_optical_melting_reference
                                                        {
                                                            sprintf (string_params[index], "int21_experimental_addition[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                                                            sprintf (string_params_human_readable[index], "int21_expadd[5'-%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(l), int_to_nuc(j));
                                                            initialize_correct_int21_expadd (i, j, k, l, m, n, o);
                                                        }
                                                        else
                                                        {
                                                            sprintf (string_params[index], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                                                            sprintf (string_params_human_readable[index], "int21[5'-%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(l), int_to_nuc(j));
                                                        }
                                                        break;
                                                    case 1:
                                                        if (similarity_rule[sim_index][0] == '\0')
                                                        {
                                                            int k_rule1, l_rule1, o_rule1;
                                                            int done = 0;
                                                            // we might have the experimental addition instead of just int11
                                                            if (apply_rule_1 (k, l, k_rule1, l_rule1) || apply_rule_1 (k, o, k_rule1, o_rule1))
                                                            {
                                                                // no mirror here
                                                                if (int21_experimental_addition[i][j][k_rule1][l_rule1][m][n][o_rule1] < INF)
                                                                {
                                                                    sprintf (similarity_rule[sim_index],
                                                                            "1 * int21_expadd[5'-%c%c%c/%c%c%c%c-3'] + ",
                                                                                int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(m),
                                                                                int_to_nuc(n), int_to_nuc(o_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                                                }
                                                                else
                                                                {
                                                                    sprintf (similarity_rule[sim_index], "1 * int21[5'-%c%c%c/%c%c%c%c-3']",
                                                                            int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(m),
                                                                                int_to_nuc(n), int_to_nuc(o_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                                                    done = 1;
                                                                }
                                                            }
                                                            if (done) break;
                                                            sprintf (similarity_rule[sim_index], "1 * internal21_initiation");
                                                            // look for AU closure
                                                            if (((i==A && j==U) || (i==U && j==A)) && ((m==A && n==U) || (m==U && n==A)))
                                                                sprintf (similarity_rule[sim_index], "%s + 2 * internal21_AU_closure", similarity_rule[sim_index]);
                                                            else if ((i==A && j==U) || (i==U && j==A))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_AU_closure", similarity_rule[sim_index]);
                                                            else if ((m==A && n==U) || (m==U && n==A))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_AU_closure", similarity_rule[sim_index]);
                                                            // look for GU closure
                                                            if (((i==G && j==U) || (i==U && j==G)) && ((m==G && n==U) || (m==U && n==G)))
                                                                sprintf (similarity_rule[sim_index], "%s + 2 * internal21_GU_closure", similarity_rule[sim_index]);
                                                            else if ((i==G && j==U) || (i==U && j==G))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_GU_closure", similarity_rule[sim_index]);
                                                            else if ((m==G && n==U) || (m==U && n==G))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_GU_closure", similarity_rule[sim_index]);
                                                            // look for AG mismatch - but not applied to 5'RA/3'YG loops
                                                            if ((k==A && l==G &&   i!=A && i!=G   &&   j!=U && j!=C) ||
                                                                (k==G && l==A) ||
                                                                (k==G && o==A &&   m!=U && m!=C   &&   n!=A && n!=G) ||
                                                                (k==A && o==G))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_AG_mismatch", similarity_rule[sim_index]);
                                                            // look for GG mismatch
                                                            if (k==G && (l==G || o==G))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_GG_mismatch", similarity_rule[sim_index]);
                                                            // look for UU mismatch
                                                            if (k==U && (l==U || o==U))
                                                                sprintf (similarity_rule[sim_index], "%s + 1 * internal21_UU_mismatch", similarity_rule[sim_index]);
                                                        }
                                                        //else    // appears in optical melting experiments
                                                        //{
                                                        //    int21_experimental_addition[i][j][k][l][m][n][o] = int21[i][j][k][l][m][n][o];
                                                            // replace the corresponding int21 in string_params
                                                        //    sprintf (string_params[index], "int21_experimental_addition[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                                                        //    sprintf (string_params_human_readable[index], "int21_expadd[5'-%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(m), int_to_nuc(n), int_to_nuc(o), int_to_nuc(l), int_to_nuc(j));
                                                        //}
                                                        break;
                                                    case 2:
                                                        if (int21_experimental_addition[i][j][k][l][m][n][o] < INF)
                                                            array[index] = int21_experimental_addition[i][j][k][l][m][n][o];
                                                        else
                                                            array[index] = int21[i][j][k][l][m][n][o];
                                                        break;
                                                    case 3:
                                                        if (int21_experimental_addition[i][j][k][l][m][n][o] < INF)
                                                            fprintf (file, "%.2lf\n", (double)int21_experimental_addition[i][j][k][l][m][n][o]/100.0);
                                                        else
                                                            fprintf (file, "%.2lf\n", (double)int21[i][j][k][l][m][n][o]/100.0);
                                                        break;
                                                    case 4:
                                                        fgets (buffer, sizeof(buffer), file);
                                                        sscanf (buffer, "%lf\n", &param);
                                                        param *= 100;
                                                        //if (int21_experimental_addition[i][j][k][l][m][n][o] < INF)
                                                        if (similarity_rule[sim_index][0] != '\0')                  // there's no symmetri here
                                                            int21_experimental_addition[i][j][k][l][m][n][o] = (PARAMTYPE) param;
                                                        else
                                                            int21[i][j][k][l][m][n][o] = (PARAMTYPE) param;
                                                        break;
                                                    case 5:
                                                        if (similarity_rule[sim_index][0] != '\0')                  // there's no symmetri here
                                                            int21_experimental_addition[i][j][k][l][m][n][o] = array[index];
                                                        else
                                                            int21[i][j][k][l][m][n][o] = array[index];

                                                        break;
                                                }
                                                index++;    sim_index++;
                                            }
                                        }
            }       // end if (!parsi_int21)
        }

        // NEXT, INT2x2
        if (parsi_int22 == T99)
        {
            // 53 params instead of all, as in Mathews et al 1999
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                        {
                            n = i;
                            m = j;
                            p = k;
                            o = l;
                            if (watson_crick(i,j) && !watson_crick(k,l))
                            {
                                // exclude duplicates
                                // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                    n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                {
                                    switch (job)
                                    {
                                        case 0:
                                            sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                            sprintf (string_params_human_readable[index], "int22[5'-%c%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(l), int_to_nuc(j));
                                            break;
                                        case 2:
                                            array[index] = int22[i][j][k][l][m][n][o][p];
                                            break;
                                        case 3:
                                            fprintf (file, "%.2lf\n", (double)int22[i][j][k][l][m][n][o][p]/100.0);
                                            break;
                                        case 4:
                                            fgets (buffer, sizeof(buffer), file);
                                            sscanf (buffer, "%lf\n", &param);
                                            param *= 100;
                                            int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
                                            // now the duplicate
                                            int22[n][m][p][o][j][i][l][k] = (PARAMTYPE) param;
                                            break;
                                        case 5:
                                            int22[i][j][k][l][m][n][o][p] = array[index];
                                            int22[n][m][p][o][j][i][l][k] = array[index];
                                            break;
                                    }
                                    index++;    sim_index++;
                                }
                            }
                        }
            // then add the 4 deltas
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_delta_same_size");
                    sprintf (string_params_human_readable[index], "int22_delta_same_size");
                    break;
                case 2:
                    array[index] = misc.internal22_delta_same_size;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_delta_same_size/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_delta_same_size = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_delta_same_size = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_delta_different_size");
                    sprintf (string_params_human_readable[index], "int22_delta_different_size");
                    break;
                case 2:
                    array[index] = misc.internal22_delta_different_size;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_delta_different_size/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_delta_different_size = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_delta_different_size = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_delta_1stable_1unstable");
                    sprintf (string_params_human_readable[index], "int22_delta_1stable_1unstable");
                    break;
                case 2:
                    array[index] = misc.internal22_delta_1stable_1unstable;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_delta_1stable_1unstable/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_delta_1stable_1unstable = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_delta_1stable_1unstable = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_delta_AC");
                    sprintf (string_params_human_readable[index], "int22_delta_AC");
                    break;
                case 2:
                    array[index] = misc.internal22_delta_AC;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_delta_AC/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_delta_AC = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_delta_AC = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_match");
                    sprintf (string_params_human_readable[index], "int22_match");
                    break;
                case 2:
                    array[index] = misc.internal22_match;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_match/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_match = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_match = array[index];
                    break;
            }
            index++;    sim_index++;
        }
        else if (parsi_int22 == PARSI || parsi_int22 == LAVISH)
        {
            // I follow the model suggested in Christiansen_Znosko_2008
            // the basic 6 parameters
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22mid_group1");
                    sprintf (string_params_human_readable[index], "int22mid_group1");
                    break;
                case 2:
                    array[index] = misc.internal22mid_group1;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22mid_group1/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22mid_group1 = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22mid_group1 = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22mid_group2");
                    sprintf (string_params_human_readable[index], "int22mid_group2");
                    break;
                case 2:
                    array[index] = misc.internal22mid_group2;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22mid_group2/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22mid_group2 = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22mid_group2 = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22mid_group3");
                    sprintf (string_params_human_readable[index], "int22mid_group3");
                    break;
                case 2:
                    array[index] = misc.internal22mid_group3;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22mid_group3/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22mid_group3 = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22mid_group3 = array[index];
                    break;
            }
            index++;    sim_index++;
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22mid_group4");
                    sprintf (string_params_human_readable[index], "int22mid_group4");
                    break;
                case 2:
                    array[index] = misc.internal22mid_group4;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22mid_group4/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22mid_group4 = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22mid_group4 = array[index];
                    break;
            }
            index++;    sim_index++;

            // also add misc.internal22_AU_closure
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_AU_closure");
                    sprintf (string_params_human_readable[index], "internal22_AU_closure");
                    break;
                case 2:
                    array[index] = misc.internal22_AU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_AU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_AU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_AU_closure = array[index];
                    break;
            }
            index++;    sim_index++;
            // also add misc.internal22_GU_closure
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_GU_closure");
                    sprintf (string_params_human_readable[index], "internal22_GU_closure");
                    break;
                case 2:
                    array[index] = misc.internal22_GU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)misc.internal22_GU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_GU_closure = (PARAMTYPE) param;
                    break;
                case 5:
                    misc.internal22_GU_closure = array[index];
                    break;
            }
            index++;    sim_index++;

            if (parsi_int22 == LAVISH)
            {
                // add the asymmetric int22
                for (i=0; i < NUCL; i++)
                    for (j=0; j < NUCL; j++)
                    {
                        if (!can_pair(i,j)) continue;
                        for (k=0; k < NUCL; k++)
                            for (l=0; l < NUCL; l++)
                            {
                                // for now, let's only include ncbp in the internal loop
                                //if (watson_crick(k,l)) continue;    // we need to include all
                                for (m=0; m < NUCL; m++)
                                    for (n=0; n < NUCL; n++)
                                    {
                                        if (!can_pair(m,n)) continue;
                                        for(o=0; o < NUCL; o++)
                                            for (p=0; p < NUCL; p++)
                                            {
                                                //if (watson_crick(o,p)) continue;    // we need to include all
                                                // exclude duplicates
                                                // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                                if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                                    n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                                {
                                                    switch (job)
                                                    {
                                                        case 0:
                                                            if (similarity_rule[sim_index][0] != '\0')
                                                            // this assumes I called the function fill_similarity_rule_with_optical_melting_reference
                                                            {
                                                                sprintf (string_params[index], "int22_experimental_addition[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                                                sprintf (string_params_human_readable[index], "int22_expadd[5'-%c%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(l), int_to_nuc(j));
                                                                initialize_correct_int22_expadd (i, j, k, l, m, n, o, p);
                                                            }
                                                            else
                                                            {
                                                                sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                                                sprintf (string_params_human_readable[index], "int22[5'-%c%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(l), int_to_nuc(j));
                                                            }
                                                            break;
                                                        case 1:
                                                            if (similarity_rule[sim_index][0] == '\0')
                                                            {
                                                                int k_rule1, l_rule1, o_rule1, p_rule1;
                                                                int done = 0;
                                                            // we might have the experimental addition instead of just int11
                                                                if (apply_rule_1 (k, l, k_rule1, l_rule1) || (apply_rule_1 (o, p, o_rule1, p_rule1)))
                                                                {
                                                                // after the conversion, it might be that we have the wrong "mirror"
                                                                    if (i*10000000 + j*1000000 + k_rule1*100000 + l_rule1*10000 + m*1000 + n*100 + o_rule1*10 + p_rule1 <=
                                                                        n*10000000 + m*1000000 + p_rule1*100000 + o_rule1*10000 + j*1000 + i*100 + l_rule1*10 + k_rule1)
                                                                    {
                                                                        if (int22_experimental_addition[i][j][k_rule1][l_rule1][m][n][o_rule1][p_rule1] < INF)
                                                                        {
                                                                            sprintf (similarity_rule[sim_index],
                                                                                    "1 * int22_expadd[5'-%c%c%c%c/%c%c%c%c-3']",
                                                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(o_rule1), int_to_nuc(m),
                                                                                    int_to_nuc(n), int_to_nuc(p_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                                                        }
                                                                        else
                                                                        {
                                                                            sprintf (similarity_rule[sim_index], "1 * int22[5'-%c%c%c%c/%c%c%c%c-3']",
                                                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(o_rule1), int_to_nuc(m),
                                                                                    int_to_nuc(n), int_to_nuc(p_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                                                            done = 1;
                                                                        }
                                                                    }
                                                                    else
                                                                    {
                                                                        if (int22_experimental_addition[n][m][p_rule1][o_rule1][j][i][l_rule1][k_rule1] < INF)
                                                                        {
                                                                            sprintf (similarity_rule[sim_index],
                                                                                    "1 * int22_expadd[5'-%c%c%c%c/%c%c%c%c-3']",
                                                                                    int_to_nuc(n), int_to_nuc(p_rule1), int_to_nuc(l_rule1), int_to_nuc(j),
                                                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(o_rule1), int_to_nuc(m));
                                                                        }
                                                                        else
                                                                        {
                                                                            sprintf (similarity_rule[sim_index], "1 * int22[5'-%c%c%c%c/%c%c%c%c-3']",
                                                                                    int_to_nuc(n), int_to_nuc(p_rule1), int_to_nuc(l_rule1), int_to_nuc(j),
                                                                                    int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(o_rule1), int_to_nuc(m));
                                                                            done = 1;
                                                                        }
                                                                    }
                                                                }
                                                                if (done) break;

                                                                char plus[5] = "";
                                                                if (similarity_rule[sim_index][0] != '\0')  strcpy (plus, " + ");
                                                                if (((i==A && j==U) || (i==U && j==A)) && ((m==A && n==U) || (m==U && n==A)))
                                                                    sprintf (similarity_rule[sim_index], "%s%s2 * internal22_AU_closure", similarity_rule[sim_index], plus);
                                                                else if ((i==A && j==U) || (i==U && j==A))
                                                                    sprintf (similarity_rule[sim_index], "%s%s1 * internal22_AU_closure", similarity_rule[sim_index], plus);
                                                                else if ((m==A && n==U) || (m==U && n==A))
                                                                    sprintf (similarity_rule[sim_index], "%s%s1 * internal22_AU_closure", similarity_rule[sim_index], plus);
                                                                // look for GU closure
                                                                if (((i==G && j==U) || (i==U && j==G)) || ((m==G && n==U) || (m==U && n==G)))
                                                                    if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                if (((i==G && j==U) || (i==U && j==G)) && ((m==G && n==U) || (m==U && n==G)))
                                                                    sprintf (similarity_rule[sim_index], "%s2 * internal22_GU_closure", similarity_rule[sim_index]);
                                                                else if ((i==G && j==U) || (i==U && j==G))
                                                                    sprintf (similarity_rule[sim_index], "%s1 * internal22_GU_closure", similarity_rule[sim_index]);
                                                                else if ((m==G && n==U) || (m==U && n==G))
                                                                    sprintf (similarity_rule[sim_index], "%s1 * internal22_GU_closure", similarity_rule[sim_index]);
                                                                if (is_int22_group_1 (k,l,o,p))
                                                                {
                                                                    if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                    sprintf (similarity_rule[sim_index], "%s1 * int22mid_group1", similarity_rule[sim_index]);
                                                                }
                                                                else if (is_int22_group_2 (k,l,o,p))
                                                                {
                                                                    if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                    sprintf (similarity_rule[sim_index], "%s1 * int22mid_group2", similarity_rule[sim_index]);
                                                                }
                                                                else if (is_int22_group_3 (k,l,o,p))
                                                                {
                                                                    if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                    sprintf (similarity_rule[sim_index], "%s1 * int22mid_group3", similarity_rule[sim_index]);
                                                                }
                                                                else if (is_int22_group_4 (k,l,o,p))
                                                                {
                                                                    if (similarity_rule[sim_index][0] != '\0')  sprintf (similarity_rule[sim_index], "%s + ", similarity_rule[sim_index]);
                                                                    sprintf (similarity_rule[sim_index], "%s1 * int22mid_group4", similarity_rule[sim_index]);
                                                                }
                                                            }
                                                            //else    // appears in optical melting experiments
                                                            //{
                                                            //    int22_experimental_addition[i][j][k][l][m][n][o][p] = int22[i][j][k][l][m][n][o][p];
                                                                // replace the corresponding int22 in string_params
                                                            //    sprintf (string_params[index], "int22_experimental_addition[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                                            //    sprintf (string_params_human_readable[index], "int22_expadd[5'-%c%c%c%c/%c%c%c%c-3']",
                                                            //        int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(l), int_to_nuc(j));
                                                            //}
                                                            break;
                                                        case 2:
                                                            if (int22_experimental_addition[i][j][k][l][m][n][o][p] < INF)
                                                                array[index] = int22_experimental_addition[i][j][k][l][m][n][o][p];
                                                            else
                                                                array[index] = int22[i][j][k][l][m][n][o][p];
                                                            break;
                                                        case 3:
                                                            if (int22_experimental_addition[i][j][k][l][m][n][o][p] < INF)
                                                                fprintf (file, "%.2lf\n", (double)int22_experimental_addition[i][j][k][l][m][n][o][p]/100.0);
                                                            else
                                                                fprintf (file, "%.2lf\n", (double)int22[i][j][k][l][m][n][o][p]/100.0);
                                                            break;
                                                        case 4:
                                                            fgets (buffer, sizeof(buffer), file);
                                                            sscanf (buffer, "%lf\n", &param);
                                                            param *= 100;
                                                            if (similarity_rule[sim_index][0] != '\0')
                                                            {
                                                                int22_experimental_addition[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
                                                                // now the duplicate
                                                                int22_experimental_addition[n][m][p][o][j][i][l][k] = int22_experimental_addition[i][j][k][l][m][n][o][p];
                                                            }
                                                            else
                                                            {
                                                                int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
                                                                // now the duplicate
                                                                int22[n][m][p][o][j][i][l][k] = int22[i][j][k][l][m][n][o][p];
                                                            }
                                                            break;
                                                        case 5:
                                                            if (similarity_rule[sim_index][0] != '\0')
                                                            {
                                                                int22_experimental_addition[i][j][k][l][m][n][o][p] = array[index];
                                                                // now the duplicate
                                                                int22_experimental_addition[n][m][p][o][j][i][l][k] = int22_experimental_addition[i][j][k][l][m][n][o][p];
                                                            }
                                                            else
                                                            {
                                                                int22[i][j][k][l][m][n][o][p] = array[index];
                                                                // now the duplicate
                                                                int22[n][m][p][o][j][i][l][k] = int22[i][j][k][l][m][n][o][p];
                                                            }
                                                            break;
                                                    }
                                                    //printf ("index=%d; int22_expadd[C][G][A][A][C][G][A][C] = %.2Lf\n", index, int22_experimental_addition[C][G][A][A][C][G][A][C]);
                                                    index++;    sim_index++;
                                                }
                                            }
                                    }
                            }
                    }
            }   // end if (!parsi_int22)
        }

    }    // end if (!simple_internal_energy)

    if (parsi_dangles == T99 || parsi_dangles == LAVISH)
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                {
                    if (dangle_top[i][j][k] < INF)
                    {
                        // no duplicates here
                        switch (job)
                        {
                            case 0:
                                sprintf (string_params[index], "dangle_top[%d][%d][%d]", i, j, k);
                                sprintf (string_params_human_readable[index], "dangle3[5'-%c/%c%c-3']", int_to_nuc(j), int_to_nuc(i), int_to_nuc(k));
                                break;
                            case 2:
                                array[index] = dangle_top[i][j][k];
                                break;
                            case 3:
                                fprintf (file, "%.2lf\n", (double)dangle_top[i][j][k]/100.0);
                                break;
                            case 4:
                                fgets (buffer, sizeof(buffer), file);
                                sscanf (buffer, "%lf\n", &param);
                                param *= 100;
                                dangle_top[i][j][k] = (PARAMTYPE) param;
                                break;
                            case 5:
                                dangle_top[i][j][k] = array[index];
                                break;
                        }
                        index++;    sim_index++;
                    }
                }
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                {
                    if (dangle_bot[i][j][k] < INF)
                    {
                        // no duplicates here
                        switch (job)
                        {
                            case 0:
                                sprintf (string_params[index], "dangle_bot[%d][%d][%d]", i, j, k);
                                sprintf (string_params_human_readable[index], "dangle5[5'-%c%c/%c-3']", int_to_nuc(k), int_to_nuc(j), int_to_nuc(i));
                                break;
                            case 2:
                                array[index] = dangle_bot[i][j][k];
                                break;
                            case 3:
                                fprintf (file, "%.2lf\n", (double)dangle_bot[i][j][k]/100.0);
                                break;
                            case 4:
                                fgets (buffer, sizeof(buffer), file);
                                sscanf (buffer, "%lf\n", &param);
                                param *= 100;
                                dangle_bot[i][j][k] = (PARAMTYPE) param;
                                break;
                            case 5:
                                dangle_bot[i][j][k] = array[index];
                                break;
                        }
                        index++;    sim_index++;
                    }
                }
    }

    /*
    #if (MODEL == EXTENDED)
    // use it as a parameter only if !parsi_length, otherwise consider it's fixed
    // ACTUALLY consider it fixed always
    if (!parsi_length)
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.param_greater30");
                sprintf (string_params_human_readable[index], "extrapolation_large_loops");
                break;
            case 2:
                array[index] = misc.param_greater30*100.0;      // multiply by 100 because all the other params are stored as value*100
                break;
            case 3:
                fprintf (file, "%.2lf\n", misc.param_greater30);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                // !! don't multiply by 100 here, but do everywhere else
                misc.param_greater30 = (PARAMTYPE) param;
                break;
        }
        index++;    sim_index++;
    }
    #endif
    */

    int end;
    if (parsi_length == T99)
    {
        if (!simple_internal_energy)    start = 4;
        else                            start = 1;
        end = MAXLOOP_I_T99;
    }
    else if (parsi_length == PARSI || parsi_length == LAVISH || parsi_length == ZL)
    {
        if (parsi_int11 == PARSI)       start = 2;
        else if (parsi_int21 == PARSI)  start = 3;
        else                            start = 4;
        if (parsi_length == PARSI || parsi_length == ZL)      end = MAXLOOP_I_PARSI;
        else                            end = MAXLOOP_I_LAVISH;
    }

    for (i=start; i <= end; i++)
    {
        if (internal_penalty_by_size[i] < INF)
        {
            // no duplicates here
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "internal_penalty_by_size[%d]", i);
                    sprintf (string_params_human_readable[index], "internal_size[%d]", i);
                    break;
                case 1:
                    if (similarity_rule[sim_index][0] == '\0')
                    {
                        sprintf (similarity_rule[sim_index], "1 * internal_size[%d] + %.4lf", i-1, misc.param_greater30*log(1.0*i/(i-1)));
                    }
                    break;
                case 2:
                    array[index] = internal_penalty_by_size[i];
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)internal_penalty_by_size[i]/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    internal_penalty_by_size[i] = (PARAMTYPE) param;
                    break;
                case 5:
                    internal_penalty_by_size[i] = array[index];
                    break;
            }
            index++;    sim_index++;
        }
    }

    // I tried to do internal loop 2D, but I don't think it makes too much sense to do it. See the code at the very end of this file

    // NOW internal asymmetries
    if (parsi_asymmetry == PARSI || parsi_asymmetry == LAVISH)
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "internal_asymmetry_initiation");
                sprintf (string_params_human_readable[index], "internal_asymmetry_initiation");
                break;
            case 2:
                array[index] = internal_asymmetry_initiation;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)internal_asymmetry_initiation/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                internal_asymmetry_initiation = (PARAMTYPE) param;
                break;
            case 5:
                internal_asymmetry_initiation = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "internal_asymmetry_slope");
                sprintf (string_params_human_readable[index], "internal_asymmetry_slope");
                break;
            case 2:
                array[index] = internal_asymmetry_slope;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)internal_asymmetry_slope/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                internal_asymmetry_slope = (PARAMTYPE) param;
                break;
            case 5:
                internal_asymmetry_slope = array[index];
                break;
        }
        index++;    sim_index++;
        if (parsi_asymmetry == LAVISH)
        {
            // have individual values for up to MAXLOOP-1
            for (i=1; i < MAXLOOP_ASYM; i++)
            {
                switch (job)
                {
                    case 0:
                        sprintf (string_params[index], "internal_asymmetry[%d]", i);
                        sprintf (string_params_human_readable[index], "internal_asymmetry[%d]", i);
                        break;
                    case 1:
                        if (similarity_rule[sim_index][0] == '\0')
                        {
                            sprintf (similarity_rule[sim_index], "1 * internal_asymmetry_initiation + %.4lf * internal_asymmetry_slope", log(i*1.0));
                        }
                        break;
                    case 2:
                        array[index] = internal_asymmetry[i];
                        break;
                    case 3:
                        fprintf (file, "%.2lf\n", (double)internal_asymmetry[i]/100.0);
                        break;
                    case 4:
                        fgets (buffer, sizeof(buffer), file);
                        sscanf (buffer, "%lf\n", &param);
                        param *= 100;
                        internal_asymmetry[i] = (PARAMTYPE) param;
                        break;
                    case 5:
                        internal_asymmetry[i] = array[index];
                        break;
                }
                index++;    sim_index++;
            }
        }
    }

    if (parsi_length == T99)
    {
        if (parsi_bulge1 == T99)     start = 1;
        else                         start = 2;
        end = MAXLOOP_B_T99;
    }
    else if (parsi_length == PARSI || parsi_length == LAVISH || parsi_length == ZL)
    {
        if (parsi_bulge1 == T99)     start = 1;
        else                         start = 2;
        if (parsi_length == PARSI || parsi_length == ZL)   end = MAXLOOP_B_PARSI;
        else                         end = MAXLOOP_B_LAVISH;
    }

    for (i=start; i <= end; i++)
    {
        if (bulge_penalty_by_size[i] < INF)
        {
            // no duplicates here
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "bulge_penalty_by_size[%d]", i);
                    sprintf (string_params_human_readable[index], "bulge_size[%d]", i);
                    break;
                case 1:
                    if (similarity_rule[sim_index][0] == '\0')
                    {
                        sprintf (similarity_rule[sim_index], "1 * bulge_size[%d] + %.4lf", i-1, misc.param_greater30*log(1.0*i/(i-1)));
                    }
                    break;
                case 2:
                    array[index] = bulge_penalty_by_size[i];
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)bulge_penalty_by_size[i]/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    bulge_penalty_by_size[i] = (PARAMTYPE) param;
                    break;
                case 5:
                    bulge_penalty_by_size[i] = array[index];
                    break;
            }
            index++;    sim_index++;
        }
    }

    // HAIRPIN LOOP LENGTH
    if (parsi_length == T99)            end = MAXLOOP_H_T99;
    else if (parsi_length == PARSI || parsi_length == ZL)     end = MAXLOOP_H_PARSI;
    else if (parsi_length == LAVISH)    end = MAXLOOP_H_LAVISH;

    for (i=1; i <= end; i++)
    {
        if (hairpin_penalty_by_size[i] < INF)
        {
            // no duplicates here
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "hairpin_penalty_by_size[%d]", i);
                    sprintf (string_params_human_readable[index], "hairpin_size[%d]", i);
                    break;
                case 1:
                    if (similarity_rule[sim_index][0] == '\0')
                    {
                        sprintf (similarity_rule[sim_index], "1 * hairpin_size[%d] + %.4lf", i-1, misc.param_greater30*log(1.0*i/(i-1)));
                    }
                    break;
                case 2:
                    array[index] = hairpin_penalty_by_size[i];
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)hairpin_penalty_by_size[i]/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    hairpin_penalty_by_size[i] = (PARAMTYPE) param;
                    break;
                case 5:
                    hairpin_penalty_by_size[i] = array[index];
                    break;
            }
            index++;    sim_index++;
        }
    }

    switch (job)
    {
        case 0:
            sprintf (string_params[index], "misc.terminal_AU_penalty");
            sprintf (string_params_human_readable[index], "terminal_AU_penalty");
            break;
        case 2:
            array[index] = misc.terminal_AU_penalty;
            break;
        case 3:
            fprintf (file, "%.2lf\n", (double)misc.terminal_AU_penalty/100.0);
            break;
        case 4:
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf\n", &param);
            param *= 100;
            misc.terminal_AU_penalty = (PARAMTYPE) param;
            break;
        case 5:
            misc.terminal_AU_penalty = array[index];
            break;
    }
    index++;    sim_index++;

    // SPECIAL HAIRPIN LOOPS
    if (parsi_special == T99 || parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_GGG");
                sprintf (string_params_human_readable[index], "hairpin_GGG");
                break;
            case 2:
                array[index] = misc.hairpin_GGG;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_GGG/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_GGG = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_GGG = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_c1");
                sprintf (string_params_human_readable[index], "hairpin_c1");
                break;
            case 2:
                array[index] = misc.hairpin_c1;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_c1/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_c1 = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_c1 = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_c2");
                sprintf (string_params_human_readable[index], "hairpin_c2");
                break;
            case 2:
                array[index] = misc.hairpin_c2;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_c2/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_c2 = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_c2 = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.hairpin_c3");
                sprintf (string_params_human_readable[index], "hairpin_c3");
                break;
            case 2:
                array[index] = misc.hairpin_c3;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.hairpin_c3/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.hairpin_c3 = (PARAMTYPE) param;
                break;
            case 5:
                misc.hairpin_c3 = array[index];
                break;
        }
        index++;    sim_index++;
    }

    // TODO
    //sprintf (string_params[index], "misc.asymmetry_penalty_max_correction");
    //sprintf (string_params[index], "misc.asymmetry_penalty_array[0]");
    //sprintf (string_params[index], "misc.asymmetry_penalty_array[1]");

    // sprintf (string_params[index], "misc.asymmetry_penalty_array[2]");
    // sprintf (string_params[index], "misc.asymmetry_penalty_array[3]");
    // Instead of these, I will just store the asymmetry for 0.5, 1, 1.5, 2, 2.5 and 3.
    //sprintf (string_params[index], "misc.asymmetry_penalty[1]");
    //sprintf (string_params[index], "misc.asymmetry_penalty[2]");
    //sprintf (string_params[index], "misc.asymmetry_penalty[3]");
    //sprintf (string_params[index], "misc.asymmetry_penalty[4]");
    //sprintf (string_params[index], "misc.asymmetry_penalty[5]");
    //sprintf (string_params[index], "misc.asymmetry_penalty[6]");

    //sprintf (string_params[index], "misc.gail_rule");
    switch (job)
    {
        case 0:
            sprintf (string_params[index], "misc.multi_offset");
            sprintf (string_params_human_readable[index], "multi_offset");
            break;
        case 2:
            array[index] = misc.multi_offset;
            break;
        case 3:
            fprintf (file, "%.2lf\n", (double)misc.multi_offset/100.0);
            break;
        case 4:
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf\n", &param);
            param *= 100;
            misc.multi_offset = (PARAMTYPE) param;
            break;
        case 5:
            misc.multi_offset = array[index];
            break;
    }
    index++;    sim_index++;
    switch (job)
    {
        case 0:
            sprintf (string_params[index], "misc.multi_helix_penalty");
            sprintf (string_params_human_readable[index], "multi_helix_penalty");
            break;
        case 2:
            array[index] = misc.multi_helix_penalty;
            break;
        case 3:
            fprintf (file, "%.2lf\n", (double)misc.multi_helix_penalty/100.0);
            break;
        case 4:
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf\n", &param);
            param *= 100;
            misc.multi_helix_penalty = (PARAMTYPE) param;
            break;
        case 5:
            misc.multi_helix_penalty = array[index];
            break;
    }
    index++;    sim_index++;
    switch (job)
    {
        case 0:
            sprintf (string_params[index], "misc.multi_free_base_penalty");
            sprintf (string_params_human_readable[index], "multi_free_base_penalty");
            break;
        case 2:
            array[index] = misc.multi_free_base_penalty;
            break;
        case 3:
            fprintf (file, "%.2lf\n", (double)misc.multi_free_base_penalty/100.0);
            break;
        case 4:
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf\n", &param);
            param *= 100;
            misc.multi_free_base_penalty = (PARAMTYPE) param;
            break;
        case 5:
            misc.multi_free_base_penalty = array[index];
            break;
    }
    index++;    sim_index++;
    switch (job)
    {
        case 0:
            sprintf (string_params[index], "misc.intermolecular_initiation");
            sprintf (string_params_human_readable[index], "intermolecular_initiation");
            break;
        case 2:
            array[index] = misc.intermolecular_initiation;
            break;
        case 3:
            fprintf (file, "%.2lf\n", (double)misc.intermolecular_initiation/100.0);
            break;
        case 4:
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf\n", &param);
            param *= 100;
            misc.intermolecular_initiation = (PARAMTYPE) param;
            break;
        case 5:
            misc.intermolecular_initiation = array[index];
            break;
    }
    index++;    sim_index++;

    if (parsi_special == T99)
    {
        for(i=0; i < nb_triloops; i++)
        {
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "triloop[%d].energy", i);
                    sprintf (string_params_human_readable[index], "triloop[%s]", triloop[i].seq);
                    break;
                case 2:
                    array[index] = triloop[i].energy;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)triloop[i].energy/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    triloop[i].energy = (PARAMTYPE) param;
                    break;
                case 5:
                    triloop[i].energy = array[index];
                    break;
            }
            index++;    sim_index++;
        }
        for(i=0; i < nb_tloops; i++)
        {
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "tloop[%d].energy", i);
                    sprintf (string_params_human_readable[index], "tetraloop[%s]", tloop[i].seq);
                    break;
                case 2:
                    array[index] = tloop[i].energy;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)tloop[i].energy/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    tloop[i].energy = (PARAMTYPE) param;
                    break;
                case 5:
                    tloop[i].energy = array[index];
                    break;
            }
            index++;    sim_index++;
        }
    }
    else if (parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        for (i=0; i < nb_special_hl; i++)
        {
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "special_hl[%d].energy", i);
                    sprintf (string_params_human_readable[index], "special_hairpin_loop[%s]", special_hl[i].seq);
                    break;
                case 2:
                    array[index] = special_hl[i].energy;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", (double)special_hl[i].energy/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    special_hl[i].energy = (PARAMTYPE) param;
                    break;
                case 5:
                    special_hl[i].energy = array[index];
                    break;
            }
            index++;    sim_index++;
        }
        // Now add the 6 internal special parameters
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_special_3GA");
                sprintf (string_params_human_readable[index], "internal_special_3GA");
                break;
            case 2:
                array[index] = misc.internal_special_3GA;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_special_3GA/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_special_3GA = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_special_3GA = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_special_2GA");
                sprintf (string_params_human_readable[index], "internal_special_2GA");
                break;
            case 2:
                array[index] = misc.internal_special_2GA;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_special_2GA/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_special_2GA = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_special_2GA = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_special_2xGA_GC");
                sprintf (string_params_human_readable[index], "internal_special_2xGA_GC");
                break;
            case 2:
                array[index] = misc.internal_special_2xGA_GC;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_special_2xGA_GC/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_special_2xGA_GC = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_special_2xGA_GC = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_special_midGA");
                sprintf (string_params_human_readable[index], "internal_special_midGA");
                break;
            case 2:
                array[index] = misc.internal_special_midGA;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_special_midGA/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_special_midGA = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_special_midGA = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_special_UG_AG");
                sprintf (string_params_human_readable[index], "internal_special_UG_AG");
                break;
            case 2:
                array[index] = misc.internal_special_UG_AG;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_special_UG_AG/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_special_UG_AG = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_special_UG_AG = array[index];
                break;
        }
        index++;    sim_index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal_special_GU_A");
                sprintf (string_params_human_readable[index], "internal_special_GU_A");
                break;
            case 2:
                array[index] = misc.internal_special_GU_A;
                break;
            case 3:
                fprintf (file, "%.2lf\n", (double)misc.internal_special_GU_A/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal_special_GU_A = (PARAMTYPE) param;
                break;
            case 5:
                misc.internal_special_GU_A = array[index];
                break;
        }
        index++;    sim_index++;
    }

    if ((job == 3) || (job == 4)) {
		fclose (file);
	}

    if ((job == 4) || (job == 5)) {
        extrapolate_parameters ();
    }

    return index;
}

int create_string_params ()
  // Mirela: Nov 23, 2003
  // writes each parameter type, excluding duplicates, in a long vector, containing the names of the parameters
  // Mirela: Dec 18, 2007: added a more human readable version of string_params
  // TODO: do I still need string_params? Or can I have one set of strings, that's human readable? Probably not, except that I use them in several places
{
    return traverse_features_and_do_work ("create_string_params", NULL, NULL);
}

void fill_data_structures_with_new_parameters (const char *filename)
  // Mirela: Dec 16, 2003
  // reads parameters from a file, and writes them in the internal data structures
  // PRE: first read the actual standard parameters, to be able to figure out which of them are
  // < INF, and maybe to also keep some old values.
{
    traverse_features_and_do_work ("fill_data_structures_with_new_parameters", NULL, filename);
}

void fill_data_structures_with_new_parameters_from_array (PARAMTYPE *array)
  // Mirela: 25 Aug 2008
  // reads parameters from and array, and writes them in the internal data structures
  // PRE: first read the actual standard parameters, to be able to figure out which of them are
  // < INF, and maybe to also keep some old values.
{
    traverse_features_and_do_work ("fill_data_structures_with_new_parameters_from_array", array, NULL);
}

int get_num_params ()
{
    num_params = create_string_params();
    return num_params;
}


int check_stability_and_size (int k, int l, int o, int p)
// helper function, to detect which delta we need for the int22 parameters
{
    // having at least one AC mismatch is the simplest, test first
    if ((k==A && l==C) || (k==C && l==A) || (o==A && p==C) || (o==C && p==A))
        return 4;

    // combination of all mismatches of equal size (purine-purine, purine-pyrimidine, and pyrimidine-pyrimidine are different sizes)
    // purine =  A, G
    // pyrimidine = C, U
    // if all purine-purines
    if ((k==A || k==G) && (l==A || l==G) && (o==A || o==G) && (p==A || p==G))
        return 1;
    // if all pyrimidine-pyrimidine
    if ((k==C || k==U) && (l==C || l==U) && (o==C || o==U) && (p==C || p==U))
        return 1;
    // if both  purine-pyrimidine
    // assume the A-C pairs have been found above
    if ( (((k==A || k==G) && (l==C || l==U)) || ((k==C || k==U) && (l==A || l==G))) &&
         (((o==A || o==G) && (p==C || p==U)) || ((o==C || o==U) && (p==A || p==G))) )
        return 1;
    // or any combination of 2 unstable mismatches except AC: AA, CC, CU, GG
    if ( ((k==A && l==A) || (k==C && l==C) || (k==C && l==U) || (k==U && l==C) || (k==G && l==G)) &&
         ((o==A && p==A) || (o==C && p==C) || (o==C && p==U) || (o==U && p==C) || (o==G && p==G)) )
        return 1;

    // two stabilizing mismatches (GU, GA, UU) of different sizes  (purine-purine, purine-pyrimidine, and pyrimidine-pyrimidine are different sizes)
    if ( (((k==G && l==U) || (k==U && l==G))   &&   ((o==G && p==A) || (o==A && p==G) || (o==U && p==U))) ||
         (((k==G && l==A) || (k==A && l==G))   &&   ((o==G && p==U) || (o==U && p==G) || (o==U && p==U))) ||
         ((k==U && l==U)                       &&   ((o==G && p==A) || (o==A && p==G) || (o==G && p==U) || (o==U && p==G))) )
        return 2;

    // one stable (GU, GA, UU) and one unstable mismatch (excluding AC) (AA, CC, CU, GG) of different sizes
    // GU
    if ( ((k==G && l==U) || (k==U && l==G)) &&
              ((o==A && p==A) || (o==C && p==C) || (o==C && p==U) || (o==U && p==C) || (o==G && p==G)) )
        return 3;
    if ( ((o==G && p==U) || (o==U && p==G)) &&
              ((k==A && l==A) || (k==C && l==C) || (k==C && l==U) || (k==U && l==C) || (k==G && l==G)) )
        return 3;
    // GA
    if ( ((k==G && l==A) || (k==A && l==G)) &&
              ((o==C && p==C) || (o==C && p==U) || (o==U && p==C)) )
        return 3;
    if ( ((o==G && p==A) || (o==A && p==G)) &&
              ((k==C && l==C) || (k==C && l==U) || (k==U && l==C)) )
        return 3;
    // UU
    if ( (k==U && l==U) &&
              ((o==A && p==A) || (o==G && p==G)) )
        return 3;
    if ( (o==U && p==U) &&
              ((k==A && l==A) || (k==G && l==G)) )
        return 3;

    return -1;

}


// OBSOLETE. Now it's part of traverse_features_and_do_work
/*
void fill_data_structures_with_new_parameters (const char *filename)
  // Mirela: Dec 16, 2003
  // reads parameters from a file, and writes them in the internal data structures
  // PRE: first read the actual standard parameters, to be able to figure out which of them are
  // < INF, and maybe to also keep some old values.
  // TODO: this should be incorporated into traverse_and_do_work
{
    #if (MODEL == EXTENDED)


    #elif (MODEL == SIMPLE)
    int index;
    int i, j, k, l, m, n, o, p;
    FILE *file;
    char buffer[100];
    double param;
    int line = 0;

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            int11[i][j][k][l][m][n] = INF;
                        }


    //printf ("FILENAME: %s\n", filename);
    if ((file = fopen (filename, "r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }

    index = 0;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
            {
                if (stack[i][j][k][l] < INF)
                {
                    // exclude duplicates
                    // stack[i][j][k][l] is the same as stack[l][k][j][i]
                    if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                    {
                        fgets (buffer, sizeof(buffer), file);
                        line++;
                        sscanf (buffer, "%lf\n", &param);
                        //printf ("\t%lf\n", param);
                        // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                        {
                            param *= 100;
                            stack[i][j][k][l] = (PARAMTYPE) param;
                            // add the duplicate too
                            stack[l][k][j][i] = (PARAMTYPE) param;
                        }
                        //sprintf (string_params[index++], "stack[%d][%d][%d][%d]", i, j, k, l);
                    }
                }
            }
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
            {
                if (tstackh[i][j][k][l] < INF)
                {
                    // no duplicates here
                    fgets (buffer, sizeof(buffer), file);
                    line++;
                    sscanf (buffer, "%lf", &param);
                    // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                    {
                        param *= 100;
                        tstackh[i][j][k][l] = (PARAMTYPE) param;
                    }
                    //sprintf (string_params[index++], "tstackh[%d][%d][%d][%d]", i, j, k, l);
                }
            }
    // use only 3 parameters, and fill the tstacki data structure by applying the Mathews 99 rules
    fgets (buffer, sizeof(buffer), file);
    line++;
    sscanf (buffer, "%lf", &param);
    param *= 100;
    misc.internal_AU_closure = (PARAMTYPE) param;

    fgets (buffer, sizeof(buffer), file);
    line++;
    sscanf (buffer, "%lf", &param);
    param *= 100;
    misc.internal_AG_mismatch = (PARAMTYPE) param;

    fgets (buffer, sizeof(buffer), file);
    line++;
    sscanf (buffer, "%lf", &param);
    param *= 100;
    misc.internal_UU_mismatch = (PARAMTYPE) param;

    // fill the tstacki data structure a bit later, after we read AU_penalty

//     for (i=0; i < NUCL; i++)
//         for (j=0; j < NUCL; j++)
//         for (k=0; k < NUCL; k++)
//             for (l=0; l < NUCL; l++)
//             {
//                 if (tstacki[i][j][k][l] < INF)
//                 {
//                     // no duplicates here
//                     fgets (buffer, sizeof(buffer), file);
//                     sscanf (buffer, "%lf", &param);
//                     // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
//                     {
//                         param *= 100;
//                         tstacki[i][j][k][l] = (PARAMTYPE) param;
//                     }
//                     //sprintf (string_params[index++], "tstacki[%d][%d][%d][%d]", i, j, k, l);
//                 }
//             }

    if (!simple_internal_energy)
    {
        // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {
                                if ( (((i==C && j==G) || (i==G && j==C)) &&
                                     ((m==C && n==G) || (m==G && n==C)) &&
                                     !can_pair(k,l)) ||
                                     (watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                {
                                    //if (int11[i][j][k][l][m][n] < INF)
                                    {
                                        // exclude duplicates
                                        // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                        if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                        {
                                            fgets (buffer, sizeof(buffer), file);
                                            line++;
                                            sscanf (buffer, "%lf", &param);
                                            param *= 100;
                                            int11[i][j][k][l][m][n] = (PARAMTYPE) param;
                                            // add the duplicate too
                                            int11[n][m][l][k][j][i] = (PARAMTYPE) param;
                                        }
                                    }
                                }
                            }
        fgets (buffer, sizeof(buffer), file);
        line++;
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal11_basic_mismatch = (PARAMTYPE) param;

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal11_GG_mismatch = (PARAMTYPE) param;


        // fill the int11 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        if ( ((i==C && j==G) || (i==G && j==C)) &&
                                             ((m==C && n==G) || (m==G && n==C)))
                                        {
                                            if (can_pair(k,l))
                                            {
                                                int11[i][j][k][l][m][n] = misc.internal11_basic_mismatch;
                                                // add the duplicate too
                                                //int11[n][m][l][k][j][i] = misc.internal11_basic_mismatch;
                                            }
                                        }

                                        else
                                        {
                                            if (!(watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                            {
                                                if (k==G && l==G)
                                                {
                                                    int11[i][j][k][l][m][n] = misc.internal11_GG_mismatch;
                                                    // add the duplicate too
                                                    //int11[n][m][l][k][j][i] = misc.internal11_GG_mismatch;
                                                }
                                                else
                                                    int11[i][j][k][l][m][n] = misc.internal11_basic_mismatch;
                                                if (has_AU_penalty(i,j))
                                                    int11[i][j][k][l][m][n] += misc.internal_AU_closure;
                                                if (has_AU_penalty(m,n))
                                                    int11[i][j][k][l][m][n] += misc.internal_AU_closure;
                                            }
                                        }

                                        // round it to match Turner parameters
                                        //if (int11[i][j][k][l][m][n] % 10 == 5) int11[i][j][k][l][m][n] += 5;
                                        //if (int11[i][j][k][l][m][n] % 10 == -5) int11[i][j][k][l][m][n] += 5;
                                    }
                }


//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//             for (k=0; k < NUCL; k++)
//                 for (l=0; l < NUCL; l++)
//                 for (m=0; m < NUCL; m++)
//                     for (n=0; n < NUCL; n++)
//                     {
//                         if (int11[i][j][k][l][m][n] < INF)
//                         {
//                             // exclude duplicates
//                             // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
//                             if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
//                             {
//                                 fgets (buffer, sizeof(buffer), file);
//                                 sscanf (buffer, "%lf", &param);
//                                 // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
//                                 {
//                                     param *= 100;
//                                     int11[i][j][k][l][m][n] = (PARAMTYPE) param;
//                                     // add the duplicate too
//                                     int11[n][m][l][k][j][i] = (PARAMTYPE) param;
//                                 }
//                                 //sprintf (string_params[index++], "int11[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
//                             }
//                         }
//                     }


        // go with few int21 parameters, as in Mathews et al 1999
        // closed by CG
        i=C; j=G; m=C; n=G;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                    if (!can_pair(k,l) && !can_pair(k,o))
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            fgets (buffer, sizeof(buffer), file);
                            sscanf (buffer, "%lf", &param);
                            param *= 100;
                            int21[i][j][k][l][m][n][o] = (PARAMTYPE) param;
                        }
        // closed by GC
        i=G; j=C; m=G; n=C;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                    if (!can_pair(k,l) && !can_pair(k,o))
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            fgets (buffer, sizeof(buffer), file);
                            sscanf (buffer, "%lf", &param);
                            param *= 100;
                            int21[i][j][k][l][m][n][o] = (PARAMTYPE) param;
                        }
        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal21_match = (PARAMTYPE) param;
        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal21_AU_closure = (PARAMTYPE) param;

        // fill the int21 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        for(o=0; o < NUCL; o++)
                                        {
                                            if ((i==C && j==G && m==C && n==G) ||  // these are already filled above, except what can pair inside
                                                (i==G && j==C && m==G && n==C))
                                            {
                                                if (can_pair(k,l) || can_pair(k,o))
                                                    int21[i][j][k][l][m][n][o] = misc.internal21_match;
                                            }
                                            else
                                            {
                                                if (can_pair(k,l) || can_pair(k,o))
                                                    int21[i][j][k][l][m][n][o] = misc.internal21_match;
                                                // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                                                else
                                                    int21[i][j][k][l][m][n][o] = (PARAMTYPE)(int21[C][G][k][l][C][G][o]/2.0) +
                                                    (PARAMTYPE)(int21[G][C][k][l][G][C][o]/2.0);
                                                if (has_AU_penalty(i,j))
                                                    int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                                if (has_AU_penalty(m,n))
                                                    int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                            }
                                            // round it to match Turner parameters - seems to be inconsistent
                                            //if (int21[i][j][k][l][m][n][o] % 10 == 5) int21[i][j][k][l][m][n][o] += 5;
                                        }
                                    }
                }

//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//             for (k=0; k < NUCL; k++)
//                 for (l=0; l < NUCL; l++)
//                 for (m=0; m < NUCL; m++)
//                     for (n=0; n < NUCL; n++)
//                     for(o=0; o < NUCL; o++)
//                         {
//                         if (int21[i][j][k][l][m][n][o] < INF)
//                             {
//                             // no duplicates here
//                             fgets (buffer, sizeof(buffer), file);
//                             sscanf (buffer, "%lf", &param);
//                             // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
//                                 {
//                                 param *= 100;
//                                 int21[i][j][k][l][m][n][o] = (PARAMTYPE) param;
//                                 }
//                             //sprintf (string_params[index++], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
//                             }
//                         }

        // go with the 53 parameters, like in Mathews et al 1999
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        n = i;
                        m = j;
                        p = k;
                        o = l;
                        if (watson_crick(i,j) && !watson_crick(k,l))
                        {
                            // exclude duplicates
                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                            {
                                fgets (buffer, sizeof(buffer), file);
                                sscanf (buffer, "%lf", &param);
                                // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                                {
                                    param *= 100;
                                    int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
                                    // add the duplicate too
                                    int22[n][m][p][o][j][i][l][k] = (PARAMTYPE) param;
                                }
                                //sprintf (string_params[index++], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                            }
                        }
                    }

//         i=C; j=G; m=C; n=G;
//         for (k=0; k < NUCL; k++)
//             for (l=0; l < NUCL; l++)
//                 for (o=0; o < NUCL; o++)
//                     for (p=0; p < NUCL; p++)
//                     {
//                         if (!watson_crick(k,l) && !watson_crick(o,p))
//                         {
//                             // exclude duplicates
//                             // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
//                             if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
//                                 n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
//                             {
//                                 fgets (buffer, sizeof(buffer), file);
//                                 sscanf (buffer, "%lf", &param);
//                                 // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
//                                 {
//                                     param *= 100;
//                                     int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
//                                     // add the duplicate too
//                                     int22[n][m][p][o][j][i][l][k] = (PARAMTYPE) param;
//                                 }
//                                 //sprintf (string_params[index++], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
//                             }
//                         }
//                     }

        // then add the 4 deltas
        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_same_size = (PARAMTYPE) param;

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_different_size = (PARAMTYPE) param;

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_1stable_1unstable = (PARAMTYPE) param;

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_AC = (PARAMTYPE) param;

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_match = (PARAMTYPE) param;

        int ii, jj, mm, nn;
        // fill the int22 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        for(o=0; o < NUCL; o++)
                                            for (p=0; p < NUCL; p++)
                                            {

//                                                 if (i==C && j==G && m==C && n==G)
//                                                 {
//                                                     if(watson_crick(k,l) || watson_crick(o,p))
//                                                     {
//                                                         int22[i][j][k][l][m][n][o][p] = misc.internal22_match;
//                                                     }
//                                                     // else do nothing, it's parameter
//                                                 }
                                                // if a closing pair is wobble, it's the same as if G would be A
                                                if (i==G && j==U)   ii = A;     else ii = i;
                                                if (i==U && j==G)   jj = A;     else jj = j;
                                                if (m==G && n==U)   mm = A;     else mm = m;
                                                if (m==U && n==G)   nn = A;     else nn = n;


                                                if (watson_crick(k,l) || watson_crick(o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] = misc.internal22_match;
                                                }
                                                else if (nn==ii && mm==jj && p==k && o==l)  // the UG closing pairs are the same as UA
                                                {
                                                    int22[i][j][k][l][m][n][o][p] = int22[ii][jj][k][l][mm][nn][o][p];
                                                }
                                                else //if (!(n==i && m==j && p==k && o==l))   // was already filled above
                                                {
                                                    int result = check_stability_and_size (k, l, o, p);
                                                    // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                                                    PARAMTYPE temp = (PARAMTYPE)(int22[ii][jj][k][l][jj][ii][l][k]/2.0) +
                                                               (PARAMTYPE)(int22[nn][mm][p][o][mm][nn][o][p]/2.0);
                                                    // rounf it to match Turner parameters
                                                    //if (temp%10 == 5) temp -= 5; if (temp%10 == -5) temp += 5;
                                                    int22[i][j][k][l][m][n][o][p] =  temp;
                                                    switch (result)
                                                    {
                                                        case 1: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_same_size; break;
                                                        case 2: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_different_size; break;
                                                        case 3: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_1stable_1unstable; break;
                                                        case 4: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_AC; break;
                                                        default: fprintf (stderr, "ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, k, l, o, p); exit(1);
                                                    }
                                                }
                                            }
                                    }
                }


//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//             for (k=0; k < NUCL; k++)
//                 for (l=0; l < NUCL; l++)
//                 for (m=0; m < NUCL; m++)
//                     for (n=0; n < NUCL; n++)
//                     for(o=0; o < NUCL; o++)
//                         for (p=0; p < NUCL; p++)
//                         {
//                             if (int22[i][j][k][l][m][n][o][p] < INF)
//                             {
//                                 // exclude duplicates
//                                 // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
//                                 if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
//                                     n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
//                                 {
//                                     fgets (buffer, sizeof(buffer), file);
//                                     sscanf (buffer, "%lf", &param);
//                                     // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
//                                     {
//                                         param *= 100;
//                                         int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
//                                         // add the duplicate too
//                                         int22[n][m][p][o][j][i][l][k] = (PARAMTYPE) param;
//                                     }
//                                     //sprintf (string_params[index++], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
//                                 }
//                             }
//                         }

    }     // end if (!simple_internal_energy)
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
            if (dangle_top[i][j][k] < INF)
                {
                // no duplicates here
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf", &param);
                // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                    {
                    param *= 100;
                    dangle_top[i][j][k] = (PARAMTYPE) param;
                    }
                //sprintf (string_params[index++], "dangle_top[%d][%d][%d]", i, j, k);
                }
            }
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        for (k=0; k < NUCL; k++)
            {
            if (dangle_bot[i][j][k] < INF)
                {
                // no duplicates here
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf", &param);
                // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                    {
                    param *= 100;
                    dangle_bot[i][j][k] = (PARAMTYPE) param;
                    }
                //sprintf (string_params[index++], "dangle_bot[%d][%d][%d]", i, j, k);
                }
            }
    int start;
    if (!simple_internal_energy)
        start = 4;
    else
        start = 1;
    for (i=start; i <= MAXLOOP_I; i++)
        {
        if (internal_penalty_by_size[i] < INF)
            {
            // no duplicates here
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf", &param);
            // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                {
                param *= 100;
                internal_penalty_by_size[i] = (PARAMTYPE) param;
                }
            //sprintf (string_params[index++], "internal_penalty_by_size[%d]", i);
            }
        }
    for (i=1; i <= MAXLOOP_B; i++)
        {
        if (bulge_penalty_by_size[i] < INF)
            {
            // no duplicates here
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf", &param);
            // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                {
                param *= 100;
                bulge_penalty_by_size[i] = (PARAMTYPE) param;
                }
            //sprintf (string_params[index++], "bulge_penalty_by_size[%d]", i);
            }
        }
    for (i=1; i <= MAXLOOP_H; i++)
        {
        if (hairpin_penalty_by_size[i] < INF)
            {
            // no duplicates here
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf", &param);
            // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                {
                param *= 100;
                hairpin_penalty_by_size[i] = (PARAMTYPE) param;
                }
            //sprintf (string_params[index++], "hairpin_penalty_by_size[%d]", i);
            }
        }
    //fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
    //    misc.param_greater30 = (PARAMTYPE) param;

    // set a fixed value to param_greater_30 for now
    misc.param_greater30 = 1.079;

    //sprintf (string_params[index++], "misc.param_greater30");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.terminal_AU_penalty = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.terminal_AU_penalty");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_GGG = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.hairpin_GGG");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_c1 = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.hairpin_c1");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_c2 = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.hairpin_c2");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_c3 = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.hairpin_c3");

    // fill the tstacki data structure, now that we also have terminal_AU_penalty - actually I removed it from here
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (!can_pair (i, j))
                        tstacki[i][j][k][l] = INF;
                    else
                    {
                        tstacki[i][j][k][l] = 0;
                        if (((i == A || i == G) && j == U) ||
                            ((j == A || j == G) && i == U))
                        {
                            tstacki[i][j][k][l] += misc.internal_AU_closure;
                        }
                        if ((k == A && l == G) ||
                            (l == A && k == G))
                        {
                            tstacki[i][j][k][l] += misc.internal_AG_mismatch;
                        }
                        if (k == U && l == U)
                        {
                            tstacki[i][j][k][l] += misc.internal_UU_mismatch;
                        }
                    }
                }


    // TODO
    // keep them fixed for now
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //fgets (buffer, sizeof(buffer), file);
    // sprintf (string_params[index++], "misc.asymmetry_penalty_max_correction");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[0]");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[1]");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[2]");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[3]");
    // Instead of these, I will just store the asymmetry for 0.5, 1, 1.5, 2, 2.5 and 3.

    // to come back!!!
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //sprintf (string_params[index++], "misc.asymmetry_penalty[1]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[2]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[3]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[4]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[5]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[6]");

    //fgets (buffer, sizeof(buffer), file);
    //sprintf (string_params[index++], "misc.gail_rule");
    // keep this fixed
    misc.gail_rule = 1;

    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.multi_offset = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.multi_offset");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.multi_helix_penalty = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.multi_helix_penalty");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.multi_free_base_penalty = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.multi_free_base_penalty");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.intermolecular_initiation = (PARAMTYPE) param;
    //sprintf (string_params[index++], "misc.intermolecular_initiation");

    for(i=0; i < nb_triloops; i++)
        {
        fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
        // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
            triloop[i].energy = (PARAMTYPE) param;
        //sprintf (string_params[index++], "triloop[%d].energy", i);
        }

    for(i=0; i < nb_tloops; i++)
        {
        fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
        // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
            tloop[i].energy = (PARAMTYPE) param;
        //sprintf (string_params[index++], "tloop[%d].energy", i);
        }

    fclose (file);
    //printf ("****** stack[1][2][0][3] = %d\n", stack[1][2][0][3]);
    #endif
}

*/




void save_paramtypes (const char *filename)
// PRE: call create_string_params ()
// save all parameter types in the given file
{
    int index;
    int i, j, k, l, m, n, o, p;
    FILE *file;

    if ((file = fopen (filename, "w")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    for (i=0; i < num_params; i++)
    {
        fprintf (file, "%s\n", string_params_human_readable[i]);
    }

    fclose (file);

}


void save_paramtypes_machine_readable (const char *filename)
// PRE: call create_string_params ()
// save all parameter types in the given file
{
    int index;
    int i, j, k, l, m, n, o, p;
    FILE *file;

    if ((file = fopen (filename, "w")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    for (i=0; i < num_params; i++)
    {
        fprintf (file, "%s\n", string_params[i]);
    }

    fclose (file);

}


void save_parameters (const char *filename)
  // Mirela: Dec 16, 2003
  // save all parameters in the given file
{
    traverse_features_and_do_work ("save_parameters", NULL, filename);
}


void save_parameters_in_array (PARAMTYPE *array)
// PRE: parameters have been read
// save all parameters in the given array
{
    traverse_features_and_do_work ("save_parameters_in_array", array, NULL);
}

double simfold_restricted_logZ (char *sequence, char *real_structure, char *restricted, double &min_energy, double &max_energy, int &actual_num_str)
// if the real_structure is not one of the suboptimal structures, add it
{
    char structure[MAXSLEN];
    double enthalpy, energy;
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];
    double Z;
    int i;
    double R, temp, beta;
    int real_str_found;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    real_str_found = 0;


    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted ();
    min_fold->return_structure (structure);
    delete min_fold;

    s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, -(int)(min_energy*100.0));
    sub_fold->set_limit (MAXSUBSTR);
    sub_fold->s_simfold_restricted (enthalpy);
    actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
    delete sub_fold;

    max_energy = tmp_energies[actual_num_str-1];
    Z = 0;
    for (i=0; i < actual_num_str; i++)
    {
        // check if the real structure is one of the suboptimal ones
        if (strcmp (tmp_structures[i], real_structure) == 0)
            real_str_found = 1;
        // recompute the free energy, i.e. with the correct dangling ends
        energy = free_energy_simfold_restricted (sequence, tmp_structures[i], restricted);
        Z += exp ((-1) * energy * beta);
    }
    if (!real_str_found)
    {
        // if not found, add it
        energy = free_energy_simfold_restricted (sequence, real_structure, restricted);
        Z += exp ((-1) * energy * beta);
    }
    return log(Z);
}


double simfold_restricted_logZ_gradient (char *sequence, char *real_structure, char *restricted, PFTYPE *logZ_gradient)
// return
// if the real_structure is not one of the suboptimal structures, add it
{
    char structure[MAXSLEN];
    double enthalpy, energy;
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];
    double min_energy, max_energy;
    int i, k;
    double R, temp, beta;
    int real_str_found, actual_num_str;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    real_str_found = 0;
    double numerator [MAXNUMPARAMS];
    double counter [MAXNUMPARAMS];
    double denominator;
    double f;

    s_min_folding *min_fold = new s_min_folding (sequence, restricted);
    min_energy = min_fold->s_simfold_restricted ();
    min_fold->return_structure (structure);
    delete min_fold;

    s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, -(int)(min_energy*100.0));
    sub_fold->set_limit (MAXSUBSTR);
    sub_fold->s_simfold_restricted (enthalpy);
    actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
    max_energy = tmp_energies[actual_num_str-1];
    delete sub_fold;

    for (i=0; i < num_params; i++)
    {
        numerator[i] = 0;
    }

    // first compute the denominator, which is the same for all parameters
    denominator = 0;
    for (k=0; k < actual_num_str; k++)
    {
        // check if the real structure is one of the suboptimal ones
        if (strcmp (tmp_structures[k], real_structure) == 0)
            real_str_found = 1;
        // recompute the free energy, i.e. with the correct dangling ends
        energy = free_energy_simfold_restricted (sequence, tmp_structures[k], restricted);
        denominator += exp ((-1) * energy * beta);
        count_each_structure_type (sequence, tmp_structures[k], restricted, counter, f, 1);
        for (i=0; i < num_params; i++)
        {
            numerator[i] += counter[i] * exp ((-1) * energy * beta);
        }
        //printf ("Sub str %d, counter[i] = %d, numerator[i] = %e\n", k, counter[7646], numerator[7646]);
    }
    if (!real_str_found)
    {
        //printf ("Real str not found, add it\n");
        // if not found, add it
        energy = free_energy_simfold_restricted (sequence, real_structure, restricted);
        denominator += exp ((-1) * energy * beta);
        count_each_structure_type (sequence, real_structure, restricted, counter, f, 1);
        for (i=0; i < num_params; i++)
        {
            numerator[i] += counter[i] * exp ((-1) * energy * beta);
        }
    }
    // now the denominator and nominator are computed

    for (i=0; i < num_params; i++)
    {
        logZ_gradient[i] = numerator[i] / denominator;
    }
    //printf ("logZ_gradient[i] = %e\n", logZ_gradient[7646]);
}


double compute_probability (double energy, double Z)
{
    double pb;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    pb = exp ((-1) * energy * beta) / Z;
    return pb;
}

int get_info_from_file (FILE *file, char *sequence, char *real_structure, char *restricted)
// return 1 if everything was ok
// return 0 instead of continue;
{
    int MLformat = 0;
    char buffer [5000];
    int i;

    if (MLformat)
    {
        if (fgets (buffer, sizeof(buffer), file) == NULL)
            return 0;
        if (strstr (buffer, "Seq: ") != NULL)
        {
            for (i=5; i < strlen (buffer)-1; i++)
                sequence[i-5] = buffer[i];
            sequence[i-5] = '\0';
            // ignore sequences longer than MAXSLEN -1
            if (strlen (sequence) > MAXSLEN-1) return 0;
        }
        else return 0;
        fgets (buffer, sizeof(buffer), file);
        if (strstr (buffer, "Str: ") != NULL)
        {
            for (i=5; i < strlen (buffer)-1; i++)
                real_structure[i-5] = buffer[i];
            real_structure[i-5] = '\0';
        }
        else
        {
            fprintf (stderr, "Str doesn't follow Seq, Seq is %s\n", sequence);
            exit(1);
        }
        fgets (buffer, sizeof(buffer), file);
        if (strstr (buffer, "Res: ") != NULL)
        {
            for (i=5; i < strlen (buffer)-1; i++)
                restricted[i-5] = buffer[i];
            restricted[i-5] = '\0';
        }
        else
        {
            fprintf (stderr, "Res doesn't follow Str\n");
            exit(1);
        }
        fgets (buffer, sizeof(buffer), file);   //====
    }
    else
    {
        if (fgets (buffer, sizeof(buffer), file) == NULL)
            return 0;
        if (buffer[0] == '>')
        {
            // next should be sequence
            fgets (buffer, sizeof(buffer), file);
            for (i=0; i < strlen (buffer)-1; i++)
                sequence[i] = buffer[i];
            sequence[i] = '\0';

            // next should be real structure
            fgets (buffer, sizeof(buffer), file);
            for (i=0; i < strlen (buffer)-1; i++)
                real_structure[i] = buffer[i];
            real_structure[i] = '\0';

            // next is restricted or empty line or predicted structure
            fgets (buffer, sizeof(buffer), file);
            if (buffer[0] == '\n')
            {
                restricted[0] = '\0';
                return 1;   // that's it, we are done
            }
            for (i=0; i < strlen (buffer)-1; i++)
                restricted[i] = buffer[i];
            restricted[i] = '\0';
            // the restricted string must contain _, otherwise it's the predicted structure
            if (strstr (restricted, "_") == NULL)   restricted[0] = '\0';
            fgets (buffer, sizeof(buffer), file);   // empty line
        }
        else return 0;
    }
    return 1;
}


PFTYPE compute_f (char *input_file)
// PRE:  the parameters have been read
// POST: computes the f function, which is sum_{i=1}^N (1/RT G_{i,nat,theta} + logZ)
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    FILE *file;
    PFTYPE logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    PFTYPE f;
    double min_energy, max_energy;
    int i, actual_num_str;

    f = 0;
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }

    int k;
    k = 0;
    int seen = 0;
    while (!feof (file))
    {
        seen++;
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;
        //printf ("      ....,....1....,....2....,....3....,....4....,....5\n");
        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
        energy = free_energy_simfold_restricted (sequence, real_structure, restricted);
        logZ = simfold_restricted_logZ (sequence, real_structure, restricted, min_energy, max_energy, actual_num_str);
        //printf ("logZ(%d) = %e\n", seen, logZ);
        f += beta*energy + logZ;
    }
    fclose (file);
    return f;
}

PFTYPE compute_likelihood_exactly (char *input_file)
// created on Oct 7, 2006
// PRE:  the parameters have been read
// POST: computes the likelihood function, which is \prod_{i=1}^N (exp (1/RT G_{i,nat,theta}) / Z)
// the partition function is computed exactly.
// the restriction string is not included so far
{

    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    FILE *file;
    PFTYPE Z;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    PFTYPE f;

    f = 1.0;
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }

    //long double pf;
    //pf = simfold_partition_function_exactly ("CAAAAGUCUGGGCUAAGCCCACUGAUGAGCCGCUGAAAUGCGGCGAAACUUUUG");
    //printf ("Exact part fun: %Le\n", pf);

    int k;
    k = 0;
    int seen = 0;
    while (!feof (file))
    {
        seen++;
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;
        //printf ("      ....,....1....,....2....,....3....,....4....,....5\n");
        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);


        // for now ignore restricted if it exists
        energy = free_energy_simfold (sequence, real_structure);

        Z = simfold_partition_function_smart (sequence);
        //printf ("Z = %Le\n", Z);
        //printf ("numerator = %e\n", exp (-1.0*beta*energy));
        f *= exp (-1.0*beta*energy) / Z;
    }
    fclose (file);
    return f;
}


PFTYPE compute_log_likelihood_smart (char *input_file)
// created on Oct 7, 2006
// PRE:  the parameters have been read
// POST: computes the likelihood function, which is \prod_{i=1}^N (exp (1/RT G_{i,nat,theta}) / Z)
// the partition function is computed exactly - with teh smart McCaskill's dynamic programming algorithm.
// the restriction string is not included so far
{

    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    FILE *file;
    PFTYPE Z;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double f;

    f = 0.0;
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }

    int k;
    k = 0;
    int seen = 0;
    while (!feof (file))
    {
        seen++;
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;
        //printf ("      ....,....1....,....2....,....3....,....4....,....5\n");
        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);


        // for now ignore restricted if it exists
        energy = free_energy_simfold (sequence, real_structure);

        Z = simfold_partition_function_smart (sequence);
        //printf ("Z = %g, en = %g, f = %g  %s\n", Z, energy, 1.0*beta*energy + log(Z), sequence);
        //printf ("numerator = %e\n", exp (-1.0*beta*energy));
        f += 1.0*beta*energy + log(Z);
    }
    fclose (file);
    return f;
}



void compute_gradient_f (char *input_file, PFTYPE *f_gradient)
// PRE:  the parameters have been read
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    char buffer [5000];
    FILE *file;
    double logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double min_energy, max_energy;
    int i, actual_num_str;
    double counter [MAXNUMPARAMS];
    PFTYPE logZ_gradient [MAXNUMPARAMS];
    double f;

    for (i=0; i < num_params; i++)
    {
        f_gradient[i] = 0;
    }
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }
    int k;
    k = 0;
    while (!feof (file))
    {
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;
        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
        count_each_structure_type (sequence, real_structure, restricted, counter, f, 1);
        simfold_restricted_logZ_gradient (sequence, real_structure, restricted, logZ_gradient);
        for (i = 0; i < num_params; i++)
        {
            f_gradient[i] += counter[i] - logZ_gradient[i];
        }
        //printf ("counter[i]=%d, logZ_gradient[i] = %e, f_gradient[i] = %e\n", counter[7621], logZ_gradient[7621], f_gradient[7621]);
    }
    fclose (file);
    for (i = 0; i < num_params; i++)
    {
        f_gradient[i] *= beta;
    }
}


void compute_gradient_f_smart (char *input_file, PFTYPE *f_gradient)
// PRE:  the parameters have been read
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))
// unrestricted for now
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    char buffer [5000];
    FILE *file;
    double logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double min_energy, max_energy;
    int i, actual_num_str;
    double counter [MAXNUMPARAMS];
    PFTYPE logZ_gradient [MAXNUMPARAMS];
    double f;

    for (i=0; i < num_params; i++)
    {
        f_gradient[i] = 0;
    }
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }
    int k;
    k = 0;
    while (!feof (file))
    {
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;

        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
        count_each_structure_type (sequence, real_structure, empty_char, counter, f, 1);
        simfold_gradient_smart (sequence, logZ_gradient);
        for (i = 0; i < num_params; i++)
        {
            f_gradient[i] += counter[i] - logZ_gradient[i];
        }
        //printf ("counter[i]=%d, logZ_gradient[i] = %e, f_gradient[i] = %e\n", counter[7621], logZ_gradient[7621], f_gradient[7621]);
    }
    fclose (file);
    for (i = 0; i < num_params; i++)
    {
        f_gradient[i] *= beta;
    }
}



PFTYPE compute_f_and_gradient_f_smart (char *input_file, PFTYPE *f_gradient)
// PRE:  the parameters have been read
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))
// returns f, which is minus log likelihood
// unrestricted for now
// DEPRICATED, another function exists in evaluate_function_smart_map
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    char buffer [5000];
    FILE *file;
    PFTYPE logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double min_energy, max_energy;
    int i, actual_num_str;
    double counter [MAXNUMPARAMS];
    double free_value;
    PFTYPE logZ_gradient [MAXNUMPARAMS];
    PFTYPE Z;
    PFTYPE neglogli;
    neglogli = 0.0;

    for (i=0; i < num_params; i++)
    {
        f_gradient[i] = 0;
    }
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }
    int k;
    k = 0;
    while (!feof (file))
    {
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;

        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
        count_each_structure_type (sequence, real_structure, empty_char, counter, free_value, 1);
        energy = free_energy_simfold (sequence, real_structure);
        Z = simfold_f_and_gradient_smart (sequence, NULL, logZ_gradient);
        neglogli += 1.0*beta*energy + log(Z);
        for (i = 0; i < num_params; i++)
        {
            f_gradient[i] += counter[i] - logZ_gradient[i];
        }
        //printf ("counter[i]=%d, logZ_gradient[i] = %e, f_gradient[i] = %e\n", counter[7621], logZ_gradient[7621], f_gradient[7621]);
    }
    fclose (file);
    for (i = 0; i < num_params; i++)
    {
        f_gradient[i] *= beta;
    }
    return neglogli;
}


PFTYPE compute_f_and_gradient_f (char *input_file, PFTYPE *f_gradient)
// PRE:  the parameters have been read
// POST: computes the f function, which is sum_{i=1}^N (1/RT G_{i,nat,theta} + logZ)
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char restricted[MAXSLEN];
    char buffer [5000];
    FILE *file;
    PFTYPE logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double min_energy, max_energy;
    int i, actual_num_str;
    double counter [MAXNUMPARAMS];
    double free_value;
    PFTYPE logZ_gradient;
    PFTYPE neglogli;
    neglogli = 0.0;
    int real_str_found;
    double numerator [MAXNUMPARAMS];

    double denominator;
    char structure[MAXSLEN];
    double enthalpy;
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];

    for (i=0; i < num_params; i++)
    {
        f_gradient[i] = 0;
    }
    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }
    int k;
    k = 0;

    int seen = 0;
    while (!feof (file))
    {
        seen++;
        // fixed bug: real_str_found was initialized with 0 outside of the while loop
        real_str_found = 0;
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;
        //printf ("Seq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
        count_each_structure_type (sequence, real_structure, restricted, counter, free_value, 1);

        // now fold sequence
        s_min_folding *min_fold = new s_min_folding (sequence, restricted);
        min_energy = min_fold->s_simfold_restricted ();
        min_fold->return_structure (structure);
        delete min_fold;

        s_sub_folding* sub_fold = new s_sub_folding(sequence, restricted, -(int)(min_energy*100.0));
        sub_fold->set_limit (MAXSUBSTR);
        sub_fold->s_simfold_restricted (enthalpy);
        actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
        max_energy = tmp_energies[actual_num_str-1];
        delete sub_fold;

        for (i=0; i < num_params; i++)
        {
            numerator[i] = 0;
        }

        // first compute the denominator, which is the same for all parameters
        denominator = 0;
        for (k=0; k < actual_num_str; k++)
        {
            // check if the real structure is one of the suboptimal ones
            if (strcmp (tmp_structures[k], real_structure) == 0)
                real_str_found = 1;
            // recompute the free energy, i.e. with the correct dangling ends
            energy = free_energy_simfold_restricted (sequence, tmp_structures[k], restricted);
            denominator += exp ((-1) * energy * beta);
            count_each_structure_type (sequence, tmp_structures[k], restricted, counter, free_value, 1);
            for (i=0; i < num_params; i++)
            {
                numerator[i] += counter[i] * exp ((-1) * energy * beta);
            }
        }
        count_each_structure_type (sequence, real_structure, restricted, counter, free_value, 1);
        energy = free_energy_simfold_restricted (sequence, real_structure, restricted);
        if (!real_str_found)
        {
            // if not found, add it
            denominator += exp ((-1) * energy * beta);
            for (i=0; i < num_params; i++)
            {
                numerator[i] += counter[i] * exp ((-1) * energy * beta);
            }
        }
        // now the denominator and nominator are computed

        neglogli += beta*energy + log(denominator);
        //printf ("logZ(%d) = %e\n", seen, log(denominator));
        for (i=0; i < num_params; i++)
        {
            logZ_gradient = numerator[i] / denominator;
            f_gradient[i] += counter[i] - logZ_gradient;
        }
        //printf ("counter[i]=%d, logZ_gradient[i] = %e, f_gradient[i] = %e\n", counter[7621], numerator[7621]/denominator, f_gradient[7621]);
        //printf ("logZ_grad[i] = %e\n", numerator[7621]/denominator);
    }
    fclose (file);
    for (i = 0; i < num_params; i++)
    {
        f_gradient[i] *= beta;
    }
    return neglogli;
}


void compute_counts_vector_LP (char *input_file, double *total_counter)
// PRE:  the parameters have been read
// POST: given p input sequence/structure pairs, compute sum_{j=1}^p (c_j, cbar_j)
//          where c_j is the counts for the known structures
//          and cbar_j is the counts for the predicted mfe structure, with the given set of parameters
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char pred_structure[MAXSLEN];
    char restricted[MAXSLEN];
    char buffer [5000];
    double counter [MAXNUMPARAMS];
    FILE *file;
    double logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double min_energy, max_energy;
    int i, actual_num_str;
    double f;

    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }
    int k;
    k = 0;

    for (i = 0; i < num_params; i++)
    {
        total_counter[i] = 0;
    }


    while (!feof (file))
    {
        if (!get_info_from_file (file, sequence, real_structure, restricted))
            continue;
        //printf ("===\nReal:\nSeq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
        count_each_structure_type (sequence, real_structure, restricted, counter, f, 1);
        for (i = 0; i < num_params; i++)  total_counter[i] += counter[i];
        // now predict the structure
        simfold_restricted (sequence, restricted, pred_structure);
        //printf ("===\nPred:\nSeq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, pred_structure, restricted);
        count_each_structure_type (sequence, pred_structure, restricted, counter, f, 1);
        for (i = 0; i < num_params; i++)  total_counter[i] -= counter[i];
    }
}

int compute_counts_matrix_LP_helper (FILE *file)
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char pred_structure[MAXSLEN];
    char restricted[MAXSLEN];
    double counter_real [MAXNUMPARAMS];
    double counter_pred [MAXNUMPARAMS];
    double f;

    if (!get_info_from_file (file, sequence, real_structure, restricted))
        return 0;
    //printf ("===\nReal:\nSeq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, real_structure, restricted);
    count_each_structure_type (sequence, real_structure, restricted, counter_real, f, 1);
    // now predict the structure
    simfold_restricted (sequence, restricted, pred_structure);
    //printf ("===\nPred:\nSeq: |%s|\nStr: |%s|\nRes: |%s|\n", sequence, pred_structure, restricted);
    count_each_structure_type (sequence, pred_structure, restricted, counter_pred, f, 1);
    // write the counters to the std output
    printf ("[");
    for (int i = 0; i < num_params-1; i++)
    {
        printf ("%.2lf, ", counter_real[i] - counter_pred[i]);
    }
    printf ("%.2lf]", counter_real[num_params-1] - counter_pred[num_params-1]);
    return 1;
}

void compute_counts_matrix_LP (char *input_file, int train_samples)
// PRE:  the parameters have been read
//       train_samples = # of training instances. -1 if all.
// POST: given p input sequence/structure pairs, compute and display the counts c_j - cbar_j
//          where c_j is the counts for the known structures
//          and cbar_j is the counts for the predicted mfe structure, with the given set of parameters
{
    char sequence[MAXSLEN];
    char real_structure[MAXSLEN];
    char pred_structure[MAXSLEN];
    char restricted[MAXSLEN];
    char buffer [5000];
    double counter_real [MAXNUMPARAMS];
    double counter_pred [MAXNUMPARAMS];
    FILE *file;
    double logZ;
    double energy;
    double R, temp, beta;
    R = 0.00198717;
    temp = 310.15;
    beta = 1/(R*temp);
    double min_energy, max_energy;
    int i, actual_num_str;

    if ((file = fopen (input_file, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", input_file);
        exit (0);
    }

    int k;
    k = 0;

    printf ("[");
    compute_counts_matrix_LP_helper (file);
    if (train_samples > 0)
    {
        for (i = 1; i < train_samples; i++)
        {
            printf (",\n ");
            if (! compute_counts_matrix_LP_helper (file))
            {
                fprintf (stderr, "Error in file\n");
                exit(1);
            }
        }
    }
    else
    {
        while (!feof (file))
        {
            printf (",\n ");
            if (! compute_counts_matrix_LP_helper (file))
            {
                fprintf (stderr, "Error in file\n");
                exit(1);
            }
        }
    }
    printf ("]\n");
    fclose (file);
}


void find_indeces_of_bbtypes (int &first, int &last, const char *bbtype, int num_params)
// PRE: the string_params are filled
// assumes the bb of type are consecutive
{
  int i;
  int met_first = 0;
  last = -1;
  for (i=0; i < num_params; i++)
    {
      if (strstr (string_params[i], bbtype) != NULL)
        {
          if (!met_first)
            {
              met_first = 1;
              first = i;
            }
        }
      else
        {
          if (met_first)  // this must be the next after last
            {
              last = i-1;
              break;
            }
        }
    }
  if (last == -1)
    last = num_params-1;
}


int choose_bbtype_randomly (int first, int last)
// PRE: the string_params are filled, the indeces first and last have been found
{
  int random_number;
  // generate a number between first and last
  random_number = first + (int) ((last-first+1.0)*rand()/(RAND_MAX+1.0));
  return random_number;
}



int generate_structure_withbb (char *sequence, char *known_structure, char *given_restricted, char *turner_structure, char structures[][MAXSLEN], double *old_counts, int threshold)
// written by Mirela in Sep, 2005
// PRE: the bbseq_left, bbstr_left, bbseq_right and bbstr_right are filled
//        i.e. call num_params = create_building_block_strings();
// First try to find if the building block bb_index can exist in any structure.
//   If found, restrict that part, and call simfold_restricted
//   Return the number of structures found
//   Look for at most threshold structures
//   num_params = total number of parameters
//   true_fe = the estimated free energy of the true structure.
//        We consider only those structures whose energies are <= true_fe
//   Update old_counts every time I consider a structure
{
    char *position_left, *position_right;
    char *position_left2, *position_right2;  // for tstacki
    char restricted[MAXSLEN];
    int seqlen, i, numstr;
    int bb_index;
    double counter_other[MAXNUMPARAMS];
    double energies[MAXSUBSTR];
    double known_fe;
    int considered;
    double f;

    seqlen = strlen(sequence);
    numstr = 0;
    bb_index = 0;

    if (given_restricted[0] == '\0')
        known_fe = free_energy_simfold (sequence, known_structure);
    else
        known_fe = free_energy_simfold_restricted (sequence, known_structure, given_restricted);
    //printf ("Known fe: %.2lf\n", known_fe);

    // first find bb_index, which is < threshold
    while (bb_index < num_params)
    {
        position_left = strstr (sequence, bbseq_left[bb_index]);
        if ((old_counts[bb_index] >= threshold) ||
            (strcmp (bbseq_left[bb_index], "") == 0) ||     // building block not implemented
            (strstr (string_params[bb_index], "dangle_") != NULL))    // not implemented yet
        {
            bb_index++;
            continue;
        }
        //printf ("Looking for parameter %s\n", string_params[bb_index]);
        considered = 0;
        while (position_left != NULL && !considered)    // include only one structure having some building block, no more
        //while (position_left != NULL && old_counts[bb_index] < threshold)
        {
            position_right = position_left + strlen (bbseq_left[bb_index])+TURN;  //?
            // check I'm not out of bounds
            if (position_right - sequence >= seqlen)
                break;
            position_right = strstr (position_right, bbseq_right[bb_index]);
            while (position_right != NULL && old_counts[bb_index] < threshold)
            {
                // found a new place for my building block
                for (i=0; i < seqlen; i++)   restricted[i] = '_'; restricted [seqlen] = '\0';
                replace_str_piece (restricted, position_left-sequence, bbstr_left[bb_index]);
                replace_str_piece (restricted, position_right-sequence, bbstr_right[bb_index]);
                // if it if type tstackh, then everything in between bbstr_left and bbstr_right must be "...."
                if (strstr (string_params[bb_index], "tstackh") != NULL)
                {
                    if (position_right - position_left > MAXLOOP)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    for (i=0; i < position_right - position_left - strlen (bbstr_left[bb_index]); i++)
                        restricted[position_left - sequence + strlen (bbstr_left[bb_index]) + i] = '.';
                }
                if (strstr (string_params[bb_index], "tstacki") != NULL)
                {
                    position_left2 = NULL;
                    position_right2 = NULL;
                    int found;
                    found = 0;
                    int bb_index2;
                    int tried = 0;
                    // we have to find another tstacki inside, at the right distance
                    // the second stacki can be any stacki, not only the current one
                    int first, last, k;
                    find_indeces_of_bbtypes (first, last, "tstacki", num_params);
                    // try about 20 times
                    for (k=0; k < 20; k++)
                    {
                        bb_index2 = choose_bbtype_randomly (first, last);
                        //printf ("first=%d, last=%d, chosen = %d\n", first, last, bb_index2);
                        if (strstr (string_params[bb_index2], "tstacki") == NULL)
                        {
                            fprintf (stderr, "bbtype was not chosen correctly\n");
                            exit(1);
                        }
                        //printf ("The second tstacki chosen was: %s\n", string_params[bb_index2]);
                        tried = 1;
                        //printf ("Looking for %s\n",  bbseq_right[bb_index2]);
                        position_left2 = strstr (position_left + strlen (bbseq_left[bb_index]), bbseq_right[bb_index2]);
                        while (position_left2 != NULL)
                        {
                            if (position_left2 > position_right -  strlen(bbseq_right[bb_index2]))   // is not inside
                                break;
                            position_right2 = strstr (position_left2 + strlen (bbseq_right[bb_index]) + TURN, bbseq_left[bb_index2]);
                            while (position_right2 != NULL)
                            {
                                if (position_right2 > position_right - strlen(bbseq_left[bb_index2]))
                                    break;
                                // now we should have both 1 and 2 somewhere
                                // now make sure the distance is right
                                if ((position_left2 - position_left == 2 && position_right - position_right2 == 2) ||
                                    (position_left2 - position_left + position_right - position_right2 > MAXLOOP))
                                    // not good, they shouldn't be 2 and 2 or longer than MAXLOOP
                                    // keep trying
                                {
                                    position_right2 = strstr (position_right2+1, bbseq_left[bb_index2]);
                                    continue;
                                }
                                // finally, found:
                                restricted[position_left2-sequence] = '.';
                                restricted[position_left2-sequence+1] = '(';
                                restricted[position_right2-sequence] = ')';
                                restricted[position_right2-sequence+1] = '.';
                                // now fill everything in between
                                for (i=position_left - sequence + strlen (bbseq_left[bb_index]); i < position_left2-sequence; i++)
                                    restricted[i] = '.';
                                for (i=position_right2 - sequence + strlen (bbseq_right[bb_index2]); i < position_right-sequence; i++)
                                    restricted[i] = '.';
                                found = 1;
                                break;
                            }  // end while position_right2
                            if (found) break;
                            position_left2 = strstr (position_left2+1, bbseq_right[bb_index2]);
                        }  // end while position_right1
                        if (found) break;
                    }  // end for bb_index2
                    if (!found)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    //printf ("\t%s\n", restricted);

                } // end if tstacki
                //printf ("\t%s\n", restricted);

                // now make sure restricted is compatible with given_restricted
                if (given_restricted[0] != '\0')
                {
                    //printf ("--%s\n++%s\n", given_restricted, restricted);
                    int comp = restricted_compatible (given_restricted, restricted);
                    if (comp)
                    {
                        //printf ("**%s\n\n", restricted);
                    }
                    else
                    {
                        //printf ("INCOMPATIBLE!\n\n");
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                }
                energies[numstr] = simfold_restricted (sequence, restricted, structures[numstr]);

                // test that this structure is not the same as the known_structure, turner_structure, or any other structure generated so far
                int to_consider = 1;
                if (strcmp (known_structure, structures[numstr]) == 0)          to_consider = 0;
                else if (strcmp (turner_structure, structures[numstr]) == 0)    to_consider = 0;
                else
                {
                    for (i=0; i < numstr; i++)
                        if (strcmp (structures[i], structures[numstr]) == 0)
                            to_consider = 0;
                }

                if (to_consider && energies[numstr] <= known_fe)    // only then consider it
                {
                    // update old_counts, by adding the new counts
                    count_each_structure_type (sequence, structures[numstr], restricted, old_counts, f, 0);
                    numstr++;
                    considered = 1;
                    if (numstr >= MAXSUBSTR)     // I have to make it stop somewhere, just in case, otherwise it can run out of bounds
                        return numstr;
                }
                position_right = strstr (position_right+1, bbseq_right[bb_index]);
            } // end while position_right
            position_left = strstr (position_left+1, bbseq_left[bb_index]);
        }    // end while position_left
        bb_index++;
    }
    return numstr;

  // stack, int11, int21, int22 are good to go - done
  // tstackh: everything in between bbstr_left and bbstr_right must be "...." - done
  // tstacki: in between left and right, there must be another tstacki, and in between, we must have "..." "..." of appropriate length - done
  // dangle_top and dangle_bot - must be in exterior loop (too complicated in multi loop) - NOT DONE YET
}



/*
int generate_structure_withoutbb (char *sequence, char *known_structure, char *given_restricted, char *turner_structure, char structures[][MAXSLEN], int *old_counts, int threshold)
// written by Mirela in Jan, 2006
// If known_structure has a "critical" building block, i.e. which does not appear in the thermodynamic set, then predict a few structures which don't have this building block
{
    char *position_left, *position_right;
    char *position_left2, *position_right2;  // for tstacki
    char restricted[MAXSLEN];
    int seqlen, i, numstr;
    int bb_index;
    int counter_known[MAXNUMPARAMS];
    int counter_other[MAXNUMPARAMS];
    double energies[MAXSUBSTR];

    double known_fe;
    int considered;

    // first see which of the critical building blocks are in the known structure
    count_each_structure_type (sequence, structure_pred, counter_known, 1);
    for (i=0; i < numparams; i++)
    {
        if (counter_known[i] > 0 && old_counts[i] < threshold)
        {
            position_left = strstr (sequence, bbseq_left[bb_index]);
        }
    }





    seqlen = strlen(sequence);
    numstr = 0;
    bb_index = 0;

    if (given_restricted[0] == '\0')
        known_fe = free_energy_simfold (sequence, known_structure);
    else
        known_fe = free_energy_simfold_restricted (sequence, known_structure, given_restricted);
    //printf ("Known fe: %.2lf\n", known_fe);

    // first find bb_index, which is < threshold
    while (bb_index < num_params)
    {
        position_left = strstr (sequence, bbseq_left[bb_index]);
        if ((old_counts[bb_index] >= threshold) ||
            (strcmp (bbseq_left[bb_index], "") == 0) ||     // building block not implemented
            (strstr (string_params[bb_index], "dangle_") != NULL))    // not implemented yet
        {
            bb_index++;
            continue;
        }
        //printf ("Looking for parameter %s\n", string_params[bb_index]);
        considered = 0;
        while (position_left != NULL && !considered)    // include only one structure having some building block, no more
        //while (position_left != NULL && old_counts[bb_index] < threshold)
        {
            position_right = position_left + strlen (bbseq_left[bb_index])+TURN;  //?
            // check I'm not out of bounds
            if (position_right - sequence >= seqlen)
                break;
            position_right = strstr (position_right, bbseq_right[bb_index]);
            while (position_right != NULL && old_counts[bb_index] < threshold)
            {
                // found a new place for my building block
                for (i=0; i < seqlen; i++)   restricted[i] = '_'; restricted [seqlen] = '\0';
                replace_str_piece (restricted, position_left-sequence, bbstr_left[bb_index]);
                replace_str_piece (restricted, position_right-sequence, bbstr_right[bb_index]);
                // if it if type tstackh, then everything in between bbstr_left and bbstr_right must be "...."
                if (strstr (string_params[bb_index], "tstackh") != NULL)
                {
                    if (position_right - position_left > MAXLOOP)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    for (i=0; i < position_right - position_left - strlen (bbstr_left[bb_index]); i++)
                        restricted[position_left - sequence + strlen (bbstr_left[bb_index]) + i] = '.';
                }
                if (strstr (string_params[bb_index], "tstacki") != NULL)
                {
                    position_left2 = NULL;
                    position_right2 = NULL;
                    int found;
                    found = 0;
                    int bb_index2;
                    int tried = 0;
                    // we have to find another tstacki inside, at the right distance
                    // the second stacki can be any stacki, not only the current one
                    int first, last, k;
                    find_indeces_of_bbtypes (first, last, "tstacki", num_params);
                    // try about 20 times
                    for (k=0; k < 20; k++)
                    {
                        bb_index2 = choose_bbtype_randomly (first, last);
                        //printf ("first=%d, last=%d, chosen = %d\n", first, last, bb_index2);
                        if (strstr (string_params[bb_index2], "tstacki") == NULL)
                        {
                            fprintf (stderr, "bbtype was not chosen correctly\n");
                            exit(1);
                        }
                        //printf ("The second tstacki chosen was: %s\n", string_params[bb_index2]);
                        tried = 1;
                        //printf ("Looking for %s\n",  bbseq_right[bb_index2]);
                        position_left2 = strstr (position_left + strlen (bbseq_left[bb_index]), bbseq_right[bb_index2]);
                        while (position_left2 != NULL)
                        {
                            if (position_left2 > position_right -  strlen(bbseq_right[bb_index2]))   // is not inside
                                break;
                            position_right2 = strstr (position_left2 + strlen (bbseq_right[bb_index]) + TURN, bbseq_left[bb_index2]);
                            while (position_right2 != NULL)
                            {
                                if (position_right2 > position_right - strlen(bbseq_left[bb_index2]))
                                    break;
                                // now we should have both 1 and 2 somewhere
                                // now make sure the distance is right
                                if ((position_left2 - position_left == 2 && position_right - position_right2 == 2) ||
                                    (position_left2 - position_left + position_right - position_right2 > MAXLOOP))
                                    // not good, they shouldn't be 2 and 2 or longer than MAXLOOP
                                    // keep trying
                                {
                                    position_right2 = strstr (position_right2+1, bbseq_left[bb_index2]);
                                    continue;
                                }
                                // finally, found:
                                restricted[position_left2-sequence] = '.';
                                restricted[position_left2-sequence+1] = '(';
                                restricted[position_right2-sequence] = ')';
                                restricted[position_right2-sequence+1] = '.';
                                // now fill everything in between
                                for (i=position_left - sequence + strlen (bbseq_left[bb_index]); i < position_left2-sequence; i++)
                                    restricted[i] = '.';
                                for (i=position_right2 - sequence + strlen (bbseq_right[bb_index2]); i < position_right-sequence; i++)
                                    restricted[i] = '.';
                                found = 1;
                                break;
                            }  // end while position_right2
                            if (found) break;
                            position_left2 = strstr (position_left2+1, bbseq_right[bb_index2]);
                        }  // end while position_right1
                        if (found) break;
                    }  // end for bb_index2
                    if (!found)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    //printf ("\t%s\n", restricted);

                } // end if tstacki
                //printf ("\t%s\n", restricted);

                // now make sure restricted is compatible with given_restricted
                if (given_restricted[0] != '\0')
                {
                    //printf ("--%s\n++%s\n", given_restricted, restricted);
                    int comp = restricted_compatible (given_restricted, restricted);
                    if (comp)
                    {
                        //printf ("**%s\n\n", restricted);
                    }
                    else
                    {
                        //printf ("INCOMPATIBLE!\n\n");
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                }
                energies[numstr] = simfold_restricted (sequence, restricted, structures[numstr]);

                // test that this structure is not the same as the known_structure, turner_structure, or any other structure generated so far
                int to_consider = 1;
                if (strcmp (known_structure, structures[numstr]) == 0)          to_consider = 0;
                else if (strcmp (turner_structure, structures[numstr]) == 0)    to_consider = 0;
                else
                {
                    for (i=0; i < numstr; i++)
                        if (strcmp (structures[i], structures[numstr]) == 0)
                            to_consider = 0;
                }

                if (to_consider && energies[numstr] <= known_fe)    // only then consider it
                {
                    // update old_counts, by adding the new counts
                    count_each_structure_type (sequence, structures[numstr], old_counts, 0);
                    numstr++;
                    considered = 1;
                    if (numstr >= MAXSUBSTR)     // I have to make it stop somewhere, just in case, otherwise it can run out of bounds
                        return numstr;
                }
                position_right = strstr (position_right+1, bbseq_right[bb_index]);
            } // end while position_right
            position_left = strstr (position_left+1, bbseq_left[bb_index]);
        }    // end while position_left
        bb_index++;
    }
    return numstr;

  // stack, int11, int21, int22 are good to go - done
  // tstackh: everything in between bbstr_left and bbstr_right must be "...." - done
  // tstacki: in between left and right, there must be another tstacki, and in between, we must have "..." "..." of appropriate length - done
  // dangle_top and dangle_bot - must be in exterior loop (too complicated in multi loop) - NOT DONE YET
}

*/

int restricted_compatible (char *given_restricted, char *restricted)
// assume given_restricted doesn't have dots
// return 1 if the two are compatible
// given_restricted is input variable
// restricted is input-output, and it will contain the joint restriction
{
    int given_p_table[MAXSLEN];
    int p_table[MAXSLEN];
    detect_original_pairs (given_restricted, given_p_table);
    detect_original_pairs (restricted, p_table);
    int len = strlen (restricted);
    int i;
    int compatible = 1;
    for (i=0; i < len; i++)
    {
        if (p_table[i] == -1)    // dot
        {
            if (given_p_table[i] >= 0)
            {
                compatible = 0;
                break;
            }
        }
        else if (p_table[i] >= 0)
        {
            if (given_p_table[i] >= 0 && given_p_table[i] != p_table[i])
            {
                compatible = 0;
                break;
            }
        }
    }
    if (compatible)
    {
        for (i=0; i < len; i++)
        {
            if (given_p_table[i] >= -1)
                restricted[i] = given_restricted[i];
        }
    }

    return compatible;
}


int generate_structure_withbb_many_thresholds (char *sequence, char structures[][MAXSLEN], double energies[], double *old_counts, int many_thresholds[], double true_fe, int num_params)
// written by Mirela in Sep, 2005
// PRE: the bbseq_left, bbstr_left, bbseq_right and bbstr_right are filled
//        i.e. call num_params = create_building_block_strings();
// First try to find if the building block bb_index can exist in any structure.
//   If found, restrict that part, and call simfold_restricted
//   Return the number of structures found
//   Look for at most threshold structures
//   num_params = total number of parameters
//   true_fe = the estimated free energy of the true structure.
//        We consider only those structures whose energies are <= true_fe
//   Update old_counts every time I consider a structure

// many_thresholds is a threshold per parameter (INF of no threshols wanted)
{
    char *position_left, *position_right;
    char *position_left2, *position_right2;  // for tstacki
    char restricted[MAXSLEN];
    int seqlen, i, numstr;
    int bb_index;
    double counter_other[MAXNUMPARAMS];
    double f;

    seqlen = strlen(sequence);
    numstr = 0;
    bb_index = 0;

    // first find bb_index, which is < threshold
    while (bb_index < num_params)
    {
        position_left = strstr (sequence, bbseq_left[bb_index]);
        if ((many_thresholds[bb_index] >= INF) ||
            (old_counts[bb_index] >= many_thresholds[bb_index]) ||
            (strcmp (bbseq_left[bb_index], "") == 0) ||     // building block not implemented
            (strstr (string_params[bb_index], "dangle_") != NULL))    // not implemented yet
        {
            bb_index++;
            continue;
        }
        while (position_left != NULL && old_counts[bb_index] < many_thresholds[bb_index])
        {
            position_right = position_left + strlen (bbseq_left[bb_index])+TURN;  //?
            // check I'm not out of bounds
            if (position_right - sequence >= seqlen)
                break;
            position_right = strstr (position_right, bbseq_right[bb_index]);
            while (position_right != NULL && old_counts[bb_index] < many_thresholds[bb_index])
            {
                // found a new place for my building block
                for (i=0; i < seqlen; i++)   restricted[i] = '_'; restricted [seqlen] = '\0';
                replace_str_piece (restricted, position_left-sequence, bbstr_left[bb_index]);
                replace_str_piece (restricted, position_right-sequence, bbstr_right[bb_index]);
                // if it if type tstackh, then everything in between bbstr_left and bbstr_right must be "...."
                if (strstr (string_params[bb_index], "tstackh") != NULL)
                {
                    if (position_right - position_left > MAXLOOP)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    for (i=0; i < position_right - position_left - strlen (bbstr_left[bb_index]); i++)
                        restricted[position_left - sequence + strlen (bbstr_left[bb_index]) + i] = '.';
                }
                if (strstr (string_params[bb_index], "tstacki") != NULL)
                {
                    position_left2 = NULL;
                    position_right2 = NULL;
                    int found;
                    found = 0;
                    int bb_index2;
                    int tried = 0;
                    // we have to find another tstacki inside, at the right distance
                    // the second stacki can be any stacki, not only the current one
                    int first, last, k;
                    find_indeces_of_bbtypes (first, last, "tstacki", num_params);
                    // try about 20 times
                    for (k=0; k < 20; k++)
                    {
                        bb_index2 = choose_bbtype_randomly (first, last);
                        //printf ("first=%d, last=%d, chosen = %d\n", first, last, bb_index2);
                        if (strstr (string_params[bb_index2], "tstacki") == NULL)
                        {
                            fprintf (stderr, "bbtype was not chosen correctly\n");
                            exit(1);
                        }
                        tried = 1;
                        position_left2 = strstr (position_left + strlen (bbseq_left[bb_index]), bbseq_right[bb_index2]);
                        while (position_left2 != NULL)
                        {
                            if (position_left2 > position_right -  strlen(bbseq_right[bb_index2]))   // is not inside
                                break;
                            position_right2 = strstr (position_left2 + strlen (bbseq_right[bb_index]) + TURN, bbseq_left[bb_index2]);
                            while (position_right2 != NULL)
                            {
                                if (position_right2 > position_right - strlen(bbseq_left[bb_index2]))
                                    break;
                                // now we should have both 1 and 2 somewhere
                                // now make sure the distance is right
                                if ((position_left2 - position_left == 2 && position_right - position_right2 == 2) ||
                                    (position_left2 - position_left + position_right - position_right2 > MAXLOOP))
                                    // not good, they shouldn't be 2 and 2 or longer than MAXLOOP
                                    // keep trying
                                {
                                    position_right2 = strstr (position_right2+1, bbseq_left[bb_index2]);
                                    continue;
                                }
                                // finally, found:
                                restricted[position_left2-sequence] = '.';
                                restricted[position_left2-sequence+1] = '(';
                                restricted[position_right2-sequence] = ')';
                                restricted[position_right2-sequence+1] = '.';
                                // now fill everything in between
                                for (i=position_left - sequence + strlen (bbseq_left[bb_index]); i < position_left2-sequence; i++)
                                    restricted[i] = '.';
                                for (i=position_right2 - sequence + strlen (bbseq_right[bb_index2]); i < position_right-sequence; i++)
                                    restricted[i] = '.';
                                found = 1;
                                break;
                            }  // end while position_right2
                            if (found) break;
                            position_left2 = strstr (position_left2+1, bbseq_right[bb_index]);
                        }  // end while position_right1
                        if (found) break;
                    }  // end for bb_index2
                    if (!found)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    //printf ("\t%s\n", restricted);

                } // end if tstacki
                //printf ("\t%s\n", restricted);
                energies[numstr] = simfold_restricted (sequence, restricted, structures[numstr]);
                if (energies[numstr] <= true_fe)    // only then consider it
                {
                    // update old_counts, by adding the new counts
                    count_each_structure_type (sequence, structures[numstr], restricted, old_counts, f, 0);
                    numstr++;
                    if (numstr >= MAXSUBSTR)     // I have to make it stop somewhere, just in case, otherwise it can run out of bounds
                        return numstr;
                }
                position_right = strstr (position_right+1, bbseq_right[bb_index]);
            } // end while position_right
            position_left = strstr (position_left+1, bbseq_left[bb_index]);
        }    // end while position_left
        bb_index++;
    }
    return numstr;

  // stack, int11, int21, int22 are good to go - done
  // tstackh: everything in between bbstr_left and bbstr_right must be "...." - done
  // tstacki: in between left and right, there must be another tstacki, and in between, we must have "..." "..." of appropriate length - done
  // dangle_top and dangle_bot - must be in exterior loop (too complicated in multi loop) - NOT DONE YET
}




void search_bb (char *sequence, double *old_counts, int threshold, int num_params)
// written by Mirela in Sep, 2005
// PRE: the bbseq_left, bbstr_left, bbseq_right and bbstr_right are filled
//        i.e. call num_params = create_building_block_strings();
// First try to find if the building block bb_index can exist in any structure.
//   If found, restrict that part, and call simfold_restricted
//   Return the number of structures found
//   Look for at most threshold structures
//   num_params = total number of parameters
//   true_fe = the estimated free energy of the true structure.
//        We consider only those structures whose energies are <= true_fe
//   Update old_counts every time I consider a structure
{
    char *position_left, *position_right;
    char *position_left2, *position_right2;  // for tstacki
    char restricted[MAXSLEN];
    int seqlen, i, numstr;
    int bb_index;
    double counter_other[MAXNUMPARAMS];

    seqlen = strlen(sequence);
    numstr = 0;
    bb_index = 0;

    // first find bb_index, which is < threshold
    while (bb_index < num_params)
    {
        position_left = strstr (sequence, bbseq_left[bb_index]);
        if ((strcmp (bbseq_left[bb_index], "") == 0) ||     // building block not implemented
            (strstr (string_params[bb_index], "dangle_") != NULL))    // not implemented yet
        {
            bb_index++;
            continue;
        }
        while (position_left != NULL)
        {
            position_right = position_left + strlen (bbseq_left[bb_index])+TURN;  //?
            // check I'm not out of bounds
            if (position_right - sequence >= seqlen)
                break;
            position_right = strstr (position_right, bbseq_right[bb_index]);
            while (position_right != NULL)
            {
                // found a new place for my building block
                for (i=0; i < seqlen; i++)   restricted[i] = '_'; restricted [seqlen] = '\0';
                replace_str_piece (restricted, position_left-sequence, bbstr_left[bb_index]);
                replace_str_piece (restricted, position_right-sequence, bbstr_right[bb_index]);
                // if it if type tstackh, then everything in between bbstr_left and bbstr_right must be "...."
                if (strstr (string_params[bb_index], "tstackh") != NULL)
                {
                    if (position_right - position_left > MAXLOOP)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    for (i=0; i < position_right - position_left - strlen (bbstr_left[bb_index]); i++)
                        restricted[position_left - sequence + strlen (bbstr_left[bb_index]) + i] = '.';
                }
                if (strstr (string_params[bb_index], "tstacki") != NULL)
                {
                    position_left2 = NULL;
                    position_right2 = NULL;
                    int found;
                    found = 0;
                    int bb_index2;
                    int tried = 0;
                    // we have to find another tstacki inside, at the right distance
                    // the second stacki can be any stacki, not only the current one
                    int first, last, k;
                    find_indeces_of_bbtypes (first, last, "tstacki", num_params);
                    // try about 20 times
                    for (k=0; k < 20; k++)
                    {
                        bb_index2 = choose_bbtype_randomly (first, last);
                        //printf ("first=%d, last=%d, chosen = %d\n", first, last, bb_index2);
                        if (strstr (string_params[bb_index2], "tstacki") == NULL)
                        {
                            fprintf (stderr, "bbtype was not chosen correctly\n");
                            exit(1);
                        }
                        tried = 1;
                        position_left2 = strstr (position_left + strlen (bbseq_left[bb_index]), bbseq_right[bb_index2]);
                        while (position_left2 != NULL)
                        {
                            if (position_left2 > position_right -  strlen(bbseq_right[bb_index2]))   // is not inside
                                break;
                            position_right2 = strstr (position_left2 + strlen (bbseq_right[bb_index]) + TURN, bbseq_left[bb_index2]);
                            while (position_right2 != NULL)
                            {
                                if (position_right2 > position_right - strlen(bbseq_left[bb_index2]))
                                    break;
                                // now we should have both 1 and 2 somewhere
                                // now make sure the distance is right
                                if ((position_left2 - position_left == 2 && position_right - position_right2 == 2) ||
                                    (position_left2 - position_left + position_right - position_right2 > MAXLOOP))
                                    // not good, they shouldn't be 2 and 2 or longer than MAXLOOP
                                    // keep trying
                                {
                                    position_right2 = strstr (position_right2+1, bbseq_left[bb_index2]);
                                    continue;
                                }
                                // finally, found:
                                restricted[position_left2-sequence] = '.';
                                restricted[position_left2-sequence+1] = '(';
                                restricted[position_right2-sequence] = ')';
                                restricted[position_right2-sequence+1] = '.';
                                // now fill everything in between
                                for (i=position_left - sequence + strlen (bbseq_left[bb_index]); i < position_left2-sequence; i++)
                                    restricted[i] = '.';
                                for (i=position_right2 - sequence + strlen (bbseq_right[bb_index2]); i < position_right-sequence; i++)
                                    restricted[i] = '.';
                                found = 1;
                                break;
                            }  // end while position_right2
                            if (found) break;
                            position_left2 = strstr (position_left2+1, bbseq_right[bb_index]);
                        }  // end while position_right1
                        if (found) break;
                    }  // end for bb_index2
                    if (!found)
                    {
                        position_right = strstr (position_right+1, bbseq_right[bb_index]);
                        continue;
                    }
                    //printf ("\t%s\n", restricted);

                } // end if tstacki
                old_counts[bb_index]++;
                position_right = strstr (position_right+1, bbseq_right[bb_index]);
            } // end while position_right
            position_left = strstr (position_left+1, bbseq_left[bb_index]);
        }    // end while position_left
        bb_index++;
    }

  // stack, int11, int21, int22 are good to go - done
  // tstackh: everything in between bbstr_left and bbstr_right must be "...." - done
  // tstacki: in between left and right, there must be another tstacki, and in between, we must have "..." "..." of appropriate length - done
  // dangle_top and dangle_bot - must be in exterior loop (too complicated in multi loop) - NOT DONE YET
}



/*
void fill_data_structures_with_new_parameters_fixed_dangles (const char *filename, char *dangfilename)
// reads all params from filename, except the dangling parameters, which reads from a different file
// assumes the fm363 model

  // Mirela: Ian 11, 2007
  // reads parameters from a file, and writes them in the internal data structures
  // PRE: first read the actual standard parameters, to be able to figure out which of them are
  // < INF, and maybe to also keep some old values.
{
    int index;
    int i, j, k, l, m, n, o, p;
    FILE *file;
    char buffer[100];
    double param;
    int line = 0;

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            int11[i][j][k][l][m][n] = INF;
                        }


    //printf ("FILENAME: %s\n", filename);
    if ((file = fopen (filename, "r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }

    index = 0;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
            {
                if (stack[i][j][k][l] < INF)
                {
                    // exclude duplicates
                    // stack[i][j][k][l] is the same as stack[l][k][j][i]
                    if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                    {
                        fgets (buffer, sizeof(buffer), file);
                        line++;
                        sscanf (buffer, "%lf\n", &param);
                        //printf ("\t%lf\n", param);
                        // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                        {
                            param *= 100;
                            stack[i][j][k][l] = (int)(round(param));
                            // add the duplicate too
                            stack[l][k][j][i] = (int)(round(param));
                        }
                        //sprintf (string_params[index++], "stack[%d][%d][%d][%d]", i, j, k, l);
                    }
                }
            }
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
            {
                if (tstackh[i][j][k][l] < INF)
                {
                    // no duplicates here
                    fgets (buffer, sizeof(buffer), file);
                    line++;
                    sscanf (buffer, "%lf", &param);
                    // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                    {
                        param *= 100;
                        tstackh[i][j][k][l] = (int)(round(param));
                    }
                    //sprintf (string_params[index++], "tstackh[%d][%d][%d][%d]", i, j, k, l);
                }
            }
    // use only 3 parameters, and fill the tstacki data structure by applying the Mathews 99 rules
    fgets (buffer, sizeof(buffer), file);
    line++;
    sscanf (buffer, "%lf", &param);
    param *= 100;
    misc.internal_AU_closure = (int)(round(param));

    fgets (buffer, sizeof(buffer), file);
    line++;
    sscanf (buffer, "%lf", &param);
    param *= 100;
    misc.internal_AG_mismatch = (int)(round(param));

    fgets (buffer, sizeof(buffer), file);
    line++;
    sscanf (buffer, "%lf", &param);
    param *= 100;
    misc.internal_UU_mismatch = (int)(round(param));

    // fill the tstacki data structure a bit later, after we read AU_penalty

    if (!simple_internal_energy)
    {
        // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {
                                if ( (((i==C && j==G) || (i==G && j==C)) &&
                                     ((m==C && n==G) || (m==G && n==C)) &&
                                     !can_pair(k,l)) ||
                                     (watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                {
                                    //if (int11[i][j][k][l][m][n] < INF)
                                    {
                                        // exclude duplicates
                                        // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                        if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                        {
                                            fgets (buffer, sizeof(buffer), file);
                                            line++;
                                            sscanf (buffer, "%lf", &param);
                                            param *= 100;
                                            int11[i][j][k][l][m][n] = (int)(round(param));
                                            // add the duplicate too
                                            int11[n][m][l][k][j][i] = (int)(round(param));
                                        }
                                    }
                                }
                            }
        fgets (buffer, sizeof(buffer), file);
        line++;
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal11_basic_mismatch = (int)(round(param));

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal11_GG_mismatch = (int)(round(param));


        // fill the int11 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        if ( ((i==C && j==G) || (i==G && j==C)) &&
                                             ((m==C && n==G) || (m==G && n==C)))
                                        {
                                            if (can_pair(k,l))
                                            {
                                                int11[i][j][k][l][m][n] = misc.internal11_basic_mismatch;
                                                // add the duplicate too
                                                //int11[n][m][l][k][j][i] = misc.internal11_basic_mismatch;
                                            }
                                        }

                                        else
                                        {
                                            if (!(watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                            {
                                                if (k==G && l==G)
                                                {
                                                    int11[i][j][k][l][m][n] = misc.internal11_GG_mismatch;
                                                    // add the duplicate too
                                                    //int11[n][m][l][k][j][i] = misc.internal11_GG_mismatch;
                                                }
                                                else
                                                    int11[i][j][k][l][m][n] = misc.internal11_basic_mismatch;
                                                if (has_AU_penalty(i,j))
                                                    int11[i][j][k][l][m][n] += misc.internal_AU_closure;
                                                if (has_AU_penalty(m,n))
                                                    int11[i][j][k][l][m][n] += misc.internal_AU_closure;
                                            }
                                        }

                                        // round it to match Turner parameters
                                        //if (int11[i][j][k][l][m][n] % 10 == 5) int11[i][j][k][l][m][n] += 5;
                                        //if (int11[i][j][k][l][m][n] % 10 == -5) int11[i][j][k][l][m][n] += 5;
                                    }
                }

        // go with few int21 parameters, as in Mathews et al 1999
        // closed by CG
        i=C; j=G; m=C; n=G;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                    if (!can_pair(k,l) && !can_pair(k,o))
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            fgets (buffer, sizeof(buffer), file);
                            sscanf (buffer, "%lf", &param);
                            param *= 100;
                            int21[i][j][k][l][m][n][o] = (int)(round(param));
                        }
        // closed by GC
        i=G; j=C; m=G; n=C;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                    if (!can_pair(k,l) && !can_pair(k,o))
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            fgets (buffer, sizeof(buffer), file);
                            sscanf (buffer, "%lf", &param);
                            param *= 100;
                            int21[i][j][k][l][m][n][o] = (int)(round(param));
                        }
        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal21_match = (int)(round(param));
        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal21_AU_closure = (int)(round(param));

        // fill the int21 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        for(o=0; o < NUCL; o++)
                                        {
                                            if ((i==C && j==G && m==C && n==G) ||  // these are already filled above, except what can pair inside
                                                (i==G && j==C && m==G && n==C))
                                            {
                                                if (can_pair(k,l) || can_pair(k,o))
                                                    int21[i][j][k][l][m][n][o] = misc.internal21_match;
                                            }
                                            else
                                            {
                                                if (can_pair(k,l) || can_pair(k,o))
                                                    int21[i][j][k][l][m][n][o] = misc.internal21_match;
                                                // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                                                else
                                                    int21[i][j][k][l][m][n][o] = (int)(int21[C][G][k][l][C][G][o]/2) +
                                                                                 (int)(int21[G][C][k][l][G][C][o]/2);
                                                if (has_AU_penalty(i,j))
                                                    int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                                if (has_AU_penalty(m,n))
                                                    int21[i][j][k][l][m][n][o] += misc.internal21_AU_closure;
                                            }
                                            // round it to match Turner parameters - seems to be inconsistent
                                            //if (int21[i][j][k][l][m][n][o] % 10 == 5) int21[i][j][k][l][m][n][o] += 5;
                                        }
                                    }
                }
        // go with the 53 parameters, like in Mathews et al 1999
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        n = i;
                        m = j;
                        p = k;
                        o = l;
                        if (watson_crick(i,j) && !watson_crick(k,l))
                        {
                            // exclude duplicates
                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                            {
                                fgets (buffer, sizeof(buffer), file);
                                sscanf (buffer, "%lf", &param);
                                // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                                {
                                    param *= 100;
                                    int22[i][j][k][l][m][n][o][p] = (int)(round(param));
                                    // add the duplicate too
                                    int22[n][m][p][o][j][i][l][k] = (int)(round(param));
                                }
                                //sprintf (string_params[index++], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                            }
                        }
                    }
        // then add the 4 deltas
        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_same_size = (int)(round(param));

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_different_size = (int)(round(param));

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_1stable_1unstable = (int)(round(param));

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_delta_AC = (int)(round(param));

        fgets (buffer, sizeof(buffer), file);
        sscanf (buffer, "%lf", &param);
        param *= 100;
        misc.internal22_match = (int)(round(param));

        int ii, jj, mm, nn;
        // fill the int22 data structure
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                if (can_pair(i,j))
                {
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                    if (can_pair(m,n))
                                    {
                                        for(o=0; o < NUCL; o++)
                                            for (p=0; p < NUCL; p++)
                                            {

//                                                 if (i==C && j==G && m==C && n==G)
//                                                 {
//                                                     if(watson_crick(k,l) || watson_crick(o,p))
//                                                     {
//                                                         int22[i][j][k][l][m][n][o][p] = misc.internal22_match;
//                                                     }
//                                                     // else do nothing, it's parameter
//                                                 }
                                                // if a closing pair is wobble, it's the same as if G would be A
                                                if (i==G && j==U)   ii = A;     else ii = i;
                                                if (i==U && j==G)   jj = A;     else jj = j;
                                                if (m==G && n==U)   mm = A;     else mm = m;
                                                if (m==U && n==G)   nn = A;     else nn = n;


                                                if (watson_crick(k,l) || watson_crick(o,p))
                                                {
                                                    int22[i][j][k][l][m][n][o][p] = misc.internal22_match;
                                                }
                                                else if (nn==ii && mm==jj && p==k && o==l)  // the UG closing pairs are the same as UA
                                                {
                                                    int22[i][j][k][l][m][n][o][p] = int22[ii][jj][k][l][mm][nn][o][p];
                                                }
                                                else //if (!(n==i && m==j && p==k && o==l))   // was already filled above
                                                {
                                                    int result = check_stability_and_size (k, l, o, p);
                                                    // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                                                    int temp = (int)(int22[ii][jj][k][l][jj][ii][l][k]/2.0) +
                                                               (int)(int22[nn][mm][p][o][mm][nn][o][p]/2.0);
                                                    // rounf it to match Turner parameters
                                                    //if (temp%10 == 5) temp -= 5; if (temp%10 == -5) temp += 5;
                                                    int22[i][j][k][l][m][n][o][p] =  temp;
                                                    switch (result)
                                                    {
                                                        case 1: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_same_size; break;
                                                        case 2: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_different_size; break;
                                                        case 3: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_1stable_1unstable; break;
                                                        case 4: int22[i][j][k][l][m][n][o][p] += misc.internal22_delta_AC; break;
                                                        default: printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, k, l, o, p); exit(1);
                                                    }
                                                }
                                            }
                                    }
                }

    }     // end if (!simple_internal_energy)

    // now get the dangling ends from a different file
    // **********************************************************************************************************
    FILE *dangfile;
    if ((dangfile = fopen (dangfilename, "r")) == NULL)
    {
        giveup ("Cannot open file", dangfilename);
    }

    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
            if (dangle_top[i][j][k] < INF)
                {
                // no duplicates here
                fgets (buffer, sizeof(buffer), dangfile);
                sscanf (buffer, "%lf", &param);
                // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                    {
                    param *= 100;
                    dangle_top[i][j][k] = (int)(round(param));
                    }
                //sprintf (string_params[index++], "dangle_top[%d][%d][%d]", i, j, k);
                }
            }
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
        for (k=0; k < NUCL; k++)
            {
            if (dangle_bot[i][j][k] < INF)
                {
                // no duplicates here
                fgets (buffer, sizeof(buffer), dangfile);
                sscanf (buffer, "%lf", &param);
                // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                    {
                    param *= 100;
                    dangle_bot[i][j][k] = (int)(round(param));
                    }
                //sprintf (string_params[index++], "dangle_bot[%d][%d][%d]", i, j, k);
                }
            }
    fclose(dangfile);
    // **********************************************************************************************************


    int start;
    if (!simple_internal_energy)
        start = 4;
    else
        start = 1;
    for (i=start; i <= MAXLOOP_I; i++)
        {
        if (internal_penalty_by_size[i] < INF)
            {
            // no duplicates here
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf", &param);
            // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                {
                param *= 100;
                internal_penalty_by_size[i] = (int)(round(param));
                }
            //sprintf (string_params[index++], "internal_penalty_by_size[%d]", i);
            }
        }
    for (i=1; i <= MAXLOOP_B; i++)
        {
        if (bulge_penalty_by_size[i] < INF)
            {
            // no duplicates here
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf", &param);
            // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                {
                param *= 100;
                bulge_penalty_by_size[i] = (int)(round(param));
                }
            //sprintf (string_params[index++], "bulge_penalty_by_size[%d]", i);
            }
        }
    for (i=1; i <= MAXLOOP_H; i++)
        {
        if (hairpin_penalty_by_size[i] < INF)
            {
            // no duplicates here
            fgets (buffer, sizeof(buffer), file);
            sscanf (buffer, "%lf", &param);
            // if (param != 0)    // put Turner's parameters if it's 0, but I can do this in the learn.pl file
                {
                param *= 100;
                hairpin_penalty_by_size[i] = (int)(round(param));
                }
            //sprintf (string_params[index++], "hairpin_penalty_by_size[%d]", i);
            }
        }
    //fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
    //    misc.param_greater30 = (int)(round(param));

    // set a fixed value to param_greater_30 for now
    misc.param_greater30 = 1.079;

    //sprintf (string_params[index++], "misc.param_greater30");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.terminal_AU_penalty = (int)(round(param));
    //sprintf (string_params[index++], "misc.terminal_AU_penalty");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_GGG = (int)(round(param));
    //sprintf (string_params[index++], "misc.hairpin_GGG");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_c1 = (int)(round(param));
    //sprintf (string_params[index++], "misc.hairpin_c1");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_c2 = (int)(round(param));
    //sprintf (string_params[index++], "misc.hairpin_c2");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.hairpin_c3 = (int)(round(param));
    //sprintf (string_params[index++], "misc.hairpin_c3");

    // fill the tstacki data structure, now that we also have terminal_AU_penalty - actually I removed it from here
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (!can_pair (i, j))
                        tstacki[i][j][k][l] = INF;
                    else
                    {
                        tstacki[i][j][k][l] = 0;
                        if (((i == A || i == G) && j == U) ||
                            ((j == A || j == G) && i == U))
                        {
                            tstacki[i][j][k][l] += misc.internal_AU_closure;
                        }
                        if ((k == A && l == G) ||
                            (l == A && k == G))
                        {
                            tstacki[i][j][k][l] += misc.internal_AG_mismatch;
                        }
                        if (k == U && l == U)
                        {
                            tstacki[i][j][k][l] += misc.internal_UU_mismatch;
                        }
                    }
                }


    // TODO
    // keep them fixed for now
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //fgets (buffer, sizeof(buffer), file);
    // sprintf (string_params[index++], "misc.asymmetry_penalty_max_correction");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[0]");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[1]");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[2]");
    // sprintf (string_params[index++], "misc.asymmetry_penalty_array[3]");
    // Instead of these, I will just store the asymmetry for 0.5, 1, 1.5, 2, 2.5 and 3.

    // to come back!!!
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //fgets (buffer, sizeof(buffer), file); fgets (buffer, sizeof(buffer), file);
    //sprintf (string_params[index++], "misc.asymmetry_penalty[1]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[2]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[3]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[4]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[5]");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[6]");

    //fgets (buffer, sizeof(buffer), file);
    //sprintf (string_params[index++], "misc.gail_rule");
    // keep this fixed
    misc.gail_rule = 1;

    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.multi_offset = (int)(round(param));
    //sprintf (string_params[index++], "misc.multi_offset");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.multi_helix_penalty = (int)(round(param));
    //sprintf (string_params[index++], "misc.multi_helix_penalty");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.multi_free_base_penalty = (int)(round(param));
    //sprintf (string_params[index++], "misc.multi_free_base_penalty");
    fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
    // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
        misc.intermolecular_initiation = (int)(round(param));
    //sprintf (string_params[index++], "misc.intermolecular_initiation");

    for(i=0; i < nb_triloops; i++)
        {
        fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
        // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
            triloop[i].energy = (int)(round(param));
        //sprintf (string_params[index++], "triloop[%d].energy", i);
        }

    for(i=0; i < nb_tloops; i++)
        {
        fgets (buffer, sizeof(buffer), file); sscanf (buffer, "%lf", &param); param *= 100;
        // if (param != 0)      // put Turner's parameters if it's 0, but I can do this in the learn.pl file
            tloop[i].energy = (int)(round(param));
        //sprintf (string_params[index++], "tloop[%d].energy", i);
        }

    fclose (file);
    //printf ("****** stack[1][2][0][3] = %d\n", stack[1][2][0][3]);
}
*/


void print_parameters_in_almost_mfold_format ()
// save the parameters in a format that is easy to parse by an external perl script, and which puts the parameters in mfold (or RNAstructure) format.
// written on Oct 23, 2007
{
    // first print the stacking parameters
    int i, j, ii, jj, ip, jp, iip, jjp;
    int oi[6] = {0,1,2,3,2,3};    // the order in which i and ip appear
    int oj[6] = {3,2,1,0,3,2};    // the order in which j and jp appear
    int index, indexp;

    printf ("#Parameters for stack.dat\n");
    for (i=0; i < NUCL; i++)
    {
        printf ("#%c*/**\n", int_to_nuc (i));
        for (ii=0; ii < NUCL; ii++)
        {
            for (j=0; j < NUCL; j++)
            {
                for (jj=0; jj < NUCL; jj++)
                {
                    if (stack[i][j][ii][jj] < INF)
                        printf ("% 3.2lf ", stack[i][j][ii][jj]/100.0);
                        // the space after % puts a space if the number is positive, so that positive and negative numbers are aligned nicely
                    else
                        printf ("   .   ");
                }
            }
            printf ("\n");
        }
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for tstackh.dat\n");
    for (i=0; i < NUCL; i++)
    {
        printf ("#%c*/**\n", int_to_nuc (i));
        for (ii=0; ii < NUCL; ii++)
        {
            for (j=0; j < NUCL; j++)
            {
                for (jj=0; jj < NUCL; jj++)
                {
                    if (tstackh[i][j][ii][jj] < INF)
                        printf ("% 3.2lf ", tstackh[i][j][ii][jj]/100.0);
                    else
                        printf ("   .   ");
                }
            }
            printf ("\n");
        }
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for tstacki.dat\n");
    for (i=0; i < NUCL; i++)
    {
        printf ("#%c*/**\n", int_to_nuc (i));
        for (ii=0; ii < NUCL; ii++)
        {
            for (j=0; j < NUCL; j++)
            {
                for (jj=0; jj < NUCL; jj++)
                {
                    if (tstacki[i][j][ii][jj] < INF)
                        printf ("% 3.2lf ", tstacki[i][j][ii][jj]/100.0);
                    else
                        printf ("   .   ");
                }
            }
            printf ("\n");
        }
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for dangle.dat\n");
    for (i=0; i < NUCL; i++)
    {
        printf ("#%c*/*\n", int_to_nuc (i));
        for (j=0; j < NUCL; j++)
        {
            for (ii=0; ii < NUCL; ii++)
            {
                if (dangle_top[i][j][ii] < INF)
                    printf ("% 3.2lf ", dangle_top[i][j][ii]/100.0);
                else
                    printf ("   .   ");

            }
        }
        printf ("\n");
    }
    for (i=0; i < NUCL; i++)
    {
        printf ("#%c/**\n", int_to_nuc (i));
        for (j=0; j < NUCL; j++)
        {
            for (ii=0; ii < NUCL; ii++)
            {
                if (dangle_bot[i][j][ii] < INF)
                    printf ("% 3.2lf ", dangle_bot[i][j][ii]/100.0);
                else
                    printf ("   .   ");

            }
        }
        printf ("\n");
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for int11.dat\n");
    for (index=0; index < 6; index++)
    {
        printf ("#%c**/%c**\n", int_to_nuc (oi[index]), int_to_nuc (oj[index]));
        for (ii=0; ii < NUCL; ii++)
        {
            for (indexp=0; indexp < 6; indexp++)
            {
                for (jj=0; jj < NUCL; jj++)
                {
                    if (int11[oi[index]][oj[index]][ii][jj][oi[indexp]][oj[indexp]] < INF)
                        printf ("% 3.2lf ", int11[oi[index]][oj[index]][ii][jj][oi[indexp]][oj[indexp]]/100.0);
                    else
                        printf ("   .   ");
                }
            }
            printf ("\n");
        }
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for int21.dat\n");
    for (index=0; index < 6; index++)
    {
        for (jjp=0; jjp < NUCL; jjp++)
        {
            printf ("#%c**/%c*%c*\n", int_to_nuc (oi[index]), int_to_nuc (oj[index]), int_to_nuc(jjp));
            for (ii=0; ii < NUCL; ii++)
            {
                for (indexp=0; indexp < 6; indexp++)
                {
                    for (jj=0; jj < NUCL; jj++)
                    {
                        if (int21[oi[index]][oj[index]][ii][jj][oi[indexp]][oj[indexp]][jjp] < INF)
                            printf ("% 3.2lf ", int21[oi[index]][oj[index]][ii][jj][oi[indexp]][oj[indexp]][jjp]/100.0);
                        else
                            printf ("   .   ");
                    }
                }
                printf ("\n");
            }
        }
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for int22.dat\n");
    for (index=0; index < 6; index++)
    {
        for (indexp=0; indexp < 6; indexp++)
        {
            printf ("#3'-%c**%c-5'/5'-%c**%c-3'\n", int_to_nuc (oi[index]), int_to_nuc (oi[indexp]), int_to_nuc (oj[index]), int_to_nuc(oj[indexp]));
            for (ii=0; ii < NUCL; ii++)
            {
                for (jj=0; jj < NUCL; jj++)
                {
                    for (iip=0; iip < NUCL; iip++)
                    {
                        for (jjp=0; jjp < NUCL; jjp++)
                        {
                            if (int22[oi[index]][oj[index]][ii][jj][oi[indexp]][oj[indexp]][iip][jjp] < INF)
                                printf ("% 3.2lf ", int22[oi[index]][oj[index]][ii][jj][oi[indexp]][oj[indexp]][iip][jjp]/100.0);
                            else
                                printf ("   .   ");
                        }
                    }
                    printf ("\n");
                }
            }
        }
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for loop.dat\n");
    // just make sure the internal_penalty_by_size of 1, 2 and 3 is INF
    internal_penalty_by_size[1] = INF;
    internal_penalty_by_size[2] = INF;
    internal_penalty_by_size[3] = INF;
    for (i=1; i <= MAXLOOP; i++)
    {
        printf ("%d", i);
        if (internal_penalty_by_size[i] < INF)
            printf ("%15.2lf", internal_penalty_by_size[i]/100.0);
        else
            printf ("            .  ");
        if (bulge_penalty_by_size[i] < INF)
            printf ("%15.2lf", bulge_penalty_by_size[i]/100.0);
        else
            printf ("            .  ");
        if (hairpin_penalty_by_size[i] < INF)
            printf ("%15.2lf", hairpin_penalty_by_size[i]/100.0);
        else
            printf ("            .  ");
        printf ("\n");
    }

    printf ("#-------------------------------\n");
    printf ("#Parameters for miscloop.dat\n");
    printf ("Miscellaneous free energy rules\n");
    printf ("-------------------------------\n\n");
    printf ("Extrapolation for large loops based on polymer theory\n");
    printf ("internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30)\n");
    printf ("--> \n");
    printf ("1.079\n\n");
    printf ("asymmetric internal loops: the ninio equation\n");
    printf ("the maximum correction\n");
    printf ("--> \n");
    printf ("%.1lf\n", misc.asymmetry_penalty_max_correction/100.0);
    printf ("\n");
    printf ("the f(m) array (see Ninio for details)\n");
    printf ("-->\n");
    printf ("%.1lf  %.1lf  %.1lf  %.1lf\n", misc.asymmetry_penalty_array[0]/100.0, misc.asymmetry_penalty_array[1]/100.0, misc.asymmetry_penalty_array[2]/100.0, misc.asymmetry_penalty_array[3]/100.0);
    printf ("\n");
    printf ("multibranched loops\n");
    printf ("  offset,  per nuc penalty,  helix penalty\n");
    printf ("--> \n");
    printf ("%.2lf    %.2lf    %.2lf\n", misc.multi_offset/100.0, misc.multi_free_base_penalty/100.0, misc.multi_helix_penalty/100.0);
    printf ("\n");
    // I don't use these, so I leave the original values
    printf ("efn2 multibranched loops\n");
    printf ("   offset,  per nuc penalty,  helix penalty\n");
    printf ("-->\n");
    printf (" 9.3         0.0              -0.9\n");
    printf ("\n");
    // I don't use this, so I leave it the way it was
    printf ("miscloop asym\n");
    printf ("-->\n");
    printf ("0.9\n");
    printf ("\n");
    printf ("multiloop strain\n");
    printf ("--> \n");
    printf ("3.1\n");
    printf ("\n");
    printf ("terminal AU penalty\n");
    printf ("--> \n");
    printf ("%.2lf\n", misc.terminal_AU_penalty/100.0);
    printf ("\n");
    printf ("bonus for GGG hairpin\n");
    printf ("--> \n");
    printf ("%.2lf\n", misc.hairpin_GGG/100.0);
    printf ("\n");
    printf ("c hairpin slope\n");
    printf ("--> \n");
    printf ("%.2lf\n", misc.hairpin_c1/100.0);
    printf ("\n");
    printf ("c hairpin intercept\n");
    printf ("--> \n");
    printf ("%.2lf\n", misc.hairpin_c2/100.0);
    printf ("\n");
    printf ("c hairpin of 3\n");
    printf ("--> \n");
    printf ("%.2lf\n", misc.hairpin_c3/100.0);
    printf ("\n");
    printf ("Intermolecular initiation free energy\n");
    printf ("--> \n");
    printf ("%.2lf\n", misc.intermolecular_initiation/100.0);
    printf ("\n");
    // I don't use this, so I leave it the way it was
    printf ("Bonus for Single C bulges adjacent to C\n");
    printf ("-->\n");
    printf ("-0.9\n");

    printf ("#-------------------------------\n");
    printf ("#Parameters for tloop.dat\n");
    for (i=0; i < nb_tloops; i++)
    {
        printf ("%s    %.2lf\n", tloop[i].seq, tloop[i].energy/100.0);
    }
    printf ("#-------------------------------\n");
    printf ("#Parameters for triloop.dat\n");
    for (i=0; i < nb_triloops; i++)
    {
        printf ("%s    %.2lf\n", triloop[i].seq, triloop[i].energy/100.0);
    }
}



/*---------------------------------------------------------------*/

void print_parameters_in_ViennaRNA_format ()
// writes the parameters in Vienna RNA format, on the screen
// this function is somewhat similar to the function write_parameter_file from the Vienna RNA library
// see format documentation at http://www.tbi.univie.ac.at/~ivo/RNA/RNAlib/Param-Files.html
// DO NOT ADD ENTHALPIES, because my method doesn't produce enthalpies, only free energies
// THE VALUES for non-standard bases or base pairs are the same as in the original Vienna RNA file
// In the Vienna RNA original file in the dangle sections, the base pair is upside down
// TABLE int11, int21, int22 - I think the headers of each 5x5 table are switched, i.e. GC instead of CG, although the values seem to be correct
// written on Nov 19, 2007
{
    int i, j, k, l, m, n;
    int base1[] = {C,G,G,U,A,U};
    int base2[] = {G,C,U,G,U,A};
    // the numbering of bases in the Vienna RNA format is ours + 1
    // i.e. A=1, C=2, G=3, U=4.
    int num_base_pairs = 6;
    const char *pnames[] = {"CG", "GC", "GU", "UG", "AU", "UA", " @"};
    char bnames[] = "@ACGU";


    printf ("## RNAfold parameter file\n");


    // FIRST BLOCK, stack energies
    printf ("\n# stack_energies\n");
    printf ("/*  CG    GC    GU    UG    AU    UA    @  */\n");
    for (i=0; i < num_base_pairs; i++)
    {
        printf (" ");
        for (j=0; j < num_base_pairs; j++)
        {
            printf ("%6d", (int)round(stack[base1[j]][base2[j]][base2[i]][base1[i]]));
            // each pair is from 5' to 3'
        }
        printf (" NST\n");
    }
    printf ("    NST   NST   NST   NST   NST   NST NST\n");


    // NEXT BLOCK is tstackh
    printf ("\n# mismatch_hairpin\n");
    for (i=0; i < num_base_pairs; i++)
    {
        printf ("     0     0     0     0     0\n");    //other base and other base as unpaired bases
        for (j=0; j < NUCL; j++)
        {
            // other base and A C G or U
            // for CG and GC, the Vienna RNA file give -90 and -70. Leave the same.
            if (i==0)        printf ("   -90");
            else if (i==1)   printf ("   -70");
            else             printf ("     0");
            for (k=0; k < NUCL; k++)
            {
                printf ("%6d", (int)round(tstackh[base1[i]][base2[i]][j][k]));
            }
            printf ("\n");
        }
    }
    for (i=0; i < NUCL+1; i++)    // the last block of zeros for nonstandard base pair
        printf ("     0     0     0     0     0\n");


    // NEXT BLOCK is tstacki
    printf ("\n# mismatch_interior\n");
    for (i=0; i < num_base_pairs; i++)
    {
        printf ("     0     0     0     0     0\n");    //other base and other base as unpaired bases
        for (j=0; j < NUCL; j++)
        {
            // other base and A C G or U
            printf ("     0");
            for (k=0; k < NUCL; k++)
            {
                printf ("%6d", (int)round(tstacki[base1[i]][base2[i]][j][k]));
            }
            printf ("\n");
        }
    }
    // values from the Vienna RNA default file, for non-standard base pair
    printf ("    90    90    90    90    90\n");
    printf ("    90    90    90    90   -20\n");
    printf ("    90    90    90    90    90\n");
    printf ("    90   -20    90    90    90\n");
    printf ("    90    90    90    90    20\n");


    // NEXT BLOCK is dangle5, i.e. dangle_bot
    printf ("\n# dangle5\n");
    printf ("/*  @     A     C     G     U   */\n");
    printf ("   INF   INF   INF   INF   INF\n");
    for (i=0; i < num_base_pairs; i++)
    {
        // other base and A C G or U
        printf ("   INF");
        for (j=0; j < NUCL; j++)
        {
            // this is read upside down - base2 first and then base1
            printf ("%6d", (int)round(dangle_bot[base2[i]][base1[i]][j]));
        }
        printf ("\n");
    }
    printf ("     0     0     0     0     0\n");


    // NEXT BLOCK is dangle3, i.e. dangle_top
    printf ("\n# dangle3\n");
    printf ("/*  @     A     C     G     U   */\n");
    printf ("   INF   INF   INF   INF   INF\n");
    for (i=0; i < num_base_pairs; i++)
    {
        // other base and A C G or U
        printf ("   INF");
        for (j=0; j < NUCL; j++)
        {
            // this is read upside down - base2 first and then base1
            printf ("%6d", (int)round(dangle_top[base2[i]][base1[i]][j]));
        }
        printf ("\n");
    }
    printf ("     0     0     0     0     0\n");


    // NEXT BLOCK is int11
    // don't print "no pair" entries for interior loop arrays
    // But first get the maximum values, to be used for the non-standard values
    int ns_CG = -INF;
    int ns_AU = -INF;
    for (i=0; i < num_base_pairs; i++)
    {
        for (j=0; j < num_base_pairs; j++)
        {
            for (k=1; k <= NUCL; k++)
            {
                for (l=1; l <= NUCL; l++)
                {
                    if (i==0 || i==1 || j==0 || j==1)
                    {
                        if ((int)round(int11[base1[i]][base2[i]][k-1][l-1][base2[j]][base1[j]]) > ns_CG)
                            ns_CG = (int)round(int11[base1[i]][base2[i]][k-1][l-1][base2[j]][base1[j]]);
                    }
                    else
                    {
                        if ((int)round(int11[base1[i]][base2[i]][k-1][l-1][base2[j]][base1[j]]) > ns_AU)
                            ns_AU = (int)round(int11[base1[i]][base2[i]][k-1][l-1][base2[j]][base1[j]]);
                    }
                }
            }
        }
    }

    //printf ("ns_AU = %d\n", ns_AU);
    //printf ("ns_CG = %d\n", ns_CG);

    printf ("\n# int11_energies\n");
    for (i=0; i <= num_base_pairs; i++)
    {
        for (j=0; j <= num_base_pairs; j++)
        {
            printf ("/* %2s..%2s */\n", pnames[i], pnames[j]);
            // looks like in the original file these headers are switched, although the values seem to be correct

            // for the non-standard base pairs or bases, it looks like in the original file there's 110
            //    if at least one of the base pairs is C-G, and 170 otherwise.
            // I think the rule is to put the maximum instead, computed above.
            for (k=0; k <= NUCL; k++)
            {
                for (l=0; l <= NUCL; l++)
                {
                    if (i==num_base_pairs || j==num_base_pairs || k==0 || l==0)
                    {
                        if (i==0 || i==1 || j==0 || j==1)
                            printf ("%6d", ns_CG);
                        else
                            printf ("%6d", ns_AU);
                    }
                    else
                        printf ("%6d", (int)round(int11[base1[i]][base2[i]][k-1][l-1][base2[j]][base1[j]]));
                }
                printf ("\n");
            }

        }
    }

    // NEXT BLOCK is int21
    // first get the maximum values, to be used for the non-standard values
    int ns = -INF;
    for (i=0; i < num_base_pairs; i++)
    {
        for (j=0; j < num_base_pairs; j++)
        {
            for (k=1; k <= NUCL; k++)
            {
                for (l=1; l <= NUCL; l++)
                {
                    for (m=1; m <= NUCL; m++)
                    {
                        if ((int)round(int21[base1[i]][base2[i]][k-1][m-1][base2[j]][base1[j]][l-1]) > ns)
                            ns = (int)round(int21[base1[i]][base2[i]][k-1][m-1][base2[j]][base1[j]][l-1]);
                    }
                }
            }
        }
    }
    //printf ("ns=%d\n", ns);
    printf  ("\n# int21_energies\n");
    for (i=0; i <= num_base_pairs; i++)
    {
        for (j=0; j <= num_base_pairs; j++)
        {
            for (k=0; k <= NUCL; k++)
            {
                printf ("/* %2s.%c..%2s */\n", pnames[i], bnames[k], pnames[j]);
                for (l=0; l <= NUCL; l++)
                {
                    for (m=0; m <= NUCL; m++)
                    {
                        // if any of the bases or base pairs is non-standard, put ns. The original value was 550
                        if (i==num_base_pairs || j==num_base_pairs || k==0 || l==0 || m==0)
                            printf ("%6d", ns);
                        else
                            printf ("%6d", (int)round(int21[base1[i]][base2[i]][k-1][m-1][base2[j]][base1[j]][l-1]));
                    }
                    printf ("\n");
                }
            }
        }
    }


    // NEXT BLOCK is int22
    // first get the maximum value
    ns = -INF;
    for (i=0; i < num_base_pairs; i++)
    {
        for (j=0; j < num_base_pairs; j++)
        {
            for (k=0; k < NUCL; k++)
            {
                for (l=0; l < NUCL; l++)
                {
                    for (m=0; m < NUCL; m++)
                    {
                        for (n=0; n < NUCL; n++)
                        {
                            if ((int)round(int22[base1[i]][base2[i]][k][n][base2[j]][base1[j]][l][m]) > ns)
                                ns = (int)round(int22[base1[i]][base2[i]][k][n][base2[j]][base1[j]][l][m]);
                        }
                    }
                }
            }
        }
    }
    //printf ("ns=%d\n", ns);
    printf ("\n# int22_energies\n");
    for (i=0; i <= num_base_pairs; i++)
    {
        for (j=0; j <= num_base_pairs; j++)
        {
            // for bases, we don't list the non-standard bases for this case
            for (k=0; k < NUCL; k++)
            {
                for (l=0; l < NUCL; l++)
                {
                    printf ("/* %2s.%c%c..%2s */\n", pnames[i], bnames[k+1], bnames[l+1], pnames[j]);
                    for (m=0; m < NUCL; m++)
                    {
                        for (n=0; n < NUCL; n++)
                        {
                            // if any of the base pairs is non-standard, put the max value. Value in default.par is 340.
                            if (i==num_base_pairs || j==num_base_pairs)
                                printf ("%6d", ns);
                            else
                                printf ("%6d", (int)round(int22[base1[i]][base2[i]][k][n][base2[j]][base1[j]][l][m]));
                        }
                        printf ("\n");
                    }
                }
            }
        }
    }


    // NEXT BLOCK is penalty by size for hairpin loops
    printf ("\n# hairpin\n");
    printf ("   INF");
    for (i=1; i <= 30; i++)
    {
        int x = (int)round(penalty_by_size(i,'H'));
        if (x >= INF)    printf ("   INF");
        else             printf ("%6d", x);
        if (i%10 == 9)  printf ("\n");
    }
    printf ("\n");


    // NEXT BLOCK is penalty by size for bulge loops
    printf ("\n# bulge\n");
    printf ("   INF");
    for (i=1; i <= 30; i++)
    {
        int x = (int)round(penalty_by_size(i,'B'));
        if (x >= INF)    printf ("   INF");
        else             printf ("%6d", x);
        if (i%10 == 9)  printf ("\n");
    }
    printf ("\n");


    // NEXT BLOCK is penalty by size for internal loops
    printf ("\n# internal_loop\n");
    printf ("   INF");
    for (i=1; i <= 30; i++)
    {
        int x = (int)round(penalty_by_size(i,'I'));
        if (x >= INF)    printf ("   INF");
        else             printf ("%6d", x);
        if (i%10 == 9)  printf ("\n");
    }
    printf ("\n");


    // NEXT BLOCK is other misc energies
    printf ("\n# ML_params\n");
    printf ("/* F = cu*n_unpaired + cc + ci*loop_degree (+TermAU) */\n");
    printf ("/*\t    cu\t    cc\t    ci\t TerminalAU */\n");
    printf ("\t%6d\t%6d\t%6d\t%6d\n", (int)round(misc.multi_free_base_penalty), (int)round(misc.multi_offset),
        (int)round(misc.multi_helix_penalty), (int)round(misc.terminal_AU_penalty));

    printf ("\n# NINIO\n");
    printf ("/* Ninio = MIN(max, m*|n1-n2| */\n");
    printf ("/*       m   max              */\n");
    printf ("\t%3d %4d\n", (int)round(misc.asymmetry_penalty_array[2]), (int)round(misc.asymmetry_penalty_max_correction));

    printf ("\n# Tetraloops\n");
    for (i=0; i < nb_tloops; i++)
        printf ("\t%.6s\t%4d\n", tloop[i].seq, (int)round(tloop[i].energy));

    printf ("\n# Triloops\n");
    for (i=0; i < nb_triloops; i++)
        printf ("\t%.5s\t%4d\n", triloop[i].seq, (int)round(triloop[i].energy));

    printf ("\n#END\n");
}

// functions to read from the thermodynamic set XML file
int get_data_from_buffer (char *buffer, const char *header, const char last_char, char *output)
// function to get the sequence, structure etc data from the XML lines
{
    char *begin;
    begin  = strstr (buffer, header);
        //printf ("%s\n", buffer);
    if (begin == NULL) { return 0; }     //printf ("Formatting error, %s not there\n", header); exit(1); }
    begin += strlen(header);
    int i = 0;
    while (1)
    {
        if (begin[i] == last_char) { output[i] = '\0'; break; }
        output[i] = begin[i++];
    }
    if (strstr (header, "sequence") != NULL || strstr (header, "structure") != NULL)
    {
        for (i=0; i < strlen(output); i++)
        {
            output[i] = toupper(output[i]);
        }
    }
    return 1;
}


void fill_similarity_rule_with_optical_melting_reference (char *xml_filename)
// PRE: create_string_params must be called
// reads data from xml_filename, and fills up the array similarity_rule with the experiment id.
// started on Mar 18, 2008.
{
    FILE *xml;
    int i, j;
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    double counter_min[MAXNUMPARAMS];
    // sometimes the initial state is not the completely unfolded sequence. In this case, sequence0 and structure0 are the initial state
    char sequence0[MAXSLEN];
    char structure0[MAXSLEN];
    double counter0[MAXNUMPARAMS];    // just in case we have sequence0 and structure0
    char buffer [10000];
    char *begin;
    char exp_id[100];
    double f;

    // first reset the similarity_rule array and counter0

    for (i=0; i < num_params; i++)
    {
        counter0[i] = 0;
        similarity_rule[i][0] = '\0';
    }

    if ((xml = fopen (xml_filename, "r")) == NULL)
    {
        fprintf (stderr, "Cannot open file %s\n", xml_filename);
        exit (0);
    }

    int k;
    k = 0;
    int delta_index = 1;
    double margin;
    int counter = 0;

    while (!feof (xml))
    {
        // the file is in XML format, among whose lines it has the following types:
        //    <EXPERIMENT sequence="CCGG CCGG" structure="(((( ))))" dG-37="-4.36" STD-dG-37="0.1"> </EXPERIMENT>
        //    <EXPERIMENT sequence="UGACCUCA UGAGGUCA" structure="(((((((( ))))))))" dG-37="-12.34" ERR-dG-37="3%"> </EXPERIMENT>

        // try including only the first paper. With everything I can't save the solution because it says "No basic solution"
        //    Problem solved: In the cplex script, I have to write: "write file.vec"
        //if (counter >= 94) break;

        fgets (buffer, sizeof(buffer), xml);
        //printf ("%s", buffer);

        sequence[0] = '\0'; structure[0] = '\0';
        begin = strstr (buffer, "<EXPERIMENT");
        if (begin == NULL) continue;

        sequence0[0] = '\0';
        structure0[0] = '\0';
        get_data_from_buffer (buffer, "id=\"", '\"', exp_id);
        get_data_from_buffer (buffer, "sequence0=\"", '\"', sequence0);
        get_data_from_buffer (buffer, "structure0=\"", '\"', structure0);
        get_data_from_buffer (buffer, "sequence=\"", '\"', sequence);
        get_data_from_buffer (buffer, "structure=\"", '\"', structure);

        if (strcmp (sequence0, "") != 0)
        {
            count_each_structure_type (sequence0, structure0, empty_char, counter0, f, 1);
        }
        //printf ("Sequence:  %s\n", sequence);
        //printf ("Structure: %s\n", structure);
        count_each_structure_type (sequence, structure, empty_char, counter_min, f, 1);
        //printf ("DONE %s\n", exp_id);

        for (i=0; i < num_params; i++)
        {
            if (counter_min[i] - counter0[i] != 0)
            {
                if (similarity_rule[i][0] == '\0')
                    strcpy (similarity_rule[i], exp_id);
            }
        }
    }
    fclose (xml);
}


void fill_similarity_rules ()
{
    traverse_features_and_do_work ("fill_similarity_rules", NULL, NULL);
}

int is_special_internal_1 (int *sequence, int i, int j, int ip, int jp)
{
    int branch1 = ip-i-1;
    int branch2 = j-jp-1;
    if (branch1 < 3 || branch2 < 3)     return 0;
    // first misc.internal_special_3GA       5'-YGGA/GAAR-3' or 5'-GGAR/YGAA-3', in loops 3x3 and larger
    if ((isY (sequence[i]) && sequence[i+1] == G && sequence[i+2] == G && sequence[i+3] == A &&
         isR (sequence[j]) && sequence[j-1] == A && sequence[j-2] == A && sequence[j-3] == G) ||
         (isR (sequence[ip]) && sequence[ip-1] == A && sequence[ip-2] == G && sequence[ip-3] == G &&
         isY (sequence[jp]) && sequence[jp+1] == G && sequence[jp+2] == A && sequence[jp+3] == A))
        return 1;
    return 0;
}

int is_special_internal_2 (int *sequence, int i, int j, int ip, int jp)
{
    // next internal_special_2GA
    // 5'-GA/GA-3' next to a closing base pair, or 5'-GG/AA-3' next to a closing base pair, for 3x3, 3x4, 4x4, and 4x5 loops;
    //  ALSO 5'-RGGA/GAAY-3' or 5'-GGAY/YGAA-3' for 3x5, 3x6 and 4x6 loops.
    // internal_2GA is USED only if internal_3GA was not used

    int branch1 = ip-i-1;
    int branch2 = j-jp-1;
    if (branch1 < 3 || branch2 < 3)     return 0;

    if (is_special_internal_1 (sequence, i, j, ip, jp))     return 0;

    if ((branch1 == 3 && (branch2==3 || branch2==4)) ||
         (branch1 == 4 && (branch2==3 || branch2==4 || branch2==5)) ||
         (branch1 == 5 && branch2 == 4))
    {
        if ((sequence[i+1]==G && sequence[i+2]==A && sequence[j-1]==A && sequence[j-2]==G) ||
             (sequence[ip-2]==G && sequence[ip-1]==A && sequence[jp+2]==A && sequence[jp+1]==G) ||
             (sequence[i+1]==G && sequence[i+2]==G && sequence[j-1]==A && sequence[j-2]==A) ||
             (sequence[ip-2]==A && sequence[ip-1]==A && sequence[jp+2]==G && sequence[jp+1]==G))
        {
            return 1;
        }
    }
    if ((branch1==3 && (branch2==5 || branch2==6)) ||
         (branch2==3 && (branch1==5 || branch1==6)) ||
         (branch1 == 4 && branch2 == 6) || (branch2 == 4 && branch1 == 6))
    {
        if ((isR (sequence[i]) && sequence[i+1] == G && sequence[i+2] == G && sequence[i+3] == A &&
             isY (sequence[j]) && sequence[j-1] == A && sequence[j-2] == A && sequence[j-3] == G) ||
             (isY (sequence[ip]) && sequence[ip-1] == A && sequence[ip-2] == G && sequence[ip-3] == G &&
             isR (sequence[jp]) && sequence[jp+1] == G && sequence[jp+2] == A && sequence[jp+3] == A))
        {
            return 1;
        }
    }
    return 0;
}


int is_special_internal_3 (int *sequence, int i, int j, int ip, int jp)
{
    int branch1 = ip-i-1;
    int branch2 = j-jp-1;
    //next internal_special_2xGA_GC;         5'-GANGC/GANGC-3' in 3x3 loops
    if (branch1==3 && branch2==3)
    {
        if (sequence[i]==G && sequence[i+1]==A && sequence[ip-1]==G && sequence[ip]==C &&
            sequence[jp]==G && sequence[jp+1]==A && sequence[j-1]==G && sequence[j]==C)
            return 1;
    }
    return 0;
}

int is_special_internal_4 (int *sequence, int i, int j, int ip, int jp)
{
    int branch1 = ip-i-1;
    int branch2 = j-jp-1;

    // next internal_special_midGA       // middle GA adjacent to RY in 3x3 loops
                                    // internal_midGA is USED only if none of internal_3GA and internal_2GA is used.
    if (is_special_internal_1 (sequence, i, j, ip, jp))     return 0;
    if (is_special_internal_2 (sequence, i, j, ip, jp))     return 0;

    if (branch1==3 && branch2==3)
    {
        if (sequence[i+2]==G && sequence[j-2]==A &&
            ((isR(sequence[i+1]) && isY(sequence[j-1])) || (isR(sequence[ip-1]) && isY(sequence[jp+1]))))
            return 1;
        if (sequence[i+2]==A && sequence[j-2]==G &&
            ((isY(sequence[i+1]) && isR(sequence[j-1])) || (isY(sequence[ip-1]) && isR(sequence[jp+1]))))
            return 1;
    }
    return 0;
}

int is_special_internal_5 (int *sequence, int i, int j, int ip, int jp)
// this returns how many times this motif appears. Can be 0, once or twice.
{
    int branch1 = ip-i-1;
    int branch2 = j-jp-1;
    int result = 0;
    // next internal_special_UG_AG;       // once or twice, for each 5'-UG/AG-3' at the terminus of loops 3x3 or larger
    if (branch1==3 && branch2==3)
    {
        if (sequence[i]==U && sequence[j]==G && sequence[i+1]==G && sequence[j-1]==A)
            result++;
        if (sequence[jp]==U && sequence[ip]==G && sequence[jp+1]==G && sequence[ip-1]==A)
            result++;
    }
    return result;
}

int is_special_internal_6 (int *sequence, int i, int j, int ip, int jp)
{
    int branch1 = ip-i-1;
    int branch2 = j-jp-1;

    // next internal_special_GU_A;        // first mismatch is GA, and U is 3' of G, for loops 3x3
    if (branch1==3 && branch2==3)
    {
        if ((sequence[i+1]==G && sequence[j-1]==A && sequence[i+2] == U) ||
             (sequence[jp+1]==G && sequence[ip-1]==A && sequence[jp+2] == U))
            return 1;
    }
    return 0;
}


PARAMTYPE special_energy_internal (int *sequence, int i, int j, int ip, int jp)
// Return the energy obtained when we consider 6 additional parameters for internal loop 3x3 and larger,
//  as described in Chen_Turner_2006b.
// the arguments are positions in sequence
{
    PARAMTYPE energy = 0;
    if (parsi_special != LAVISH && parsi_special != T99_LAVISH)  return 0;

    if (is_special_internal_1 (sequence, i, j, ip, jp))
        energy += misc.internal_special_3GA;

    if (is_special_internal_2 (sequence, i, j, ip, jp))
        energy += misc.internal_special_2GA;

    if (is_special_internal_3 (sequence, i, j, ip, jp))
        energy += misc.internal_special_2xGA_GC;

    if (is_special_internal_4 (sequence, i, j, ip, jp))
        energy += misc.internal_special_midGA;

    energy += misc.internal_special_UG_AG * is_special_internal_5 (sequence, i, j, ip, jp);

    if (is_special_internal_6 (sequence, i, j, ip, jp))
        energy += misc.internal_special_GU_A;
    return energy;
}

PARAMTYPE count_special_internal (double *counter, int *sequence, int i, int j, int ip, int jp)
// Return the energy and counts obtained when we consider 6 additional parameters for internal loop 3x3 and larger,
//  as described in Chen_Turner_2006b.
// the arguments are positions in sequence
{
    PARAMTYPE energy = 0;
    if (parsi_special != LAVISH && parsi_special != T99_LAVISH)  return 0;

    int index;
    if (is_special_internal_1 (sequence, i, j, ip, jp))
    {
        energy += misc.internal_special_3GA;
        index = structure_type_index ("misc.internal_special_3GA");
        counter[index]++;
    }

    if (is_special_internal_2 (sequence, i, j, ip, jp))
    {
        energy += misc.internal_special_2GA;
        index = structure_type_index ("misc.internal_special_2GA");
        counter[index]++;
    }

    if (is_special_internal_3 (sequence, i, j, ip, jp))
    {
        energy += misc.internal_special_2xGA_GC;
        index = structure_type_index ("misc.internal_special_2xGA_GC");
        counter[index]++;
    }

    if (is_special_internal_4 (sequence, i, j, ip, jp))
    {
        energy += misc.internal_special_midGA;
        index = structure_type_index ("misc.internal_special_midGA");
        counter[index]++;
    }

    int sp5 = is_special_internal_5 (sequence, i, j, ip, jp);
    if (sp5 > 0)
    {
        energy += misc.internal_special_UG_AG * sp5;
        index = structure_type_index ("misc.internal_special_UG_AG");
        counter[index] += sp5;
    }

    if (is_special_internal_6 (sequence, i, j, ip, jp))
    {
        energy += misc.internal_special_GU_A;
        index = structure_type_index ("misc.internal_special_GU_A");
        counter[index]++;
    }
    return energy;
}



/////////////////////////////////////////////////////////////////////////////////////
// Pieces of traverse_features_and_do_work that I tried and later abandoned.


//     for (i=1; i < MAXLOOP_I; i++)
//     {
//         for (j=MAX(i,3); j < MAXLOOP_I; j++)
//         {
//             if (internal_penalty_by_size_2D[i][j]  < INF)
//             {
//                 // no duplicates here
//                 switch (job)
//                 {
//                     case 0:
//                         sprintf (string_params[index], "internal_penalty_by_size_2D[%d][%d]", i, j);
//                         sprintf (string_params_human_readable[index], "internal_size_2D[%d,%d]", i, j);
//                         break;
//                     case 1:
//                         if (similarity_rule[sim_index][0] == '\0')
//                         {
//                             // the idea is: similarity (i,j) = JS(x,y), where |i-j| >+ |x-y| first, and i+j >+ x+y
//                             //  (where by >+ I mean greater, but as close as possible)
//                             if (i==j)       // x=i-1, y=j-1
//                             {
//                                 // apply Jackobson_Stockmeyer on i-1, j-1
//                                 sprintf (similarity_rule[sim_index], "1 * internal_size_2D[%d,%d] + %.4lf * extrapolation_large_loops",
//                                     i-1, j-1, log(1.0*(i+j)/(i-1+j-1)));
//                             }
//                             else    // x=i, y=j-1
//                             {
//                                 // apply Jackobson_Stockmeyer on i-1, j-1
//                                 sprintf (similarity_rule[sim_index], "1 * internal_size_2D[%d,%d] + %.4lf * extrapolation_large_loops",
//                                     i, j-1, log(1.0*(i+j)/(i+j-1)));
//                             }
//
//                             // THIS I TRIED TO DO BEFORE, but I CHANGED MY MIND
//                             //sprintf (similarity_rule[sim_index], "1 * internal_symmetry[%d]", i-1);
//                             //internal_penalty_by_size_2D[i][j] = internal_penalty_by_size[i+j] + asymmetry_penalty(i,j);
//                         }
//                         break;
//                     case 2:
//                         array[index] = internal_penalty_by_size_2D[i][j];
//                         break;
//                     case 3:
//                         fprintf (file, "%.2lf\n", internal_penalty_by_size_2D[i][j]/100.0);
//                         break;
//                 }
//                 index++;
//
//             }
//         }
//     }

//     for (i=3; i <= MAXLOOP_I/2; i++)
//     {
//         if (internal_symmetry[i] < INF)
//         {
//             // no duplicates here
//             switch (job)
//             {
//                 case 0:
//                     sprintf (string_params[index], "internal_symmetry[%d]", i);
//                     sprintf (string_params_human_readable[index], "internal_symmetry[%d]", i);
//                     break;
//                 case 1:
//                     if (similarity_rule[sim_index][0] == '\0')
//                     {
//                         sprintf (similarity_rule[sim_index], "1 * internal_symmetry[%d]", i-1);
//                     }
//                     break;
//                 case 2:
//                     array[index] = internal_symmetry[i];
//                     break;
//                 case 3:
//                     fprintf (file, "%.2lf\n", internal_symmetry[i]/100.0);
//                     break;
//             }
//             index++;
//         }
//     }
//     for (i=1; i <= MAXLOOP_I-2; i++)
//     {
//         if (internal_asymmetry[i] < INF)
//         {
//             // no duplicates here
//             switch (job)
//             {
//                 case 0:
//                     sprintf (string_params[index], "internal_asymmetry[%d]", i);
//                     sprintf (string_params_human_readable[index], "internal_asymmetry[%d]", i);
//                     break;
//                 case 1:
//                     if (similarity_rule[sim_index][0] == '\0')
//                     {
//                         //sprintf (similarity_rule[sim_index], "1 * internal_penalty_by_size[%d] + %.4lf * extrapolation_large_loops", i-1, log(1.0*i/(i-1)));
//                     }
//                     break;
//                 case 2:
//                     array[index] = internal_asymmetry[i];
//                     break;
//                 case 3:
//                     fprintf (file, "%.2lf\n", internal_asymmetry[i]/100.0);
//                     break;
//             }
//             index++;
//         }
//     }




/*
int create_building_block_strings ()
// Mirela: Sep 20, 2005
// For most of the building blocks, write which sequence(s) correspond to each building block
// also fill string_params
// DEPRICATED
{
  int index;
  int i, j, k, l, m, n, o, p;
  index = 0;
  for (i=0; i < NUCL; i++)
    for (j=0; j < NUCL; j++)
      for (k=0; k < NUCL; k++)
        for (l=0; l < NUCL; l++)
          {
            if (stack[i][j][k][l] < INF)
              {
                // exclude duplicates
                // stack[i][j][k][l] is the same as stack[l][k][j][i]
                // TODO: to come back to duplicates
                if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                  {
                    sprintf (string_params[index], "stack[%d][%d][%d][%d]", i, j, k, l);
                    bbseq_left [index][0] = int_to_nuc (i);
                    bbseq_left [index][1] = int_to_nuc (k);
                    bbseq_left [index][2] = '\0';
                    strcpy (bbstr_left[index], "(("); bbstr_left[index][2] = '\0';
                    bbseq_right[index][0] = int_to_nuc (l);
                    bbseq_right[index][1] = int_to_nuc (j);
                    bbseq_right[index][2] = '\0';
                    strcpy (bbstr_right[index], "))"); bbstr_right[index][2] = '\0';
                    index++;
                  }
              }
          }

  for (i=0; i < NUCL; i++)
    for (j=0; j < NUCL; j++)
      for (k=0; k < NUCL; k++)
        for (l=0; l < NUCL; l++)
          {
            if (tstackh[i][j][k][l] < INF)
              {
                // no duplicates here
                sprintf (string_params[index], "tstackh[%d][%d][%d][%d]", i, j, k, l);
                bbseq_left [index][0] = int_to_nuc (i);
                bbseq_left [index][1] = int_to_nuc (k);
                bbseq_left [index][2] = '\0';
                strcpy (bbstr_left[index], "(."); bbstr_left[index][2] = '\0';
                bbseq_right[index][0] = int_to_nuc (l);
                bbseq_right[index][1] = int_to_nuc (j);
                bbseq_right[index][2] = '\0';
                strcpy (bbstr_right[index], ".)"); bbstr_right[index][2] = '\0';
                index++;
              }
          }

    // replace all tstacki parameters by 3 parameters, as described in Mathews 1999
    sprintf (string_params[index], "misc.internal_AU_closure");
    bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
    sprintf (string_params[index], "misc.internal_AG_mismatch");
    bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
    sprintf (string_params[index], "misc.internal_UU_mismatch");
    bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;



//   for (i=0; i < NUCL; i++)
//     for (j=0; j < NUCL; j++)
//       for (k=0; k < NUCL; k++)
//         for (l=0; l < NUCL; l++)
//           {
//             if (tstacki[i][j][k][l] < INF)
//               {
//                 // no duplicates here
//                 sprintf (string_params[index], "tstacki[%d][%d][%d][%d]", i, j, k, l);
//                 // here we have to be careful: this is only half of an internal loop
//                 bbseq_left [index][0] = int_to_nuc (i);
//                 bbseq_left [index][1] = int_to_nuc (k);
//                 bbseq_left [index][2] = '\0';
//                 strcpy (bbstr_left[index], "(."); bbstr_left[index][2] = '\0';
//                 bbseq_right[index][0] = int_to_nuc (l);
//                 bbseq_right[index][1] = int_to_nuc (j);
//                 bbseq_right[index][2] = '\0';
//                 strcpy (bbstr_right[index], ".)"); bbstr_right[index][2] = '\0';
//                 index++;
//               }
//           }

    if (!simple_internal_energy)
    {

        // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {
                                if ( (((i==C && j==G) || (i==G && j==C)) &&
                                     ((m==C && n==G) || (m==G && n==C)) &&
                                     !can_pair(k,l)) ||
                                     (watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                {
                                    if (int11[i][j][k][l][m][n] < INF)
                                    {
                                        // exclude duplicates
                                        // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                        if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                        {
                                            sprintf (string_params[index], "int11[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
                                            bbseq_left [index][0] = int_to_nuc (i);
                                            bbseq_left [index][1] = int_to_nuc (k);
                                            bbseq_left [index][2] = int_to_nuc (m);
                                            bbseq_left [index][3] = '\0';
                                            strcpy (bbstr_left[index], "(.("); bbstr_left[index][3] = '\0';
                                            bbseq_right[index][0] = int_to_nuc (n);
                                            bbseq_right[index][1] = int_to_nuc (l);
                                            bbseq_right[index][2] = int_to_nuc (j);
                                            bbseq_right[index][3] = '\0';
                                            strcpy (bbstr_right[index], ").)"); bbstr_right[index][3] = '\0';
                                            index++;
                                        }
                                    }
                                }
                            }
        sprintf (string_params[index], "misc.internal11_basic_mismatch");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        sprintf (string_params[index], "misc.internal11_GG_mismatch");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;



//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//             for (k=0; k < NUCL; k++)
//                 for (l=0; l < NUCL; l++)
//                 for (m=0; m < NUCL; m++)
//                     for (n=0; n < NUCL; n++)
//                     {
//                         if (int11[i][j][k][l][m][n] < INF)
//                         {
//                             // exclude duplicates
//                             // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
//                             if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
//                             {
//                                 sprintf (string_params[index], "int11[%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n);
//                                 bbseq_left [index][0] = int_to_nuc (i);
//                                 bbseq_left [index][1] = int_to_nuc (k);
//                                 bbseq_left [index][2] = int_to_nuc (m);
//                                 bbseq_left [index][3] = '\0';
//                                 strcpy (bbstr_left[index], "(.("); bbstr_left[index][3] = '\0';
//                                 bbseq_right[index][0] = int_to_nuc (n);
//                                 bbseq_right[index][1] = int_to_nuc (l);
//                                 bbseq_right[index][2] = int_to_nuc (j);
//                                 bbseq_right[index][3] = '\0';
//                                 strcpy (bbstr_right[index], ").)"); bbstr_right[index][3] = '\0';
//                                 index++;
//                             }
//                         }
//                     }


        // go with few parameters, as in Mathews et al 1999
        i=C; j=G; m=C; n=G;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                {
                    if (!can_pair(k,l) && !can_pair(k,o))
                    {
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            sprintf (string_params[index], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                            bbseq_left [index][0] = int_to_nuc (i);
                            bbseq_left [index][1] = int_to_nuc (k);
                            bbseq_left [index][2] = int_to_nuc (m);
                            bbseq_left [index][3] = '\0';
                            strcpy (bbstr_left[index], "(.("); bbstr_left[index][3] = '\0';
                            bbseq_right[index][0] = int_to_nuc (n);
                            bbseq_right[index][1] = int_to_nuc (o);
                            bbseq_right[index][2] = int_to_nuc (l);
                            bbseq_right[index][3] = int_to_nuc (j);
                            bbseq_right[index][4] = '\0';
                            strcpy (bbstr_right[index], ")..)"); bbstr_right[index][4] = '\0';
                            index++;
                        }
                    }
                }
        i=G; j=C; m=G; n=C;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                {
                    if (!can_pair(k,l) && !can_pair(k,o))
                    {
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            sprintf (string_params[index], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
                            bbseq_left [index][0] = int_to_nuc (i);
                            bbseq_left [index][1] = int_to_nuc (k);
                            bbseq_left [index][2] = int_to_nuc (m);
                            bbseq_left [index][3] = '\0';
                            strcpy (bbstr_left[index], "(.("); bbstr_left[index][3] = '\0';
                            bbseq_right[index][0] = int_to_nuc (n);
                            bbseq_right[index][1] = int_to_nuc (o);
                            bbseq_right[index][2] = int_to_nuc (l);
                            bbseq_right[index][3] = int_to_nuc (j);
                            bbseq_right[index][4] = '\0';
                            strcpy (bbstr_right[index], ")..)"); bbstr_right[index][4] = '\0';
                            index++;
                        }
                    }
                }
        sprintf (string_params[index], "misc.internal21_match");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        sprintf (string_params[index], "misc.internal21_AU_closure");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;



//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//             for (k=0; k < NUCL; k++)
//                 for (l=0; l < NUCL; l++)
//                 for (m=0; m < NUCL; m++)
//                     for (n=0; n < NUCL; n++)
//                     for(o=0; o < NUCL; o++)
//                         {
//                         if (int21[i][j][k][l][m][n][o] < INF)
//                             {
//                             // no duplicates here
//                             sprintf (string_params[index], "int21[%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o);
//                             bbseq_left [index][0] = int_to_nuc (i);
//                             bbseq_left [index][1] = int_to_nuc (k);
//                             bbseq_left [index][2] = int_to_nuc (m);
//                             bbseq_left [index][3] = '\0';
//                             strcpy (bbstr_left[index], "(.("); bbstr_left[index][3] = '\0';
//                             bbseq_right[index][0] = int_to_nuc (n);
//                             bbseq_right[index][1] = int_to_nuc (o);
//                             bbseq_right[index][2] = int_to_nuc (l);
//                             bbseq_right[index][3] = int_to_nuc (j);
//                             bbseq_right[index][4] = '\0';
//                             strcpy (bbstr_right[index], ")..)"); bbstr_right[index][4] = '\0';
//                             index++;
//                             }
//                         }


        // go with 53 parameters, like in Mathews et al 1999
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        n = i;
                        m = j;
                        p = k;
                        o = l;
                        if (watson_crick(i,j) && !watson_crick(k,l))
                        {
                            // exclude duplicates
                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                            {
                                sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                bbseq_left [index][0] = int_to_nuc (i);
                                bbseq_left [index][1] = int_to_nuc (k);
                                bbseq_left [index][2] = int_to_nuc (o);
                                bbseq_left [index][3] = int_to_nuc (m);
                                bbseq_left [index][4] = '\0';
                                strcpy (bbstr_left[index], "(..("); bbstr_left[index][4] = '\0';
                                bbseq_right[index][0] = int_to_nuc (n);
                                bbseq_right[index][1] = int_to_nuc (p);
                                bbseq_right[index][2] = int_to_nuc (l);
                                bbseq_right[index][3] = int_to_nuc (j);
                                bbseq_right[index][4] = '\0';
                                strcpy (bbstr_right[index], ")..)"); bbstr_right[index][4] = '\0';
                                index++;
                            }
                        }
                    }

//         i=C; j=G; m=C; n=G;
//         for (k=0; k < NUCL; k++)
//             for (l=0; l < NUCL; l++)
//                 for (o=0; o < NUCL; o++)
//                     for (p=0; p < NUCL; p++)
//                     {
//                         if (!watson_crick(k,l) && !watson_crick(o,p))
//                         {
//                             // exclude duplicates
//                             // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
//                             if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
//                                 n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
//                             {
//                                 sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
//                                 bbseq_left [index][0] = int_to_nuc (i);
//                                 bbseq_left [index][1] = int_to_nuc (k);
//                                 bbseq_left [index][2] = int_to_nuc (o);
//                                 bbseq_left [index][3] = int_to_nuc (m);
//                                 bbseq_left [index][4] = '\0';
//                                 strcpy (bbstr_left[index], "(..("); bbstr_left[index][4] = '\0';
//                                 bbseq_right[index][0] = int_to_nuc (n);
//                                 bbseq_right[index][1] = int_to_nuc (p);
//                                 bbseq_right[index][2] = int_to_nuc (l);
//                                 bbseq_right[index][3] = int_to_nuc (j);
//                                 bbseq_right[index][4] = '\0';
//                                 strcpy (bbstr_right[index], ")..)"); bbstr_right[index][4] = '\0';
//                                 index++;
//                             }
//                         }
//                     }

        // then add the 4 deltas
        sprintf (string_params[index], "misc.internal22_delta_same_size");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        sprintf (string_params[index], "misc.internal22_delta_different_size");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        sprintf (string_params[index], "misc.internal22_delta_1stable_1unstable");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        sprintf (string_params[index], "misc.internal22_delta_AC");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        sprintf (string_params[index], "misc.internal22_match");
        bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;


//         for (i=0; i < NUCL; i++)
//             for (j=0; j < NUCL; j++)
//             for (k=0; k < NUCL; k++)
//                 for (l=0; l < NUCL; l++)
//                 for (m=0; m < NUCL; m++)
//                     for (n=0; n < NUCL; n++)
//                     for(o=0; o < NUCL; o++)
//                         for (p=0; p < NUCL; p++)
//                         {
//                             if (int22[i][j][k][l][m][n][o][p] < INF)
//                             {
//                                 // exclude duplicates
//                                 // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
//                                 if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
//                                     n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
//                                 {
//                                     sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
//                                     bbseq_left [index][0] = int_to_nuc (i);
//                                     bbseq_left [index][1] = int_to_nuc (k);
//                                     bbseq_left [index][2] = int_to_nuc (o);
//                                     bbseq_left [index][3] = int_to_nuc (m);
//                                     bbseq_left [index][4] = '\0';
//                                     strcpy (bbstr_left[index], "(..("); bbstr_left[index][4] = '\0';
//                                     bbseq_right[index][0] = int_to_nuc (n);
//                                     bbseq_right[index][1] = int_to_nuc (p);
//                                     bbseq_right[index][2] = int_to_nuc (l);
//                                     bbseq_right[index][3] = int_to_nuc (j);
//                                     bbseq_right[index][4] = '\0';
//                                     strcpy (bbstr_right[index], ")..)"); bbstr_right[index][4] = '\0';
//                                     index++;
//                                 }
//                             }
//                         }

    }    // end if (!simple_internal_energy)
  for (i=0; i < NUCL; i++)
    for (j=0; j < NUCL; j++)
      for (k=0; k < NUCL; k++)
        {
          if (dangle_top[i][j][k] < INF)
            {
              // no duplicates here
              sprintf (string_params[index], "dangle_top[%d][%d][%d]", i, j, k);
              bbseq_left [index][0] = int_to_nuc (j);
              bbseq_left [index][1] = '\0';
              strcpy (bbstr_left[index], "("); bbstr_left[index][1] = '\0';
              bbseq_right[index][0] = int_to_nuc (i);
              bbseq_right[index][1] = int_to_nuc (k);
              bbseq_right[index][2] = '\0';
              strcpy (bbstr_right[index], ")."); bbstr_right[index][2] = '\0';
              index++;
            }
        }
  for (i=0; i < NUCL; i++)
    for (j=0; j < NUCL; j++)
      for (k=0; k < NUCL; k++)
        {
          if (dangle_bot[i][j][k] < INF)
            {
              // no duplicates here
              sprintf (string_params[index], "dangle_bot[%d][%d][%d]", i, j, k);
              bbseq_left [index][0] = int_to_nuc (k);
              bbseq_left [index][1] = int_to_nuc (j);
              bbseq_left [index][2] = '\0';
              strcpy (bbstr_left[index], ".("); bbstr_left[index][2] = '\0';
              bbseq_right[index][0] = int_to_nuc (i);
              bbseq_right[index][1] = '\0';
              strcpy (bbstr_right[index], ")"); bbstr_right[index][1] = '\0';
              index++;
            }
        }
    int start;
    if (!simple_internal_energy)
        start = 4;
    else
        start = 1;
    for (i=start; i <= MAXLOOP_I; i++)
    {

      if (internal_penalty_by_size[i] < INF)
        {
          // no duplicates here
          sprintf (string_params[index], "internal_penalty_by_size[%d]", i);
          bbseq_left [index][0] = '\0';
          bbseq_right[index][0] = '\0';
          bbstr_left [index][0] = '\0';
          bbstr_right[index][0] = '\0';
          index++;
        }
    }
  for (i=1; i <= MAXLOOP_B; i++)
    {
      if (bulge_penalty_by_size[i] < INF)
        {
          // no duplicates here
          sprintf (string_params[index], "bulge_penalty_by_size[%d]", i);
          bbseq_left [index][0] = '\0';
          bbseq_right[index][0] = '\0';
          bbstr_left [index][0] = '\0';
          bbstr_right[index][0] = '\0';
          index++;
        }
    }
  for (i=1; i <= MAXLOOP_H; i++)
    {
      if (hairpin_penalty_by_size[i] < INF)
        {
          // no duplicates here
          sprintf (string_params[index], "hairpin_penalty_by_size[%d]", i);
          bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
        }
    }

 //sprintf (string_params[index], "misc.param_greater30");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.terminal_AU_penalty");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.hairpin_GGG");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.hairpin_c1");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.hairpin_c2");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.hairpin_c3");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 // TODO: TO ADD the following
 //sprintf (string_params[index], "misc.asymmetry_penalty_max_correction");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
 //sprintf (string_params[index], "misc.asymmetry_penalty_array[0]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
 //sprintf (string_params[index], "misc.asymmetry_penalty_array[1]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 // the next ones are never used
 // sprintf (string_params[index++], "misc.asymmetry_penalty_array[2]");
 // sprintf (string_params[index++], "misc.asymmetry_penalty_array[3]");

 // Instead of these, I will just store the asymmetry for 0.5, 1, 1.5, 2, 2.5 and 3.

 //sprintf (string_params[index], "misc.asymmetry_penalty[1]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 //sprintf (string_params[index], "misc.asymmetry_penalty[2]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 //sprintf (string_params[index], "misc.asymmetry_penalty[3]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 //sprintf (string_params[index], "misc.asymmetry_penalty[4]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 //sprintf (string_params[index], "misc.asymmetry_penalty[5]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 //sprintf (string_params[index], "misc.asymmetry_penalty[6]");
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 //sprintf (string_params[index], "misc.gail_rule");    // this is fixed
 //bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.multi_offset");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.multi_helix_penalty");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.multi_free_base_penalty");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 sprintf (string_params[index], "misc.intermolecular_initiation");
 bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;

 for(i=0; i < nb_triloops; i++)
   {
     // TODO: I can do this
     sprintf (string_params[index], "triloop[%d].energy", i);
     bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
   }

 for(i=0; i < nb_tloops; i++)
   {
     sprintf (string_params[index], "tloop[%d].energy", i);
     bbseq_left [index][0] = '\0'; bbseq_right[index][0] = '\0'; bbstr_left [index][0] = '\0'; bbstr_right[index][0] = '\0'; index++;
   }
 return index;
}
*/







/*
        // I had started to use a crazy model for int22, but I abandoned it.
        // I follow a combination between the model suggested in Christiansen_Znosko_2008 and the one suggested by Mathews_Turner_1999

        // add the middle parts of asymmetric int2x2 (int22mid)
        int num_group1 = 0;
        int num_group2 = 0;
        int num_group3 = 0;
        int num_group4 = 0;
        // strings for the int22mid_groupx
        char sgroup1[200][50];
        char sgroup2[200][50];
        char sgroup3[200][50];
        char sgroup4[200][50];


        if (! parsi_int22)
        {
            // first consider all the symmetric 2x2 internal loops
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                {
                    // i and j must pair
                    if (!can_pair(i,j)) continue;
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                        {
                            // don't know if I should also include the case when k and l pair, I think Dave says yes
                            n = i;
                            m = j;
                            p = k;
                            o = l;
                            //if (!watson_crick(k,l))   // I should include all of them
                            {
                                // exclude duplicates
                                // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                    n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                {
                                    switch (job)
                                    {
                                        case 0:
                                            sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                            sprintf (string_params_human_readable[index], "int22[5'-%c%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(l), int_to_nuc(j));
                                            break;
                                        case 1:
                                            if (similarity_rule[sim_index][0] == '\0')
                                            {
                                                int k_rule1, l_rule1, o_rule1, p_rule1;
                                                apply_rule_1 (k, l, k_rule1, l_rule1);
                                                apply_rule_1 (o, p, o_rule1, p_rule1);
                                                sprintf (similarity_rule[sim_index], "1 * int22[5'-%c%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k_rule1), int_to_nuc(o_rule1), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p_rule1), int_to_nuc(l_rule1), int_to_nuc(j));
                                            }
                                            break;
                                        case 2:
                                            array[index] = int22[i][j][k][l][m][n][o][p];
                                            break;
                                        case 3:
                                            fprintf (file, "%.2lf\n", int22[i][j][k][l][m][n][o][p]/100.0);
                                            break;
                                        case 4:
                                            fgets (buffer, sizeof(buffer), file);
                                            sscanf (buffer, "%lf\n", &param);
                                            param *= 100;
                                            int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
                                            // now the duplicate
                                            int22[n][m][p][o][j][i][l][k] = (PARAMTYPE) param;
                                            break;
                                    }
                                    index++;
                                }
                            }
                        }
                }

            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                    for (k=0;  k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                        {
                            // exclude the ones that correspond to symmetric 2x2 internal loops
                            // Actually we also want the symmetric ones!
                            //if (i==l && j==k)   continue;
                            // for now, exclude the ones that form canonical base pairs.
                            //  However, according to Dave Mathews, they should be there for the partition function calculation.
                            // I guess I can include them directly in the int22 features
                            if (can_pair(i,j))  continue;
                            if (can_pair(k,l))  continue;
                            // they are symmetric, so store only half
                            if (i*1000 + j*100 + k*10 + l <= l*1000+ k*100 + j*10 + i)
                            {
                                switch (job)
                                {
                                    case 0:
                                        sprintf (string_params[index], "int22mid[%d][%d][%d][%d]", i, j, k, l);
                                        sprintf (string_params_human_readable[index], "int22mid[5'-%c%c/%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j));
                                        break;
                                    case 1:
                                        if (similarity_rule[sim_index][0] == '\0')
                                        {
                                            // check which asymmetric 2x2 internal loops have this int22mid, and average
                                            int i_bp, j_bp, ip_bp, jp_bp;
                                            int num_covered_asymmetric = 0;
                                            char s[100];
                                            int index_s;
                                            int there_is_exp = 0;   // is 1 if there is at least one experiment involving this mid
                                            // first traverse to get the counts for correct averaging of the remaining int22mid
                                            // ALSO: count # occurences in each group.
                                            for (i_bp=0; i_bp < NUCL; i_bp++)
                                                for (j_bp=0; j_bp < NUCL; j_bp++)
                                                {
                                                    if (!can_pair(i_bp,j_bp)) continue;
                                                    for (ip_bp=0;  ip_bp < NUCL; ip_bp++)
                                                        for (jp_bp=0; jp_bp < NUCL; jp_bp++)
                                                        {
                                                            if (!can_pair(ip_bp,jp_bp)) continue;
                                                            // exclude the int22 sequence symmetric ones
                                                            if (i_bp==jp_bp && j_bp==ip_bp && i==l && j==k) continue;
                                                            // these are symmetric, so I'm only looking at half of them
                                                            if (i_bp*10000000 + j_bp*1000000 + i*100000 + j*10000 + ip_bp*1000 + jp_bp*100 + k*10 + l <=
                                                                jp_bp*10000000 + ip_bp*1000000 + l*100000 + k*10000 + j_bp*1000 + i_bp*100 + j*10 + i)
                                                            {
                                                                sprintf (s, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i_bp, j_bp, i, j, ip_bp, jp_bp, k, l);
                                                            }
                                                            else        // switch for the other half
                                                            {
                                                                sprintf (s, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", jp_bp, ip_bp, l, k, j_bp, i_bp, j, i);
                                                            }
                                                            //printf ("%s\t%s\n", string_params[index], s);
                                                            int index_s = structure_type_index (s);
                                                            if (similarity_rule[index_s][0] != '\0')
                                                            {
                                                                num_covered_asymmetric++;
                                                                there_is_exp = 1;
                                                            }

                                                        }
                                                }
                                            if (there_is_exp)
                                            {
                                                if (is_int22_group_1 (i,j,k,l))
                                                {
                                                    strcpy (sgroup1[num_group1], string_params_human_readable[index]);
                                                    num_group1++;
                                                }
                                                else if (is_int22_group_2 (i,j,k,l))
                                                {
                                                    strcpy (sgroup2[num_group2], string_params_human_readable[index]);
                                                    num_group2++;
                                                }
                                                else if (is_int22_group_3 (i,j,k,l))
                                                {
                                                    strcpy (sgroup3[num_group3], string_params_human_readable[index]);
                                                    num_group3++;
                                                }
                                                else if (is_int22_group_4 (i,j,k,l))
                                                {
                                                    strcpy (sgroup4[num_group4], string_params_human_readable[index]);
                                                    num_group4++;
                                                }
                                            }
                                            // next, traverse again to fill up similarity rule
                                            // IT MUST BE THE SAME TRAVERSAL AS ABOVE
                                            for (i_bp=0; i_bp < NUCL; i_bp++)
                                                for (j_bp=0; j_bp < NUCL; j_bp++)
                                                {
                                                    if (!can_pair(i_bp,j_bp)) continue;
                                                    for (ip_bp=0;  ip_bp < NUCL; ip_bp++)
                                                        for (jp_bp=0; jp_bp < NUCL; jp_bp++)
                                                        {
                                                            if (!can_pair(ip_bp,jp_bp)) continue;
                                                            // exclude the int22 sequence symmetric ones
                                                            if (i_bp==jp_bp && j_bp==ip_bp && i==l && j==k) continue;
                                                            // these are symmetric, so I'm only looking at half of them
                                                            if (i_bp*10000000 + j_bp*1000000 + i*100000 + j*10000 + ip_bp*1000 + jp_bp*100 + k*10 + l <=
                                                                jp_bp*10000000 + ip_bp*1000000 + l*100000 + k*10000 + j_bp*1000 + i_bp*100 + j*10 + i)
                                                            {
                                                                sprintf (s, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i_bp, j_bp, i, j, ip_bp, jp_bp, k, l);
                                                            }
                                                            else        // switch for the other half
                                                            {
                                                                sprintf (s, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", jp_bp, ip_bp, l, k, j_bp, i_bp, j, i);
                                                            }
                                                            //printf ("%s\t%s\n", string_params[index], s);
                                                            int index_s = structure_type_index (s);
                                                            if (similarity_rule[index_s][0] != '\0')
                                                            {
                                                                // get the indexes of the first and second symmetric int22
                                                                int index_sym1, index_sym2;
                                                                char s_sym1[100];
                                                                char s_sym2[100];
                                                                sprintf (s_sym1, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i_bp, j_bp, i, j, j_bp, i_bp, j, i);
                                                                sprintf (s_sym2, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", jp_bp, ip_bp, l, k, ip_bp, jp_bp, k, l);
                                                                index_sym1 = structure_type_index (s_sym1);
                                                                index_sym2 = structure_type_index (s_sym2);
                                                                // just add plus and a space if sth was added before
                                                                if (similarity_rule[sim_index][0] != '\0')
                                                                    sprintf (similarity_rule[sim_index], "%s + ",
                                                                        similarity_rule[sim_index]);
                                                                sprintf (similarity_rule[sim_index], "%s%g * %s + -%g * %s + -%g * %s",
                                                                    similarity_rule[sim_index],
                                                                    1.0/num_covered_asymmetric,
                                                                    string_params_human_readable[index_s],
                                                                    1.0/num_covered_asymmetric/2.0,
                                                                    string_params_human_readable[index_sym1],
                                                                    1.0/num_covered_asymmetric/2.0,
                                                                    string_params_human_readable[index_sym2]);
                                                            }
                                                        }
                                                }
                                        }
                                        if (similarity_rule[sim_index][0] == '\0')  // it wasn't filled above
                                        {
                                            int i_rule1, j_rule1, k_rule1, l_rule1;
                                            // use the corresponding group, according to Christiansen_Znosko_2008
                                            if (is_int22_group_1 (i,j,k,l))
                                                sprintf (similarity_rule[sim_index], "1 * int22mid_group1");
                                            else if (is_int22_group_2 (i,j,k,l))
                                                sprintf (similarity_rule[sim_index], "1 * int22mid_group2");
                                            else if (is_int22_group_3 (i,j,k,l))
                                                sprintf (similarity_rule[sim_index], "1 * int22mid_group3");
                                            else if (is_int22_group_4 (i,j,k,l))
                                                sprintf (similarity_rule[sim_index], "1 * int22mid_group4");
                                            // follow RULE 1
                                            else if (apply_rule_1 (i,j, i_rule1, j_rule1) || apply_rule_1 (k,l, k_rule1, l_rule1))
                                            {
                                                // int22mid is symmetric, so make sure I use the good one
                                                if (i_rule1*1000 + j_rule1*100 + k_rule1*10 + l_rule1 <= l_rule1*1000+ k_rule1*100 + j_rule1*10 + i_rule1)
                                                    sprintf (similarity_rule[sim_index], "1 * int22mid[5'-%c%c/%c%c-3']",
                                                        int_to_nuc(i_rule1), int_to_nuc(k_rule1), int_to_nuc(l_rule1), int_to_nuc(j_rule1));
                                                else
                                                    sprintf (similarity_rule[sim_index], "1 * int22mid[5'-%c%c/%c%c-3']",
                                                        int_to_nuc(l_rule1), int_to_nuc(j_rule1), int_to_nuc(i_rule1), int_to_nuc(k_rule1));
                                            }
                                        }
                                        break;
                                    case 2:
                                        array[index] = int22mid[i][j][k][l];
                                        break;
                                    case 3:
                                        fprintf (file, "%.2lf\n", int22mid[i][j][k][l]/100.0);
                                        break;
                                    case 4:
                                        fgets (buffer, sizeof(buffer), file);
                                        sscanf (buffer, "%lf\n", &param);
                                        param *= 100;
                                        int22mid[i][j][k][l] = (PARAMTYPE) param;
                                        // now the duplicate
                                        int22mid[l][k][j][i] = (PARAMTYPE) param;
                                        break;
                                }
                                index++;
                            }
                        }
        }   // end if (! parsi_int22)
        // add the four groups

        // add these no matter what
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal22mid_group1");
                sprintf (string_params_human_readable[index], "int22mid_group1");
                break;
            case 1:
                if (num_group1 > 0)
                    sprintf (similarity_rule[sim_index], "%g * %s", 1.0/num_group1, sgroup1[0]);
                for (i=1; i < num_group1; i++)
                    sprintf (similarity_rule[sim_index], "%s + %g * %s", similarity_rule[sim_index], 1.0/num_group1, sgroup1[i]);
                break;
            case 2:
                array[index] = misc.internal22mid_group1;
                break;
            case 3:
                fprintf (file, "%.2lf\n", misc.internal22mid_group1/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal22mid_group1 = (PARAMTYPE) param;
                break;
        }
        index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal22mid_group2");
                sprintf (string_params_human_readable[index], "int22mid_group2");
                break;
            case 1:
                if (num_group2 > 0)
                    sprintf (similarity_rule[sim_index], "%g * %s", 1.0/num_group2, sgroup2[0]);
                for (i=1; i < num_group2; i++)
                    sprintf (similarity_rule[sim_index], "%s + %g * %s", similarity_rule[sim_index], 1.0/num_group2, sgroup2[i]);
                break;
            case 2:
                array[index] = misc.internal22mid_group2;
                break;
            case 3:
                fprintf (file, "%.2lf\n", misc.internal22mid_group2/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal22mid_group2 = (PARAMTYPE) param;
                break;
        }
        index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal22mid_group3");
                sprintf (string_params_human_readable[index], "int22mid_group3");
                break;
            case 1:
                if (num_group3 > 0)
                    sprintf (similarity_rule[sim_index], "%g * %s", 1.0/num_group3, sgroup3[0]);
                for (i=1; i < num_group3; i++)
                    sprintf (similarity_rule[sim_index], "%s + %g * %s", similarity_rule[sim_index], 1.0/num_group3, sgroup3[i]);
                break;
            case 2:
                array[index] = misc.internal22mid_group3;
                break;
            case 3:
                fprintf (file, "%.2lf\n", misc.internal22mid_group3/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal22mid_group3 = (PARAMTYPE) param;
                break;
        }
        index++;
        switch (job)
        {
            case 0:
                sprintf (string_params[index], "misc.internal22mid_group4");
                sprintf (string_params_human_readable[index], "int22mid_group4");
                break;
            case 1:
                if (num_group4 > 0)
                    sprintf (similarity_rule[sim_index], "%g * %s", 1.0/num_group4, sgroup4[0]);
                for (i=1; i < num_group4; i++)
                    sprintf (similarity_rule[sim_index], "%s + %g * %s", similarity_rule[sim_index], 1.0/num_group4, sgroup4[i]);
                break;
            case 2:
                array[index] = misc.internal22mid_group4;
                break;
            case 3:
                fprintf (file, "%.2lf\n", misc.internal22mid_group4/100.0);
                break;
            case 4:
                fgets (buffer, sizeof(buffer), file);
                sscanf (buffer, "%lf\n", &param);
                param *= 100;
                misc.internal22mid_group4 = (PARAMTYPE) param;
                break;
        }
        index++;

        if (parsi_int22)
        {
            // also add misc.internal22_AU_closure
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_AU_closure");
                    sprintf (string_params_human_readable[index], "internal22_AU_closure");
                    break;
                case 2:
                    array[index] = misc.internal22_AU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", misc.internal22_AU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_AU_closure = (PARAMTYPE) param;
                    break;
            }
            index++;
            // also add misc.internal22_GU_closure
            switch (job)
            {
                case 0:
                    sprintf (string_params[index], "misc.internal22_GU_closure");
                    sprintf (string_params_human_readable[index], "internal22_GU_closure");
                    break;
                case 2:
                    array[index] = misc.internal22_GU_closure;
                    break;
                case 3:
                    fprintf (file, "%.2lf\n", misc.internal22_GU_closure/100.0);
                    break;
                case 4:
                    fgets (buffer, sizeof(buffer), file);
                    sscanf (buffer, "%lf\n", &param);
                    param *= 100;
                    misc.internal22_GU_closure = (PARAMTYPE) param;
                    break;
            }
            index++;
        }

        if (!parsi_int22)
        {
            // add the asymmetric int22
            for (i=0; i < NUCL; i++)
                for (j=0; j < NUCL; j++)
                {
                    if (!can_pair(i,j)) continue;
                    for (k=0; k < NUCL; k++)
                        for (l=0; l < NUCL; l++)
                        {
                            // for now, let's only include ncbp in the internal loop
                            //if (watson_crick(k,l)) continue;    // we need to include all
                            for (m=0; m < NUCL; m++)
                                for (n=0; n < NUCL; n++)
                                {
                                    if (!can_pair(m,n)) continue;
                                    for(o=0; o < NUCL; o++)
                                        for (p=0; p < NUCL; p++)
                                        {
                                            //if (watson_crick(o,p)) continue;    // we need to include all
                                            if (!(n==i && m==j && p==k && o==l))       // i.e. asymmetric
                                            {
                                                // exclude duplicates
                                                // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                                if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                                    n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                                {
                                                    switch (job)
                                                    {
                                                        case 0:
                                                            sprintf (string_params[index], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                                            sprintf (string_params_human_readable[index], "int22[5'-%c%c%c%c/%c%c%c%c-3']", int_to_nuc(i), int_to_nuc(k), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(l), int_to_nuc(j));
                                                            break;
                                                        case 1:
                                                            if (similarity_rule[sim_index][0] == '\0')
                                                            {
                                                                int k_rule1, l_rule1, o_rule1, p_rule1;
                                                                apply_rule_1 (k, l, k_rule1, l_rule1);
                                                                apply_rule_1 (o, p, o_rule1, p_rule1);
                                                                char smid[10];
                                                                if (k_rule1*1000 + l_rule1*100 + o_rule1*10 + p_rule1 <= p_rule1*1000 + o_rule1*100 + l_rule1*10 + k_rule1)
                                                                    sprintf (smid, "%c%c/%c%c", int_to_nuc(k_rule1), int_to_nuc(o_rule1), int_to_nuc(p_rule1), int_to_nuc(l_rule1));
                                                                else
                                                                    sprintf (smid, "%c%c/%c%c", int_to_nuc(p_rule1), int_to_nuc(l_rule1), int_to_nuc(k_rule1), int_to_nuc(o_rule1));
                                                                sprintf (similarity_rule[sim_index], "0.5 * int22[5'-%c%c%c%c/%c%c%c%c-3'] + 0.5 * int22[5'-%c%c%c%c/%c%c%c%c-3'] + 1 * int22mid[5'-%s-3']",
                                                                    int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j), int_to_nuc(i), int_to_nuc(k), int_to_nuc(l), int_to_nuc(j),
                                                                    int_to_nuc(n), int_to_nuc(p), int_to_nuc(o), int_to_nuc(m), int_to_nuc(n), int_to_nuc(p), int_to_nuc(o), int_to_nuc(m),
                                                                    smid);
                                                            }
                                                            break;
                                                        case 2:
                                                            array[index] = int22[i][j][k][l][m][n][o][p];
                                                            break;
                                                        case 3:
                                                            fprintf (file, "%.2lf\n", int22[i][j][k][l][m][n][o][p]/100.0);
                                                            break;
                                                        case 4:
                                                            fgets (buffer, sizeof(buffer), file);
                                                            sscanf (buffer, "%lf\n", &param);
                                                            param *= 100;
                                                            int22[i][j][k][l][m][n][o][p] = (PARAMTYPE) param;
                                                            // now the duplicate
                                                            int22[n][m][p][o][j][i][l][k] = (PARAMTYPE) param;
                                                            break;
                                                    }
                                                    index++;
                                                }
                                            }
                                        }
                                }
                        }
                }
        }   // end if (!parsi_int22)
*/



