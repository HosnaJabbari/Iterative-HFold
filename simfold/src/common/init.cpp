/***************************************************************************
                          init.cpp  -  description
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


// this file contains functions to read the thermodynamic parameters from files
//    into internal data structures

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "constants.h"
#include "structs.h"
#include "globals.h"
#include "common.h"
#include "params.h"





int ascii_to_int (char *string)
// PRE:  string is either in float format or it is a '.'
// POST: convert in infinity (INF) if it is '.', in a float otherwise
{
    char *ptr;
    double en;
    if (strcmp(string,".") == 0)
        return INF;
    en = strtod(string, &ptr);
    if (en < 0) en -= EPSILON;
    else en += EPSILON;
    return (int) (en*100);
}

double ascii_to_double (char *string)
// PRE:  string is either in float format or it is a '.'
// POST: convert in infinity (INF) if it is '.', in a float otherwise
{
    char *ptr;
    double en;
    if (strcmp(string,".") == 0)
        return INF;
    en = strtod(string, &ptr);
    return en;
    // took out the *100.0 on Nov 20, 2007, because we don't want this for misc.param_greater30
}

PARAMTYPE ascii_to_param_type (char *string)
{

    #ifdef DOUBLEPARAMS
    return 100.0*ascii_to_double (string);
    // added the *100.0 on Nov 20, 2007. We don't want it inside ascii_to_double because we don't want this for misc.param_greater30
    #else
    return ascii_to_int (string);
    #endif
}


void read_configuration_file (char *filename)
// PRE:  None
// POST: read the configuration file, which must be in the standard format
//      - see the documentation
{
    char buffer[256], token1[50], token2[5], token3[50];
    FILE *file;
    int i;
    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '[' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s", token1, token2, token3);
        for (i=0; i<nb_params; i++)
        {
            if (strcmp(token1, par_name[i]) == 0)
            {
                strcpy (par_value[i], token3);
                break;
            }
        }
        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);
}


void read_loop_file (char *filename, PARAMTYPE internal_penalty_by_size[],
        PARAMTYPE bulge_penalty_by_size[], PARAMTYPE hairpin_penalty_by_size[])
// PRE:  None
// POST: Read information from loop.dat file - internal, bulge, hairpin
{
    char buffer[256];
    FILE *file;
    int size;
    char internal[10], bulge[10], hairpin[10];
    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%d%s%s%s", &size, internal, bulge, hairpin);

        internal_penalty_by_size[size] = ascii_to_param_type(internal);
        bulge_penalty_by_size[size] = ascii_to_param_type(bulge);
        hairpin_penalty_by_size[size] = ascii_to_param_type(hairpin);

        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);
}

void read_tloop_file (char *filename, hairpin_tloop tloop[], int &nb_loops)
// PRE:  None
// POST: Read information from tetra-loop file - hairpin
{
    char buffer[256], energy[10];
    FILE *file;
    int i,j;
    i = 0;
    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s", tloop[i].seq, energy);
        for (j=0; j < strlen(tloop[i].seq); j++)
            tloop[i].seq[j] = toupper(tloop[i].seq[j]);
        tloop[i].energy = ascii_to_param_type(energy);
        fgets (buffer, sizeof(buffer), file);
        i++;    // may be improved - hash-table??
    }
    nb_loops = i;
    fclose (file);
}



void read_int11_file (char *filename, PARAMTYPE int11[][NUCL][NUCL][NUCL][NUCL][NUCL])
// PRE:  filename is the file that contains the energy for 1-1 internal loops
// POST: the values are stored in the 6-D array int11
{
    char buffer[256];
    FILE *file;
    char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
    char  v9[10], v10[10], v11[10], v12[10], v13[10], v14[10], v15[10], v16[10];
    char v17[10], v18[10], v19[10], v20[10], v21[10], v22[10], v23[10], v24[10];
    int ii, jj;                            // ii = i+1; jj = j-1
    int oi[6] = {0,1,2,3,2,3}, oii, oiip;    // the order in which i and ip appear
    int oj[6] = {3,2,1,0,3,2}, ojj, ojjp;    // the order in which j and jp appear
    // In the following configuration:
    //        X
    //      A  C
    //      T  G
    //        Y
    // i = A; j = T; ii = X; jj = Y; ip = C; jp = G
    // int11[i][j][ii][jj][ip][jp]

    // initialize int11

    ii=0; jj=0;                            // ii = i+1; jj = j-1
    oii=0; oiip=0;    // the order in which i and ip appear
    ojj=0; ojjp=0;    // the order in which j and jp appear

    for (int i=0; i < NUCL; i++)
        for (int j=0; j < NUCL; j++)
            for (int ii=0; ii < NUCL; ii++)
                for (int jj=0; jj < NUCL; jj++)
                    for (int ip=0; ip < NUCL; ip++)
                        for (int jp=0; jp < NUCL; jp++)
                        {
                            int11[i][j][ii][jj][ip][jp] = INF;
                            int11_experimental_addition[i][j][ii][jj][ip][jp] = INF;
                        }


    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,
                v18,v19,v20,v21,v22,v23,v24);
        if (ii == NUCL)
        {
            ii = 0; oii++; ojj++;
        }
        // v1-v4 - jj changes
        // ip and jp are always in the following order:
        // AT, CG, GC, TA, GT, TG  - that means
        // 03, 12, 21, 30, 23, 32
        jj = 0; oiip = 0; ojjp = 0;

        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v1);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v2);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v3);
        int11[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v4);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v5);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v6);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v7);
        int11[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v8);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v9);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v10);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v11);
        int11[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v12);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v13);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v14);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v15);
        int11[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v16);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v17);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v18);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v19);
        int11[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v20);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v21);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v22);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v23);
        int11[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]] = ascii_to_param_type(v24);
        // go to the next line: grow ii
        ii++;

        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);
}


void read_int21_file (char *filename, PARAMTYPE int21[][NUCL][NUCL][NUCL][NUCL][NUCL][NUCL])
// PRE:  filename is the file that contains the energy for 2-1 internal loops
// POST: the values are stored in the 7-D array int21
{
    char buffer[256];
    FILE *file;
    char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
    char  v9[10], v10[10], v11[10], v12[10], v13[10], v14[10], v15[10], v16[10];
    char v17[10], v18[10], v19[10], v20[10], v21[10], v22[10], v23[10], v24[10];
    int ii, jj, jjp;                    // ii = i+1; jj = j-1
    int oi[6] = {0,1,2,3,2,3}, oii, oiip;    // the order in which i and ip appear
    int oj[6] = {3,2,1,0,3,2}, ojj, ojjp;    // the order in which j and jp appear
    // In the following configuration:
    //        X
    //      A  C
    //      T  G
    //       YZ
    // i = A; j = T; ii = X; jj = Y; ip = C; jp = G; jjp = Z
    // int21[i][j][ii][jj][ip][jp][jjp]

    ii=0; jj=0; jjp = 0;                    // ii = i+1; jj = j-1
    oii=0; oiip=0;    // the order in which i and ip appear
    ojj=0; ojjp=0;

    // initialize int21

    for (int i=0; i < NUCL; i++)
        for (int j=0; j < NUCL; j++)
            for (int ii=0; ii < NUCL; ii++)
                for (int jj=0; jj < NUCL; jj++)
                    for (int ip=0; ip < NUCL; ip++)
                        for (int jp=0; jp < NUCL; jp++)
                            for (int jjp=0; jjp < NUCL; jjp++)
                            {
                                int21[i][j][ii][jj][ip][jp][jjp] = INF;
                                int21_experimental_addition[i][j][ii][jj][ip][jp][jjp] = INF;
                            }


    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,
                v18,v19,v20,v21,v22,v23,v24);
        if (ii == NUCL)
        {
            ii = 0; jjp++;
        }
        if (jjp == NUCL)
        {
            jjp = 0; oii++; ojj++;
        }
        // v1-v4 - jj changes
        // ip and jp are always in the following order:
        // AT, CG, GC, TA, GT, TG  - that means
        // 03, 12, 21, 30, 23, 32
        jj = 0; oiip = 0; ojjp = 0;
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v1);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v2);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v3);
        int21[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v4);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v5);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v6);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v7);
        int21[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v8);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v9);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v10);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v11);
        int21[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v12);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v13);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v14);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v15);
        int21[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v16);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v17);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v18);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v19);
        int21[oi[oii]][oj[ojj]][ii][jj  ][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v20);
        // j grows, jj starts again from 0
        jj = 0; oiip++; ojjp++;
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v21);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v22);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v23);
        int21[oi[oii]][oj[ojj]][ii][jj++][oi[oiip]][oj[ojjp]][jjp] = ascii_to_param_type(v24);
        // go to the next line: grow ii
        ii++;

        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);
}


void read_int22_file (char *filename, PARAMTYPE int22[][NUCL][NUCL][NUCL][NUCL][NUCL][NUCL][NUCL])
// PRE:  filename is the file that contains the energy for 1-1 internal loops
// POST: the values are stored in the 6-D array int11
{
    char buffer[256];
    FILE *file;
    char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
    char  v9[10], v10[10], v11[10], v12[10], v13[10], v14[10], v15[10], v16[10];
    int ii, jj, iip, jjp;                // ii = i+1; jj = j-1, iip=i'-1, jjp=j'+1
    int oi[6] = {0,1,2,3,2,3};
    int i, ip;    // the order in which i and ip appear
    int oj[6] = {3,2,1,0,3,2};;
    int j, jp;    // the order in which j and jp appear
    // In the following configuration:
    //        XV
    //      A   C
    //      T   G
    //        YZ
    // i = A; j = T; ii = X; jj = Y; ip = C; jp = G; iip = V; jjp = Z
    // int22[i][j][ii][jj][ip][jp][iip][jjp]

    ii=0; jj=0; iip=0; jjp=0;                // ii = i+1; jj = j-1, iip=i'-1, jjp=j'+1


    // initialize int22

    for (i=0; i < NUCL; i++)
        for (int j=0; j < NUCL; j++)
            for (int ii=0; ii < NUCL; ii++)
                for (int jj=0; jj < NUCL; jj++)
                    for (int ip=0; ip < NUCL; ip++)
                        for (int jp=0; jp < NUCL; jp++)
                            for (int iip=0; iip < NUCL; iip++)
                                for (int jjp=0; jjp < NUCL; jjp++)
                                {
                                    int22[i][j][ii][jj][ip][jp][iip][jjp] = INF;
                                    int22_experimental_addition[i][j][ii][jj][ip][jp][iip][jjp] = INF;
                                }

    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    i = 0; ip = 0; j = 0; jp = 0;
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16);
        if (jj == NUCL)
        {
            jj = 0; ii++;
        }
        if (ii == NUCL)
        {
            ii = 0; ip++; jp++;
        }
        if (ip == 6)
        {
            ip = 0; jp = 0; i++; j++;
        }
        // v1-v4 - jj changes
        // ip and jp are always in the following order:
        // AT, CG, GC, TA, GT, TG  - that means
        // 03, 12, 21, 30, 23, 32
        iip = 0; jjp = 0;
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v1);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v2);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v3);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp  ] = ascii_to_param_type(v4);
        // j grows, jj starts again from 0
        iip++; jjp = 0;
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v5);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v6);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v7);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp  ] = ascii_to_param_type(v8);
        // j grows, jj starts again from 0
        iip++; jjp = 0;
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v9);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v10);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v11);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp  ] = ascii_to_param_type(v12);
        // j grows, jj starts again from 0
        iip++; jjp = 0;
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v13);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v14);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v15);
        int22[oi[i]][oj[j]][ii][jj][oi[ip]][oj[jp]][iip][jjp++] = ascii_to_param_type(v16);
        // go to the next line: grow jj
        jj++;

        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);

}


void read_stack_file (char *filename, PARAMTYPE stack[][NUCL][NUCL][NUCL])
// PRE:  filename is the file that contains the stacked pairs energy data for DNA
// POST: the values are stored in the 4-D array stack
{
    char buffer[256];
    FILE *file;
    char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
    char  v9[10], v10[10], v11[10], v12[10], v13[10], v14[10], v15[10], v16[10];
    int i, j, ii, jj;       // ii = i+1; jj = j-1

    i=0; j=0; ii=0; jj=0;
    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16);
        if (ii == NUCL)
        {
            ii = 0; i++;
        }
        // v1-v4 - jj changes
        jj = 0; j = 0;

        stack[i][j][ii][jj++] = ascii_to_param_type(v1);
        stack[i][j][ii][jj++] = ascii_to_param_type(v2);
        stack[i][j][ii][jj++] = ascii_to_param_type(v3);
        stack[i][j][ii][jj  ] = ascii_to_param_type(v4);
        // j grows, jj starts again from 0
        jj = 0; j++;
        stack[i][j][ii][jj++] = ascii_to_param_type(v5);
        stack[i][j][ii][jj++] = ascii_to_param_type(v6);
        stack[i][j][ii][jj++] = ascii_to_param_type(v7);
        stack[i][j][ii][jj  ] = ascii_to_param_type(v8);
        // j grows, jj starts again from 0
        jj = 0; j++;
        stack[i][j][ii][jj++] = ascii_to_param_type(v9);
        stack[i][j][ii][jj++] = ascii_to_param_type(v10);
        stack[i][j][ii][jj++] = ascii_to_param_type(v11);
        stack[i][j][ii][jj  ] = ascii_to_param_type(v12);
        // j grows, jj starts again from 0
        jj = 0; j++;
        stack[i][j][ii][jj++] = ascii_to_param_type(v13);
        stack[i][j][ii][jj++] = ascii_to_param_type(v14);
        stack[i][j][ii][jj++] = ascii_to_param_type(v15);
        stack[i][j][ii][jj++] = ascii_to_param_type(v16);
        // go to the next line: grow ii
        ii++;

        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);
}


void read_dangle_file (char *filename, PARAMTYPE dangle_top[][NUCL][NUCL], PARAMTYPE dangle_bot[][NUCL][NUCL])
// PRE:  filename is the file that contains the stacked pairs energy data for DNA
// POST: the values are stored in the 4-D array stack
{
    char buffer[256];
    FILE *file;
    char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
    char  v9[10], v10[10], v11[10], v12[10], v13[10], v14[10], v15[10], v16[10];
    int i, j, ii, jj;       // ii = i+1; jj = j-1

    i=0; j=0; ii=0; jj=0;
    int counter;
    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    i = 0;
    counter = 0;
    // read the top dangling pairs
    while (counter < 4 && !feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16);
        ii = 0; j = 0;

        dangle_top[i][j][ii++] = ascii_to_param_type(v1);
        dangle_top[i][j][ii++] = ascii_to_param_type(v2);
        dangle_top[i][j][ii++] = ascii_to_param_type(v3);
        dangle_top[i][j][ii  ] = ascii_to_param_type(v4);
        // j grows, ii starts again from 0
        ii = 0; j++;
        dangle_top[i][j][ii++] = ascii_to_param_type(v5);
        dangle_top[i][j][ii++] = ascii_to_param_type(v6);
        dangle_top[i][j][ii++] = ascii_to_param_type(v7);
        dangle_top[i][j][ii  ] = ascii_to_param_type(v8);
        // j grows, ii starts again from 0
        ii = 0; j++;
        dangle_top[i][j][ii++] = ascii_to_param_type(v9);
        dangle_top[i][j][ii++] = ascii_to_param_type(v10);
        dangle_top[i][j][ii++] = ascii_to_param_type(v11);
        dangle_top[i][j][ii  ] = ascii_to_param_type(v12);
        // j grows, ii starts again from 0
        ii = 0; j++;
        dangle_top[i][j][ii++] = ascii_to_param_type(v13);
        dangle_top[i][j][ii++] = ascii_to_param_type(v14);
        dangle_top[i][j][ii++] = ascii_to_param_type(v15);
        dangle_top[i][j][ii  ] = ascii_to_param_type(v16);

        fgets (buffer, sizeof(buffer), file);
        i++;
        counter++;
    }
    i = 0;
    counter = 0;
    // read the bottom dangling pairs
    while (counter < 4 && !feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
                v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16);
        jj = 0; j = 0;
        dangle_bot[i][j][jj++] = ascii_to_param_type(v1);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v2);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v3);
        dangle_bot[i][j][jj  ] = ascii_to_param_type(v4);
        // j grows, ii starts again from 0
        jj = 0; j++;
        dangle_bot[i][j][jj++] = ascii_to_param_type(v5);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v6);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v7);
        dangle_bot[i][j][jj  ] = ascii_to_param_type(v8);
        // j grows, ii starts again from 0
        jj = 0; j++;
        dangle_bot[i][j][jj++] = ascii_to_param_type(v9);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v10);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v11);
        dangle_bot[i][j][jj  ] = ascii_to_param_type(v12);
        // j grows, ii starts again from 0
        jj = 0; j++;
        dangle_bot[i][j][jj++] = ascii_to_param_type(v13);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v14);
        dangle_bot[i][j][jj++] = ascii_to_param_type(v15);
        dangle_bot[i][j][jj  ] = ascii_to_param_type(v16);

        fgets (buffer, sizeof(buffer), file);
        i++;
        counter++;
    }

    fclose (file);
}


void read_sequence_file (char *filename, char *sequence)
// PRE:  filename is the input file with the pool of input words
// POST: Read data from file and put them in the global variable pool
{
    char buffer[MAXSLEN];
    FILE *file;

    // make sure the sequence is empty
    empty_string (sequence);
    if ((file = fopen (filename, "r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }
    fgets (buffer, sizeof(buffer), file);
    while (!feof(file))
    {
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        if (buffer[strlen(buffer)-1] == '\n')
            buffer[strlen(buffer)-1] = '\0';
        if (strcmp (sequence, "") == 0)
            strcpy (sequence, buffer);
        else
            strcat (sequence, buffer);
        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);
}


void read_miscloop_file (char *filename, miscinfo &misc)
// read the information from the miscloop file
{
    char buffer[256];
    FILE *file;
    char  v1[10],  v2[10],  v3[10],  v4[10];
    if ((file = fopen (filename,"r")) == NULL)
    {
        giveup ("Cannot open file", filename);
    }

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.param_greater30 = ascii_to_double (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.asymmetry_penalty_max_correction = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s%s%s%s", v1, v2, v3, v4);
    misc.asymmetry_penalty_array[0] = ascii_to_param_type (v1);
    misc.asymmetry_penalty_array[1] = ascii_to_param_type (v2);
    misc.asymmetry_penalty_array[2] = ascii_to_param_type (v3);
    misc.asymmetry_penalty_array[3] = ascii_to_param_type (v4);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s%s%s", v1, v2, v3);
    misc.multi_offset = ascii_to_param_type (v1);
    misc.multi_free_base_penalty = ascii_to_param_type (v2);
    misc.multi_helix_penalty = ascii_to_param_type (v3);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    // efn2 multibranched loops - not in use
    sscanf (buffer, "%s%s%s", v1, v2, v3);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.terminal_AU_penalty = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.hairpin_GGG = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.hairpin_c1 = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.hairpin_c2 = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.hairpin_c3 = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.intermolecular_initiation = ascii_to_param_type (v1);

    fgets (buffer, sizeof(buffer), file);
    while (buffer[0] == '#' || buffer[0] == '\n')
        fgets (buffer, sizeof(buffer), file);
    sscanf (buffer, "%s", v1);
    misc.gail_rule = atoi (v1);
	fclose (file);
}

int round_double (double number)
{
    int intrep = (int)number;
    if (number < 0)
    {
        if (number-intrep <= -0.5) return intrep-1;
        else return intrep;
    }
    else
    {
        if (number-intrep >= 0.5) return intrep+1;
        else return intrep;
    }
}

void formula (double temp, PARAMTYPE &energy, PARAMTYPE enthalpy)
{
    if (energy < INF)
    {
        // TODO: this is used only for temperature != 37, so it's ok for now, but maybe it shouldn't be converted to int
        energy = (int)(round_double(((double)(enthalpy) - (temp + 273.15)*((double)(enthalpy-energy))/310.15)));
    }
}


void calculate_energies (double temp)
{
    int i, j, k, l, m, n, o, p;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                formula (temp, dangle_top[i][j][k], enthalpy_dangle_top[i][j][k]);
                formula (temp, dangle_bot[i][j][k], enthalpy_dangle_bot[i][j][k]);
                for (l=0; l < NUCL; l++)
                {
                    formula (temp, stack[i][j][k][l], enthalpy_stack[i][j][k][l]);
                    formula (temp, tstackh[i][j][k][l], enthalpy_tstackh[i][j][k][l]);
                    formula (temp, tstacki[i][j][k][l], enthalpy_tstacki[i][j][k][l]);
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            formula (temp, int11[i][j][k][l][m][n], enthalpy_int11[i][j][k][l][m][n]);
                            for (o=0; o < NUCL; o++)
                            {
                                formula (temp, int21[i][j][k][l][m][n][o], enthalpy_int21[i][j][k][l][m][n][o]);
                                for (p=0; p < NUCL; p++)
                                    formula (temp, int22[i][j][k][l][m][n][o][p], enthalpy_int22[i][j][k][l][m][n][o][p]);
                            }
                        }
                }
            }
    for (i=0; i < MAXLOOP+1; i++)
    {
        formula (temp, internal_penalty_by_size[i], enthalpy_internal_penalty_by_size[i]);
        formula (temp, bulge_penalty_by_size[i], enthalpy_bulge_penalty_by_size[i]);
        formula (temp, hairpin_penalty_by_size[i], enthalpy_hairpin_penalty_by_size[i]);
    }

    formula (temp, misc.terminal_AU_penalty, enthalpy_misc.terminal_AU_penalty);
    formula (temp, misc.hairpin_GGG, enthalpy_misc.hairpin_GGG);
    formula (temp, misc.hairpin_c1, enthalpy_misc.hairpin_c1);
    formula (temp, misc.hairpin_c2, enthalpy_misc.hairpin_c2);
    formula (temp, misc.hairpin_c3, enthalpy_misc.hairpin_c3);
    formula (temp, misc.asymmetry_penalty_max_correction, enthalpy_misc.asymmetry_penalty_max_correction);
    formula (temp, misc.asymmetry_penalty_array[0], enthalpy_misc.asymmetry_penalty_array[0]);
    formula (temp, misc.asymmetry_penalty_array[1], enthalpy_misc.asymmetry_penalty_array[1]);
    formula (temp, misc.asymmetry_penalty_array[2], enthalpy_misc.asymmetry_penalty_array[2]);
    formula (temp, misc.asymmetry_penalty_array[3], enthalpy_misc.asymmetry_penalty_array[3]);
    formula (temp, misc.multi_offset, enthalpy_misc.multi_offset);
    formula (temp, misc.multi_helix_penalty, enthalpy_misc.multi_helix_penalty);
    formula (temp, misc.multi_free_base_penalty, enthalpy_misc.multi_free_base_penalty);
    formula (temp, misc.intermolecular_initiation, enthalpy_misc.intermolecular_initiation);

    // added enthalpy_nb_tloops on August 6th
    for (i=0; i < enthalpy_nb_triloops; i++)
        formula (temp, triloop[i].energy, enthalpy_triloop[i].energy);

    for (i=0; i < enthalpy_nb_tloops; i++)
        formula (temp, tloop[i].energy, enthalpy_tloop[i].energy);
}

//void init_data (char *config_file, int what, double temperature)
void init_data(const char *arg, char *config_file, int what, double temperature)
// the function that must be called by the main program to read data files
// PRE:  None
// POST: Read all data and configuration files
{
    char stack_energy37_filename[200], loop_energy37_filename[200];
    char tloop_energy37_filename[200], tstackh_energy37_filename[200];
    char tstacki_energy37_filename[200], int11_energy37_filename[200];
    char int21_energy37_filename[200], int22_energy37_filename[200];
    char miscloop_energy37_filename[200], dangle_energy37_filename[200];
    char triloop_energy37_filename[200];

    char stack_enthalpy_filename[200], loop_enthalpy_filename[200];
    char tloop_enthalpy_filename[200], tstackh_enthalpy_filename[200];
    char tstacki_enthalpy_filename[200], int11_enthalpy_filename[200];
    char int21_enthalpy_filename[200], int22_enthalpy_filename[200];
    char miscloop_enthalpy_filename[200], dangle_enthalpy_filename[200];
    char triloop_enthalpy_filename[200];
    char conf[200];

    // Added on Mar 19, 2008
    char special_hl_energy37_filename[1000];
    char special_t99_l_energy37_filename[1000];

    // configuration file
    // first find the path to params
    int i, j, k, l, m, n, o, p, ip, jp;
    int len;
    int index;
    char path[200];
    char separator;
    char config_dir[200];

    len = strlen(arg);
    index = -1;
    separator = '/';

    strcpy (path, "");
    for (i = len; i >=0; i--)
    {
        // make it work on Linux
        if (arg[i] == '/')
        {
            separator = '/';
            index = i;
            break;
        }
        // make it work on Windows
        else if (arg[i] == '\\')
        {
            separator = '\\';
            index = i;
            break;
        }
    }
    if (index > -1)
	{
		for (i=0; i < index+1; i++)
			path[i] = arg[i];
		path[i] = '\0';
	}

        strncpy (path, arg, index+1);

    // get the path of the configuration directory
    strcpy (conf, config_file);

    len = strlen(conf);
    index = -1;
    separator = '/';

    strcpy (config_dir, "");
    for (i = len; i >=0; i--)
    {
        // make it work on Linux
        if (conf[i] == '/')
        {
            separator = '/';
            index = i;
            break;
        }
        // make it work on Windows
        else if (conf[i] == '\\')
        {
            separator = '\\';
            index = i;
            break;
        }
    }

    if (index > -1)
	{
		for (i=0; i < index+1; i++)
			config_dir[i] = conf[i];
		config_dir[i] = '\0';
	}
    //printf ("config_dir: %s\n", config_dir);
    //if (separator == '/')
    strcat (path, config_dir);
    //else if (separator == '\\')
    //    strcat (path, "params\\");

    if (what != RNA && what != DNA)
    {
        fprintf (stderr, "Please specify what to fold: RNA or DNA\n");
        exit(1);
    }
    if (temperature < 0 || temperature > 100)
    {
        fprintf (stderr, "Temperature must be between 0 and 100 degrees Celsius\n");
        exit(1);
    }

    read_configuration_file (conf);
    strcpy (std_dir_par, config_dir);

    // if temperature is 37, there's no need to read the enthalpies
    if (what == RNA)
    {
        if (temperature == 37)
        {
            strcpy (stack_energy37_filename, std_dir_par);
            strcat (stack_energy37_filename, rna_stack_energy37_v31_par);
            strcpy (loop_energy37_filename, std_dir_par);
            strcat (loop_energy37_filename, rna_loop_energy37_v31_par);
            strcpy (triloop_energy37_filename, std_dir_par);
            strcat (triloop_energy37_filename, rna_triloop_energy37_v31_par);
            strcpy (tloop_energy37_filename, std_dir_par);
            strcat (tloop_energy37_filename, rna_tloop_energy37_v31_par);
            strcpy (tstackh_energy37_filename, std_dir_par);
            strcat (tstackh_energy37_filename, rna_tstackh_energy37_v31_par);
            strcpy (tstacki_energy37_filename, std_dir_par);
            strcat (tstacki_energy37_filename, rna_tstacki_energy37_v31_par);
            strcpy (int11_energy37_filename, std_dir_par);
            strcat (int11_energy37_filename, rna_int11_energy37_v31_par);
            strcpy (int21_energy37_filename, std_dir_par);
            strcat (int21_energy37_filename, rna_int21_energy37_v31_par);
            strcpy (int22_energy37_filename, std_dir_par);
            strcat (int22_energy37_filename, rna_int22_energy37_v31_par);
            strcpy (miscloop_energy37_filename, std_dir_par);
            strcat (miscloop_energy37_filename, rna_miscloop_energy37_v31_par);
            strcpy (dangle_energy37_filename, std_dir_par);
            strcat (dangle_energy37_filename, rna_dangle_energy37_v31_par);
            strcpy (special_hl_energy37_filename, std_dir_par);
            strcat (special_hl_energy37_filename, rna_special_hl_energy_par);
            strcpy (special_t99_l_energy37_filename, std_dir_par);
            strcat (special_t99_l_energy37_filename, rna_special_t99_l_energy_par);
        }
        else
        {
            strcpy (stack_energy37_filename, std_dir_par);
            strcat (stack_energy37_filename, rna_stack_energy37_v23_par);
            strcpy (loop_energy37_filename, std_dir_par);
            strcat (loop_energy37_filename, rna_loop_energy37_v23_par);
            strcpy (triloop_energy37_filename, std_dir_par);
            strcat (triloop_energy37_filename, rna_triloop_energy37_v23_par);
            strcpy (tloop_energy37_filename, std_dir_par);
            strcat (tloop_energy37_filename, rna_tloop_energy37_v23_par);
            strcpy (tstackh_energy37_filename, std_dir_par);
            strcat (tstackh_energy37_filename, rna_tstackh_energy37_v23_par);
            strcpy (tstacki_energy37_filename, std_dir_par);
            strcat (tstacki_energy37_filename, rna_tstacki_energy37_v23_par);
            strcpy (int11_energy37_filename, std_dir_par);
            strcat (int11_energy37_filename, rna_int11_energy37_v23_par);
            strcpy (int21_energy37_filename, std_dir_par);
            strcat (int21_energy37_filename, rna_int21_energy37_v23_par);
            strcpy (int22_energy37_filename, std_dir_par);
            strcat (int22_energy37_filename, rna_int22_energy37_v23_par);
            strcpy (miscloop_energy37_filename, std_dir_par);
            strcat (miscloop_energy37_filename, rna_miscloop_energy37_v23_par);
            strcpy (dangle_energy37_filename, std_dir_par);
            strcat (dangle_energy37_filename, rna_dangle_energy37_v23_par);
        }

        strcpy (stack_enthalpy_filename, std_dir_par);
        strcat (stack_enthalpy_filename, rna_stack_enthalpy_par);
        strcpy (loop_enthalpy_filename, std_dir_par);
        strcat (loop_enthalpy_filename, rna_loop_enthalpy_par);
        strcpy (triloop_enthalpy_filename, std_dir_par);
        strcat (triloop_enthalpy_filename, rna_triloop_enthalpy_par);
        strcpy (tloop_enthalpy_filename, std_dir_par);
        strcat (tloop_enthalpy_filename, rna_tloop_enthalpy_par);
        strcpy (tstackh_enthalpy_filename, std_dir_par);
        strcat (tstackh_enthalpy_filename, rna_tstackh_enthalpy_par);
        strcpy (tstacki_enthalpy_filename, std_dir_par);
        strcat (tstacki_enthalpy_filename, rna_tstacki_enthalpy_par);
        strcpy (int11_enthalpy_filename, std_dir_par);
        strcat (int11_enthalpy_filename, rna_int11_enthalpy_par);
        strcpy (int21_enthalpy_filename, std_dir_par);
        strcat (int21_enthalpy_filename, rna_int21_enthalpy_par);
        strcpy (int22_enthalpy_filename, std_dir_par);
        strcat (int22_enthalpy_filename, rna_int22_enthalpy_par);
        strcpy (miscloop_enthalpy_filename, std_dir_par);
        strcat (miscloop_enthalpy_filename, rna_miscloop_enthalpy_par);
        strcpy (dangle_enthalpy_filename, std_dir_par);
        strcat (dangle_enthalpy_filename, rna_dangle_enthalpy_par);
    }

    else if (what == DNA)
    {
        strcpy (stack_energy37_filename, std_dir_par);
        strcat (stack_energy37_filename, dna_stack_energy37_par);
        strcpy (loop_energy37_filename, std_dir_par);
        strcat (loop_energy37_filename, dna_loop_energy37_par);
        strcpy (triloop_energy37_filename, std_dir_par);
        strcat (triloop_energy37_filename, dna_triloop_energy37_par);
        strcpy (tloop_energy37_filename, std_dir_par);
        strcat (tloop_energy37_filename, dna_tloop_energy37_par);
        strcpy (tstackh_energy37_filename, std_dir_par);
        strcat (tstackh_energy37_filename, dna_tstackh_energy37_par);
        strcpy (tstacki_energy37_filename, std_dir_par);
        strcat (tstacki_energy37_filename, dna_tstacki_energy37_par);
        strcpy (int11_energy37_filename, std_dir_par);
        strcat (int11_energy37_filename, dna_int11_energy37_par);
        strcpy (int21_energy37_filename, std_dir_par);
        strcat (int21_energy37_filename, dna_int21_energy37_par);
        strcpy (int22_energy37_filename, std_dir_par);
        strcat (int22_energy37_filename, dna_int22_energy37_par);
        strcpy (miscloop_energy37_filename, std_dir_par);
        strcat (miscloop_energy37_filename, dna_miscloop_energy37_par);
        strcpy (dangle_energy37_filename, std_dir_par);
        strcat (dangle_energy37_filename, dna_dangle_energy37_par);

        // write the enthalpy files
        strcpy (stack_enthalpy_filename, std_dir_par);
        strcat (stack_enthalpy_filename, dna_stack_enthalpy_par);
        strcpy (loop_enthalpy_filename, std_dir_par);
        strcat (loop_enthalpy_filename, dna_loop_enthalpy_par);
        strcpy (triloop_enthalpy_filename, std_dir_par);
        strcat (triloop_enthalpy_filename, dna_triloop_enthalpy_par);
        strcpy (tloop_enthalpy_filename, std_dir_par);
        strcat (tloop_enthalpy_filename, dna_tloop_enthalpy_par);
        strcpy (tstackh_enthalpy_filename, std_dir_par);
        strcat (tstackh_enthalpy_filename, dna_tstackh_enthalpy_par);
        strcpy (tstacki_enthalpy_filename, std_dir_par);
        strcat (tstacki_enthalpy_filename, dna_tstacki_enthalpy_par);
        strcpy (int11_enthalpy_filename, std_dir_par);
        strcat (int11_enthalpy_filename, dna_int11_enthalpy_par);
        strcpy (int21_enthalpy_filename, std_dir_par);
        strcat (int21_enthalpy_filename, dna_int21_enthalpy_par);
        strcpy (int22_enthalpy_filename, std_dir_par);
        strcat (int22_enthalpy_filename, dna_int22_enthalpy_par);
        strcpy (miscloop_enthalpy_filename, std_dir_par);
        strcat (miscloop_enthalpy_filename, dna_miscloop_enthalpy_par);
        strcpy (dangle_enthalpy_filename, std_dir_par);
        strcat (dangle_enthalpy_filename, dna_dangle_enthalpy_par);
    }

    read_stack_file (stack_energy37_filename, stack);
    read_loop_file (loop_energy37_filename, internal_penalty_by_size, bulge_penalty_by_size, hairpin_penalty_by_size);
    if (parsi_special == T99)
    {
        read_tloop_file (triloop_energy37_filename, triloop, nb_triloops);
        read_tloop_file (tloop_energy37_filename, tloop, nb_tloops);
    }
    else if (parsi_special == LAVISH)
    {
        read_tloop_file (special_hl_energy37_filename, special_hl, nb_special_hl);
    }
    else if (parsi_special == T99_LAVISH)
    {
        read_tloop_file (special_t99_l_energy37_filename, special_hl, nb_special_hl);
        //printf ("nb_special_hl = %d\n", nb_special_hl);
    }
    read_stack_file (tstackh_energy37_filename, tstackh);
    read_stack_file (tstacki_energy37_filename, tstacki);
    read_miscloop_file (miscloop_energy37_filename, misc);

    // now fill the 3 misc.internal parameters, maybe we need them
    // this is a good place, because we have read tstacki and misc.terminal_AU_penalty
    // Actually better to get the values from Chen_Turner_2006b
//     misc.internal_AU_closure = 65;
//     misc.internal_AG_mismatch = tstacki[C][G][A][G];
//     misc.internal_GA_mismatch = tstacki[C][G][G][A];
//     misc.internal_GG_mismatch = tstacki[C][G][G][G];
//     misc.internal_UU_mismatch = tstacki[C][G][U][U];

    misc.internal_AU_closure = 73;   // does not need terminal_AU_penalty to be added to it
    misc.internal_GA_AG_mismatch = -91;  // this is just for the simple model
    //#if (MODEL == EXTENDED)
    // use GA and AG separately, and add GG, according to Schroeder_Turner_2000
    misc.internal_AG_mismatch = -39;
    misc.internal_GA_mismatch = -120;
    misc.internal_GG_mismatch = -74;
    if (parsi_special)
    {
        misc.internal_special_3GA = 0;
        misc.internal_special_2GA = 0;
        misc.internal_special_2xGA_GC = 0;
        misc.internal_special_midGA = 0;
        misc.internal_special_UG_AG = 0;
        misc.internal_special_GU_A = 0;
    }
    else
    {
        misc.internal_special_3GA = -236;
        misc.internal_special_2GA = -118;
        misc.internal_special_2xGA_GC = -96;
        misc.internal_special_midGA = -91;
        misc.internal_special_UG_AG = -95;
        misc.internal_special_GU_A = 96;
    }
    //#endif
    misc.internal_UU_mismatch = -34;

    read_dangle_file (dangle_energy37_filename, dangle_top, dangle_bot);
    read_int11_file (int11_energy37_filename, int11);
    //#if (MODEL == SIMPLE)
    misc.internal11_basic_mismatch = 40;
    misc.internal11_GG_mismatch = -170;
    //#elif (MODEL == EXTENDED)
    // values from Davis_Znosko_2007
    misc.internal11_AU_closure = 120;
    misc.internal11_GU_closure = 120;
    misc.internal11_AG_mismatch = -40;
    misc.internal11_GG_mismatch = -210;
    misc.internal11_UU_mismatch = -30;
    misc.internal11_5YRR_5YRR = 70;
    misc.internal11_5RYY_5RYY = -50;
    misc.internal11_5YYR_5YYR = 40;
    misc.internal11_5YRY_5RYR = -40;
    misc.internal11_5RRY_5RYY = -100;
    //#endif

    read_int21_file (int21_energy37_filename, int21);
    //#if (MODEL == SIMPLE)
    misc.internal21_match = 400;
    misc.internal21_AU_closure = 75;
    //#elif (MODEL == EXTENDED)
    // values from Badhwar_Znosko_2007
    misc.internal21_initiation = 220;
    misc.internal21_AU_closure = 70;
    misc.internal21_GU_closure = 70;
    misc.internal21_AG_mismatch = -110;
    misc.internal21_GG_mismatch = -110;
    misc.internal21_UU_mismatch = -100;
    //#endif

    read_int22_file (int22_energy37_filename, int22);

    // fill the misc.internal22 parameters
    //misc.internal22_delta_same_size = int22[G][C][A][A][G][C][A][G] - (int22[G][C][A][A][C][G][A][A] + int22[C][G][G][A][G][C][A][G])/2;
    //misc.internal22_delta_different_size = int22[G][C][G][U][G][C][A][G] - (int22[G][C][G][U][C][G][U][G] + int22[C][G][G][A][G][C][A][G])/2;
    // the above don't give the right numbers, for some reason. Just use the numbers from Mathews et al 1999
    misc.internal22_delta_same_size = 0;
    misc.internal22_delta_different_size = 180;
    misc.internal22_delta_1stable_1unstable = 100;
    misc.internal22_delta_AC = 0;
    misc.internal22_match = 200;

    //#if (MODEL == EXTENDED)
    // fill up initial values according to Christiansen_Znosko_2008
    misc.internal22mid_group1 = 110;
    misc.internal22mid_group2 = -120;
    misc.internal22mid_group3 = 80;
    misc.internal22mid_group4 = -30;
    misc.internal22_AU_closure = 50;
    misc.internal22_GU_closure = 120;

    // internal asymmetries are obtained using matlab's nlinfit on a logarithmic function, see MultiRNAFold/FIT_CURVES
    internal_asymmetry_initiation = 18;
    internal_asymmetry_slope = 113;
    if (!parsi_asymmetry)
    {
        for (i=1; i < MAXLOOP; i++)
        {
            internal_asymmetry[i] = (PARAMTYPE) (internal_asymmetry_initiation + internal_asymmetry_slope * log(i));
        }
    }

    // for bulges, the default value is 0
    bulgeA = bulge_penalty_by_size[1];
    bulgeC = bulge_penalty_by_size[1];
    bulgeG = bulge_penalty_by_size[1];
    bulgeU = bulge_penalty_by_size[1];
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



    misc.hairpin_AU_closure = 60;  // according to Serra_Temel_1997
    misc.hairpin_AG_mismatch = (PARAMTYPE) (tstackh[C][G][A][G]/2.0 + tstackh[G][C][A][G]/2.0);
    misc.hairpin_GA_mismatch = (PARAMTYPE) (tstackh[C][G][G][A]/2.0 + tstackh[G][C][G][A]/2.0);
    misc.hairpin_UU_mismatch = (PARAMTYPE) (tstackh[C][G][U][U]/2.0 + tstackh[G][C][U][U]/2.0);

    //#endif

    // read enthalpy files
    read_stack_file (stack_enthalpy_filename, enthalpy_stack);
    read_loop_file (loop_enthalpy_filename, enthalpy_internal_penalty_by_size, enthalpy_bulge_penalty_by_size, enthalpy_hairpin_penalty_by_size);
    read_tloop_file (triloop_enthalpy_filename, enthalpy_triloop, enthalpy_nb_triloops);
    read_tloop_file (tloop_enthalpy_filename, enthalpy_tloop, enthalpy_nb_tloops);
    read_stack_file (tstackh_enthalpy_filename, enthalpy_tstackh);
    read_stack_file (tstacki_enthalpy_filename, enthalpy_tstacki);
    read_miscloop_file (miscloop_enthalpy_filename, enthalpy_misc);
    read_dangle_file (dangle_enthalpy_filename, enthalpy_dangle_top, enthalpy_dangle_bot);
    read_int11_file (int11_enthalpy_filename, enthalpy_int11);
    read_int21_file (int21_enthalpy_filename, enthalpy_int21);
    read_int22_file (int22_enthalpy_filename, enthalpy_int22);


    extrapolate_parameters();

    // OBSOLETE
    // add internal_penalty_by_size_2D
//     for (i=1; i < MAXLOOP_I; i++)
//     {
//         for (j=MAX(i,3); j < MAXLOOP_I; j++)
//         {
//             internal_penalty_by_size_2D[i][j] = internal_penalty_by_size[i+j] + asymmetry_penalty(i,j);
//         }
//     }

    // in the simple model this was considered 0
    // we start from 3 because we have full tables for int1x1 and int2x2
//     for (i=3; i <= MAXLOOP_I/2; i++)
//         internal_symmetry[i] = 0;
    // we initialize with the Ninio values for the simple model
//     for (i=1; i <= 5; i++)
//         internal_asymmetry[i] = i*100.0/2.0;
//     for (i=6; i <= MAXLOOP_I-2; i++)
//         internal_asymmetry[i] = 300;
//     #endif

    if (temperature != 37)
    {
        calculate_energies (temperature);
    }
}
