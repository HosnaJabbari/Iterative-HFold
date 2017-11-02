/***************************************************************************
                          s_specific_functions.h  -  description
                             -------------------
    begin                : Thu Sep 5 2002
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
 
 
#ifndef S_SPECIFIC_FUNCTIONS_H
#define S_SPECIFIC_FUNCTIONS_H

#include <stdio.h>
#include "structs.h"


PARAMTYPE s_dangling_energy (int *sequence, char *structure, int i1, int i2, int i3, int i4);
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling energy between i2 and i3


PARAMTYPE s_dangling_energy_left (int *sequence, char *structure, int i1, int i2, int i3, int i4);
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i1 and i3


PARAMTYPE s_dangling_energy_right (int *sequence, char *structure, int i1, int i2, int i3, int i4);
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2


PARAMTYPE s_dangling_enthalpy (int *sequence, int i1, int i2, int i3, int i4);
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling enthalpy between i2 and i3


PARAMTYPE s_dangling_enthalpy_left (int *sequence, int i1, int i2, int i3, int i4);
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling enthalpy between i1 and i3


PARAMTYPE s_dangling_enthalpy_right (int *sequence, int i1, int i2, int i3, int i4);
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling enthalpy between i1 and i3


PARAMTYPE s_calculate_energy (int *sequence, char *csequence, char *structure, str_features *f, str_features *fres);
// PRE:  the structure features have been determined
// POST: return the free energy (as integer)


PARAMTYPE s_calculate_enthalpy (int *sequence, char *csequence, str_features *f);
// PRE:  the structure features have been determined
// POST: return the enthalpy (as integer)


#endif
