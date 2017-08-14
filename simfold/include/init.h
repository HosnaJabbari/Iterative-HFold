/***************************************************************************
                          init.h  -  description
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


#ifndef INIT_H
#define INIT_H
void init_data(const char *arg, char *config_file, int what, double temperature);
// the function that must be called by the main program to read data files
// PRE:  None
// POST: Read all data and configuration files

#endif

