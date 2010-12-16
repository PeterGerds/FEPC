/*
 * FEPC
 * Copyright (C) 2010 Stefan Handschuh (handschu@mis.mpg.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __FEPCCP
#define __FEPCCP


#include "fepc_easy.h"


typedef struct {
	int rank;
	int dimension;
	func_t * functions;
} func_cp;

func_cp *
func_cp_new(int rank, int dimension, int * maxlevels);

void
func_cp_init(func_cp * function, int rank, int dimension, int * maxlevels);


/*func_add
 * Each summand has a different interval
 */
func_cp *
func_cp_new_blockstructure(int rank, int dimension, Funcimpl * functions, interval_t * intervals, int * maxlevels);

/*
 * Each summand has the same interval
 */
func_cp *
func_cp_new_cp(int rank, int dimension, Funcimpl * functions, interval_t * interval, int * maxlevels);

func_t *
func_cp_extract(int current_rank, int current_dimension, func_cp * func_cp);

void
func_cp_del(func_cp * func_cp);

func_cp *
func_cp_faltung(func_cp * function1, func_cp * function2, func_t * resulting_structure, fepc_real_t h);

func_cp *
func_cp_multi(func_cp * function1, func_cp * function2, fepc_real_t h);

void
func_cp_print(func_cp * function);

int *
int_array_new(int length);


#endif // __FEPCCP
