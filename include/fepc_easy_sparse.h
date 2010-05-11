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
 
#ifndef __FEPCEASYSPARSE
#define __FEPCEASYSPARSE

#include "fepc_easy.h"



typedef struct {
	int  summands;														
	int  dimension;											
	func_p * factors;
} func_sparse_t;

typedef func_sparse_t * func_sparse_p;


func_sparse_p 
func_sparse_new(int summands, int dimension);

void
func_sparse_del(func_sparse_p func_sparse);

func_p 
func_sparse_to_func(func_sparse_p);

func_sparse_p
func_sparse_multi(func_sparse_p f, func_sparse_p g, fepc_real_t stepping);

func_sparse_p
func_sparse_reflect(func_sparse_p f);


func_sparse_p
func_sparse_convolute(func_sparse_p f, func_sparse_p g, func_p * result_structure, fepc_real_t stepping);












#endif