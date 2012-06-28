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
 
#include "fepc_easy_sparse.h"
#include "fepc_easy_helper.h"

func_sparse_p 
func_sparse_new(int summands, int dimension) {
	func_sparse_p result = (func_sparse_p) malloc(sizeof(func_sparse_t));
	
	result->summands = summands;
	result->dimension = dimension;
	result->factors = (func_p*) malloc(sizeof(func_p)*summands*dimension);
	return result;
}


void
func_sparse_del(func_sparse_p func_sparse) {
	int n , i = func_sparse->summands*func_sparse->dimension;
	
	for (n = 0; n < i; n++) {
		func_del(func_sparse->factors[n]);
	}
	free(func_sparse->factors);
	free(func_sparse);
}


func_p 
func_sparse_to_func(func_sparse_p func_sparse) {
	func_p result = func_new(1, func_sparse->dimension);
	
	return result;
}


func_sparse_p
func_sparse_multi(func_sparse_p f, func_sparse_p g, fepc_real_t stepping) {
	int n, i, k, dimension = f->dimension;
	
	func_sparse_p result = func_sparse_new(f->summands * g->summands, dimension);
	
	for (n = 0; n < f->summands; n++) {
		for (i = 0; i < g->summands; i++) {
			for (k = 0; k < dimension; k++) {
				result->factors[(n*g->summands+i)*dimension+k] = func_multi(f->factors[n*dimension+k], g->factors[i*dimension+k], stepping);
			}	
		}	
	}
	return result;	
}

func_sparse_p
func_sparse_reflect(func_sparse_p f) {
	int n, i;
	
	func_sparse_p result = func_sparse_new(f->summands, f->dimension);
	
	for (n = 0, i = f->summands*f->dimension; n < i; n++) {
		result->factors[n] = func_reflect(f->factors[n]);	
	} 
	return result;
}

func_sparse_p
func_sparse_convolute(func_sparse_p f, func_sparse_p g, func_p * result_structure, fepc_real_t stepping) {
	int n, i, k, dimension = f->dimension;
	
	func_sparse_p result = func_sparse_new(f->summands * g->summands, dimension);
	
	for (n = 0; n < f->summands; n++) {
		for (i = 0; i < g->summands; i++) {
			for (k = 0; k < dimension; k++) {
				result->factors[(n*g->summands+i)*dimension+k] = faltung_fepc(f->factors[n*dimension+k], g->factors[i*dimension+k], result_structure[k], stepping);
			}	
		}	
	}
	return result;	
}


