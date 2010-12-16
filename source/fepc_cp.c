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

#include "fepc_cp.h"
#include "fepc_easy_helper.h"

func_cp *
func_cp_new(int rank, int dimension, int * maxlevels) {
	func_cp * result = (func_cp*) malloc(sizeof(func_cp));

	func_cp_init(result, rank, dimension, maxlevels);
	return result;
}

void
func_cp_init(func_cp * function, int rank, int dimension, int * maxlevels) {
	int n, k = rank*dimension;
	function->rank = rank;
	function->dimension = dimension;
	function->functions = (func_t *) malloc(sizeof(func_t)*k);

	for (n = 0; n < k; n++) {
		func_init(&(function->functions[n]), maxlevels[n], 1);
	}
}

func_cp *
func_cp_new_blockstructure(int rank, int dimension, Funcimpl * functions, interval_p intervals, int * maxlevels) {
	func_cp * result = func_cp_new(rank, dimension, maxlevels);


	return result;
}

func_cp *
func_cp_new_cp(int rank, int dimension, Funcimpl * functions, interval_p interval, int * maxlevels) {
	func_cp * result = func_cp_new(rank, dimension, maxlevels);


	return result;
}

func_t *
func_cp_extract(int current_rank, int current_dimension, func_cp * func_cp) {
	ASSERT(current_rank < func_cp->rank);
	ASSERT(current_dimension < func_cp->dimension);
	return &(func_cp->functions[current_rank*func_cp->dimension+current_dimension]);
}

func_cp *
func_cp_faltung(func_cp * function1, func_cp * function2, func_t * resulting_structure, fepc_real_t h) {
	ASSERT(function1->dimension == function2->dimension);

	int * maxlevels = int_array_new(100); // FIXME

	func_cp * result = func_cp_new(function1->rank*function2->rank, function1->dimension, maxlevels);

	int n, k, d;

	for (n = 0; n < function1->rank; n++) {
		for (k = 0; k < function2->rank; k++) {
			for (d = 0; d < result->dimension; d++) {
				faltung_fepc_overwrite(&(result->functions[(n*function2->rank + k)*result->dimension+d]), func_cp_extract(n, d, function1), func_cp_extract(k, d, function2), &resulting_structure[d], h);
			}
		}
	}
	free(maxlevels);
	return result;
}

void
func_cp_del(func_cp * func_cp) {
	funcs_del_type(func_cp->functions, func_cp->dimension*func_cp->rank);
	free(func_cp);
}

func_cp *
func_cp_multi(func_cp * function1, func_cp * function2, fepc_real_t h) {
	ASSERT(function1->dimension == function2->dimension);

	int * maxlevels = int_array_new(100); // FIXME

	func_cp * result = func_cp_new(function1->rank*function2->rank, function1->dimension, maxlevels);

	int n, k, d;

	for (n = 0; n < function1->rank; n++) {
		for (k = 0; k < function2->rank; k++) {
			for (d = 0; d < result->dimension; d++) {
				func_multi_overwrite(&(result->functions[(n*function2->rank + k)*result->dimension+d]), func_cp_extract(n, d, function1), func_cp_extract(k, d, function2), h);
			}
		}
	}
	free(maxlevels);
	return result;
}

void
func_cp_print(func_cp * cp) {
	int n, k;

	printf("CP function with rank %i:\n", cp->rank);
	for (n = 0; n < cp->rank; n++) {
		for (k = 0; k < cp->dimension; k++) {
			func_print(func_cp_extract(n, k, cp), 3);
		}
		if (n +1 < cp->rank) {
			printf("\n+");
		} else {
			printf("\n\n");
		}
	}
}

int *
int_array_new(int length) {
	int n;

	int * result = (int*) malloc(sizeof(int)*length);

	for (n = 0; n < length; n++) {
		result[n] = 0;
	}
	return result;
}


