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

fepc_real_t
id(vec_real_p vector) {
	return vector->array[0];
}

int
main() {
	fepc_real_t h = 0.2;

	interval_p * intervals = (interval_p*) malloc(sizeof(interval_p)*1);

	intervals[0] = interval_new(1);
	intervals[0]->start[0] = -1.;
	intervals[0]->end[0] = 1.;

	int rank1 = 3, rank2 = 4, dimension = 2;

	int n, k;

	func_cp_p cp1 = func_cp_new(rank1, dimension);

	for (n = 0; n < rank1; n++) {
		for (k = 0; k < dimension; k++) {
			cp1->functions[n*dimension+k] = setup_fepc_structure(id, intervals, 1, 0, h);
		}
	}

	func_cp_p cp2 = func_cp_new(rank2, dimension);

	for (n = 0; n < rank2; n++) {
		for (k = 0; k < dimension; k++) {
			cp2->functions[n*dimension+k] = setup_fepc_structure(id, intervals, 1, 0, h);
		}
	}

	func_p * result_intervals = (func_p*) malloc(sizeof(func_p)*dimension);

	for (k = 0; k < dimension; k++) {
		result_intervals[k] = func_new(0, 1);
		set_gridstructure(result_intervals[k], intervals, h);
	}

	func_cp_p result = func_cp_faltung(cp1, cp2, result_intervals, h);


	func_cp_del(cp1);
	func_cp_del(cp2);

	func_cp_print(result);

	funcs_del(result_intervals, dimension);
	intervals_del(intervals, 1);
	func_cp_del(result);


	return 0;
}
