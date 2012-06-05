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

	interval_t * intervals = (interval_t*) malloc(sizeof(interval_t)*1);

	interval_init(&intervals[0], 1);
	intervals[0].start[0] = -1.;
	intervals[0].end[0] = 1.;

	int rank1 = 3, rank2 = 4, dimension = 2;

	int n, k, interval_count = 1;

	int * maxlevels_cp1 = int_array_new(rank1*dimension);

	func_cp * cp1 = func_cp_new(rank1, dimension, maxlevels_cp1);

	free(maxlevels_cp1);

	for (n = 0; n < rank1; n++) {
		for (k = 0; k < dimension; k++) {
			setup_fepc_structure(&(cp1->functions[n*dimension+k]), id, &intervals, 1, 0, h);
		}
	}

	int * maxlevels_cp2 = int_array_new(rank2*dimension);

	func_cp * cp2 = func_cp_new(rank2, dimension, maxlevels_cp2);
	free(maxlevels_cp2);
	for (n = 0; n < rank2; n++) {
		for (k = 0; k < dimension; k++) {
			setup_fepc_structure(&(cp2->functions[n*dimension+k]), id, &intervals, 1, 0, h);
		}
	}

	func_t * result_intervals = (func_t*) malloc(sizeof(func_t)*dimension);

	for (k = 0; k < dimension; k++) {
		func_init(&result_intervals[k], 0, 1);
		set_gridstructure(&result_intervals[k], &intervals, h);
	}

	func_cp * result = func_cp_faltung(cp1, cp2, result_intervals, h);


	func_cp_del(cp1);
	func_cp_del(cp2);

	func_cp_print(result);

	funcs_del_type(result_intervals, dimension);
	intervals_del_type(intervals, 1);
	func_cp_del(result);
#ifdef HAS_FFTW3
	fftw_cleanup();
#endif
	return 0;
}
