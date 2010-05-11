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
 
#include "discont.h"


int
main(void) {
    /* Example */
    
    fepc_real_t stepping = 0.1;

    int steps, n;

    func_p f1, f2, w, convolution_result;

    steps = 2;

    discont_function_p function, result;

    function = discont_function_new(steps);

    fepc_real_t y11[] = {0.1, 0.0, 0.3, 3.0, 4.2, 0.3, 5.8, 7.3, 1.0, 2.4}; // the stepcount is 10 in the first level because h_0 is 0.1
    fepc_real_t y21[] = {0.2, 1.3, 0.1, 6.5, 4.2, 3.3, 4.1, 0.3, 9.3, 2.2};


    // first step
    discont_function_setup_points(function, 0, 0.0, 1.0, y11, y21, stepping); // interval [0, 1] at level 0

    fepc_real_t y12[] = {3.4, 2.6, 0.3, 4.7, 0.0, 7.5}; // the stepcount is 6 in the second level because h_1 is 0.05
    fepc_real_t y22[] = {9.4, 2.5, 3.0, 1.9, 4.2, 0.5};


    // second step
    discont_function_setup_points(function, 1, 0.2, 0.5, y12, y22, stepping); // interval [0.2, 0.5] at level 1

    discont_function_print(function);
    
    f1 = convert_discont_function(function, stepping);
    discont_function_del(function);

    func_print(f1, 3);
	

    // Create f2 in the same way!


    // set up the structure w of the result

    function = discont_function_new(steps);
    discont_function_setup_points(function, 0, 0.0, 1.0, NULL, NULL, stepping); // interval [0, 1] at level 0
    discont_function_setup_points(function, 1, 0.3, 0.6, NULL, NULL, stepping); // interval [0.3, 0.6] at level 1
    w = func_new(steps-1, 1);

    set_gridstructure(w, function->intervals, stepping);

    convolution_result = faltung_fepc(f1, f2, w, stepping);

    func_del(f1);
    func_del(f2);
    func_del(w);

    result = convert_func(convolution_result, function->intervals, stepping);
    discont_function_del(function);

    func_del(convolution_result);


    // the result is stored in "result" as a discont-function
    
    return 0;
}
