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


fepc_real_t 
_round(fepc_real_t value) {
    if (value - ((int) value) < 0.5) {
        return (fepc_real_t) (int) value;
    } else {
        return (fepc_real_t) (int) value +1.0;
    }
}

linear_function_p
linear_function_new(fepc_real_t y_0, fepc_real_t slope) {
    linear_function_p result = (linear_function_p) malloc(sizeof(linear_function_t));
    
    result->y_0 = y_0;
    result->slope = slope;
    return result;
}


linear_function_p
linear_function_new_points(fepc_real_t y_0, fepc_real_t y_1, fepc_real_t x_0, fepc_real_t x_1) {
	fepc_real_t slope = (y_1-y_0) / (x_1-x_0);

    return linear_function_new(y_0-slope*x_0, slope);
}


linear_function_set_p
linear_function_set_new(int count) {
    linear_function_set_p result = (linear_function_set_p) malloc(sizeof(linear_function_set_t));
    
    result->count = count;
    result->functions = (linear_function_p*) malloc(sizeof(linear_function_p)*count);
    return result;
}


discont_function_p
discont_function_new(int stepcount) {
    discont_function_p result = (discont_function_p) malloc(sizeof(discont_function_t));
    
    result->stepcount = stepcount;
    result->intervals = (interval_p*) malloc(sizeof(interval_p)*stepcount);
    result->function_sets = (linear_function_set_p*) malloc(sizeof(linear_function_set_p)*stepcount);
    return result;
}

void
discont_function_del(discont_function_p function) {
    int n;
    
    for (n = 0; n < function->stepcount; n++) {
        interval_del(function->intervals[n]);
        linear_function_set_del(function->function_sets[n]);
    }
    free(function->intervals);
    free(function->function_sets);
    free(function);
}

void
linear_function_set_del(linear_function_set_p function_set) {
    int n;
    

    for (n = 0; n < function_set->count; n++) {
    	if (function_set->functions[n] != NULL) {
    		free(function_set->functions[n]);
    	}
    }
    free(function_set->functions);
    free(function_set);
}

void
discont_function_setup(discont_function_p function, int step, fepc_real_t start, fepc_real_t end, linear_function_set_p function_set) {
    function->intervals[step] = interval_new(1);
    function->intervals[step]->start[0] = start;
    function->intervals[step]->end[0] = end;
    function->function_sets[step] = function_set;
}

void
discont_function_setup_points(discont_function_p function, int step, fepc_real_t start, fepc_real_t end, fepc_real_t * y1, fepc_real_t * y2, fepc_real_t stepping) {
    function->intervals[step] = interval_new(1);
    function->intervals[step]->start[0] = start;
    function->intervals[step]->end[0] = end;
    
    fepc_real_t h_l = get_h_l(step, stepping);
    
    int count, n;
    
    count = _round((end-start)/h_l); // this has to be equal to the array sizes of y1 and y2

    function->function_sets[step] = linear_function_set_new(count);

    if (y1 != NULL && y2 != NULL) {
        for (n = 0; n < count; n++) {
            function->function_sets[step]->functions[n] = linear_function_new_points(y1[n], y2[n], start+n*h_l, start+(n+1)*h_l );
        }
    } else {
    	for (n = 0; n < count; n++) {
            function->function_sets[step]->functions[n] = NULL;
        }
    }
}


void
discont_function_print(discont_function_p function) {
    int n;
    
    printf("Discontinuous function with %i levels\n=====================================\n", function->stepcount);
    
    for (n = 0; n < function->stepcount; n++) {

        printf("Level %i\n", n);
        print_interval(function->intervals[n]);
        linear_function_set_print(function->function_sets[n]);
    }
}


void
linear_function_set_print(linear_function_set_p function_set) {
    int n;
    
    printf("Linear function set\n=====\n");
    for (n = 0; n < function_set->count; n++) {
        printf("   %f*x + %f\n", function_set->functions[n]->slope, function_set->functions[n]->y_0);
    }
}


func_p 
convert_discont_function(discont_function_p function, fepc_real_t stepping) {
    func_p result = func_new(function->stepcount-1, 1);
    
    set_degree(result, 1); // sets the polynomial degree to 1    
    set_gridstructure(result, function->intervals, stepping); // set the correct intervals
    add_folgenentries_discont(result, function, stepping); // add folgenentries
    return result;
}


discont_function_p
convert_func(func_p function, interval_p * intervals, fepc_real_t stepping) {
    discont_function_p result = discont_function_new(function->maxlevel+1);
    
    fepc_real_t h_l, temp_x1, temp_x2, temp_y1, temp_y2, slope;
    
    int n, k, stepcount, start;
    
    vec_real_p x;
    
    for (n = 0; n < result->stepcount; n++) {
        h_l = get_h_l(n, stepping);
        result->intervals[n] = intervals[n];
        stepcount = _round((intervals[n]->end[0] - intervals[n]->start[0])/h_l);
        result->function_sets[n] = linear_function_set_new(stepcount);
        start = function->hierarchie[n]->vektor[0]->start->array[0];
        for (k = 0; k < stepcount; k++) {
            temp_x1 = (ONE_THIRD+start+k)*h_l;
            temp_x2 = (TWO_THIRD+start+k)*h_l;
            x = vec_real_new(1);
            x->array[0] = temp_x1;
            temp_y1 = get_value_at_step(function, x, n, stepping);
            x->array[0] = temp_x2;
            temp_y2 = get_value_at_step(function, x, n, stepping);
            vec_real_del(x);
            slope = 3*(temp_y2-temp_y1)/h_l;
            result->function_sets[n]->functions[k] = linear_function_new(temp_y1 - slope*temp_x1, slope);
        }
    }
    return result;
}


fepc_real_t
integrate_coeff_discont(discont_function_p function, int position, int v, int p, int step, fepc_real_t stepping) {
    // integrate from v[0]*h_l till (v[0]+1)*h_l
	
    fepc_real_t h_l, slope, y_0, sqrt_h_l;
    h_l = get_h_l(step, stepping);
	
    slope = function->function_sets[step]->functions[position]->slope;
    y_0 = function->function_sets[step]->functions[position]->y_0;

    sqrt_h_l = sqrt(h_l);
    if (p == 0) { // legendre(0) = sqrt(1/h_l)
        return sqrt_h_l*(y_0 + h_l*(slope/2.0));
    } else { // p == 1 --> legendre(1) = sqrt(12)(x-(v+0.5)*h_l)/(h_l^1.5)
        return h_l*sqrt_h_l*slope/SQRT_12;
    }
}


void
add_folgenentries_discont(func_p function, discont_function_p discont_function, fepc_real_t stepping) {

    int step, n, k, interval_element_count, position, grad_count, pos_2;

    vec_p v, x;

    for (step = 0; step <= function->maxlevel; step++) {
        grad_count = get_degree_count(function->hierarchie[step]->grad);
        for (k = 0; k < grad_count; k++) {
            for (n = 0, interval_element_count = get_interval_element_count(function->hierarchie[step]->vektor[k]); n < interval_element_count; n++) {
                v = generate_folgenvector(function->hierarchie[step]->vektor[k], n);
                pos_2 = v->array[0];
                position = entry_d2one(v, function->hierarchie[step]->vektor[k]->lang);
                vec_add2(v, function->hierarchie[step]->vektor[k]->start); // needed bc v starts from 0
                if (step == function->maxlevel || !is_in_latter_interval(v, function->hierarchie[step+1]->vektor[k])) { // only set the vector if it is in the valid interval
                    x = entry_one2d_sloppy(k, function->hierarchie[step]->grad); // don't care about 0 as entries
                    function->hierarchie[step]->vektor[k]->glied[position] = integrate_coeff_discont(discont_function, pos_2 ,v->array[0], x->array[0], step, stepping);
                    vec_del(x);
                }
                vec_del(v);
            }
        }
    }
}


