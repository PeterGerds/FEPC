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
 
#include "fepc_easy.h"


/**
 * This method converts standard d-dimensional intervals to the form that is needed by the fepc algorithm.
 * Only called for debug. The transformation is also done inside of set_gridstructure.
 * 
 * @param steps     defines the number of steps that are represented by the invervals; equals the number of intervals
 * @param intervals array of intervals that should be converted
 */
interval_p * convert_interval(int steps, interval_p * intervals, fepc_real_t stepping) {
    interval_p* result = (interval_p*) malloc(sizeof(interval_p) * steps);
    
    fepc_real_t pow, h_l;
    
    pow = 2.;
    
    int n, k;
    
    for (n = 0; n < steps; n++) {
        pow *= 0.5;
        h_l = pow*stepping;
        result[n] = interval_new(intervals[n]->dimension);
        for (k = 0; k < result[n]->dimension; k++) {
            result[n]->start[k] = intervals[n]->start[k] / h_l; 
            result[n]->end[k] = intervals[n]->end[k] / h_l - 1; 
        }
    }
    
    return result;
}

/**
 * Prints out all elements of the set and its results.
 * 
 * @param set set to print out
 */
void print_vec_real_set(vec_real_set_p set) {
    int n, k;
    printf("   Elements \n"); 
    for (n = 0; n < set->count; n++) {
        for (k = 0; k < set->elements[n]->dim; k++) {
            printf("%.5f  ", set->elements[n]->array[k]);
        }
        printf("%.5f\n", set->element_results[n]);
    }
}

/**
 * Creates a new set with allocated memory for the entries.
 *
 * @param count number of entries
 *
 * @return new set
 */
vec_real_set_p vec_real_set_new(int count) {
    vec_real_set_p result = (vec_real_set_p) malloc(sizeof(vec_real_set_t));
    
    result->count = count;
    result->elements = (vec_real_p*) malloc(sizeof(vec_real_p)*count);
    result->element_results = (fepc_real_t*) malloc(sizeof(fepc_real_t)*count);
    
    int n;
    
    for (n = 0; n < count; n++) {
        result->element_results[n] = 0.0;
    }
    return result;
}


void vec_real_set_del(vec_real_set_p set) {
    int n;
    
    for (n = 0; n < set->count; n++) {
        vec_real_del(set->elements[n]);
    }
    free(set->elements);
    free(set->element_results);
    free(set);
}


fepc_real_t get_h_l(int step, fepc_real_t stepping) {
    return pow(.5, step)*stepping;
}

/**
 * Returns the value of the legendre-polynomial with the given degree at the given point.
 *
 * @param degree degree of the legendre-polynomial
 * @param x      point where the polynomial should be evaluated
 *
 * @return P_degree(x)
 */
fepc_real_t legendre(int degree, fepc_real_t x) {
    if (degree == 0) {
        return 1;
    } else if (degree == 1) {
        return x;
    } else {
        return ((2*degree-1)*x*legendre(degree-1, x) - (degree-1)*legendre(degree-2, x)) / degree;
    }
}

/**
 * This method calculates the basis function.
 *
 * @param step current step
 * @param v    vector that describes the integral-bounds; integral is taken from v[n]*h_l upto (v[n]+1)*h_l
 * @param p    degree-vector of the basis function
 * @param x    point at which the basis function is evaluated
 *
 * @return phi_{v,p}(x)
 */
fepc_real_t phi_l(int step, vec_p v, vec_p p, vec_real_p x, fepc_real_t stepping) {
    fepc_real_t h_l, result;
    
    h_l = get_h_l(step, stepping);
    result = 1.;
    
    int n;
    
    for (n = 0; n < v->dim; n++) { 
        if (v->array[n]*h_l > x->array[n] || (v->array[n]+1)*h_l < x->array[n]) {
            return 0;
        }
        result *= sqrt((2*p->array[n] + 1)/h_l)*legendre(p->array[n], 2.0*x->array[n]/h_l - 2.0*v->array[n] - 1.0); 
    } 
    return result;
}

/**
 * Retransforms the resulting points from S(M) into R^n.
 * This method should be used at the end of the calculation to regain the actual R^n coodinates of the resulting points.
 * 
 * @param function        function in S(M); this should be a result of fecp
 * @param generate_points a function that generates the points where the R^n values should be taken
 *
 * @return a set of all valid points
 */
vec_real_set_p get_value(func_p function, Funcimpl_vec_step generate_points, fepc_real_t stepping) {   
    vec_p r, x;
    
    int n, k, m, l, intervals_element_count, position, degree_count;
    
    int* interval_element_counts = (int*) malloc(sizeof(int)*(function->maxlevel+1));
    
    // get the total count of interval elements as well as the separate interval lengths
    for (n = 0, intervals_element_count = 0; n <= function->maxlevel; n++) {
        interval_element_counts[n] = get_interval_element_count(function->hierarchie[n]->vektor[0]); // intervals have all the same length, so we can use any "vektor"
        intervals_element_count += interval_element_counts[n];
    }
    vec_real_set_p result = vec_real_set_new(intervals_element_count);
    for (n = 0, m = 0; n <= function->maxlevel; n++) {
        degree_count = get_degree_count(function->hierarchie[n]->grad);
        for (k = 0; k < interval_element_counts[n]; k++, m++) {
            for (l = 0; l < degree_count; l++) {
                r = generate_folgenvector(function->hierarchie[n]->vektor[l], k);
                position = entry_d2one(r, function->hierarchie[n]->vektor[l]->lang);
                vec_add2(r, function->hierarchie[n]->vektor[l]->start);           
                // only take the point if it is inside of the correct A-interval; if not in interval -> reduce size of set
                if (n == function->maxlevel || !is_in_latter_interval(r, function->hierarchie[n+1]->vektor[l])) {
                    x = entry_one2d_sloppy(l, function->hierarchie[n]->grad);
                    result->elements[m] = generate_points(r, x, n, stepping);
                    result->element_results[m] += function->hierarchie[n]->vektor[l]->glied[position]* phi_l(n, r, x, result->elements[m], stepping);
                    vec_del(x);
                } else {
                    m--;
                    result->count--;
                    // TODO free also memory?
                }
                vec_del(r);
            } 
        }
    }
    free(interval_element_counts);
    return result;
}

fepc_real_t get_value_at_step(func_p function, vec_real_p x, int step, fepc_real_t stepping) {
    ASSERT(function != NULL && step <= function->maxlevel);
    
    int n, l, k, degree_count, interval_element_count, position;
    
    vec_p r, p;
    
    fepc_real_t result = 0.;
    
    degree_count = get_degree_count(function->hierarchie[step]->grad);
    interval_element_count = get_interval_element_count(function->hierarchie[step]->vektor[0]); // intervals have all the same length, so we can use any "vektor"
    for (n = 0; n < interval_element_count; n++) {
        for (k = 0; k < degree_count; k++) {
            r = generate_folgenvector(function->hierarchie[step]->vektor[k], n);
            position = entry_d2one(r, function->hierarchie[step]->vektor[k]->lang);
            vec_add2(r, function->hierarchie[step]->vektor[k]->start);           
            // only take the point if it is inside of the correct A-interval
            if (step == function->maxlevel || !is_in_latter_interval(r, function->hierarchie[step+1]->vektor[k])) {
                p = entry_one2d_sloppy(k, function->hierarchie[step]->grad);
                result += function->hierarchie[step]->vektor[k]->glied[position]* phi_l(step, r, p, x, stepping);
                vec_del(p);
            }
            vec_del(r);
        } 
    }
    return result;

}

fepc_real_t get_value_at(func_p function, vec_real_p x, fepc_real_t stepping) {
    int n;
    
    for (n = 0; n < function->maxlevel; n++) {
        if (!is_in_latter_interval_real(x, function->hierarchie[n+1]->vektor[0], n+1, stepping)) { // the intervals are at the level constant, independent of the degree        
            return get_value_at_step(function, x, n, stepping);//*(function->maxlevel - n +1);
        }
    }
    return get_value_at_step(function, x, function->maxlevel, stepping);
}

double round(double x) {
    fepc_real_t rest = x - (int) x;
    
    if (rest < 0.5) {
        return (int) x;
    } else {
        return 1 + (int) x;
    }    
}

/**
 * Sets the correct grid-structure that is needed by the fepc algorithm out of common human-readable intervals.
 *
 * @param function  function whose grid structure should be set
 * @param intervals human-readable intervals; note the the inverval count has to be function->maxlevel+1
 */
void set_gridstructure(func_p function, interval_p* intervals, fepc_real_t stepping) {
    vec_p start, lang;
   
    int n, dimension, k, steps, grad_count;

    steps = function->maxlevel+1;
    dimension = function->dim;
    
    fepc_real_t pow, h_l;
    
    pow = 2.;    
    
    folge_p folge;
    
    for (n = 0; n < steps; n++) {
        start = vec_new(dimension);
        lang = vec_new(dimension);
        
        pow *= 0.5;
        h_l = pow*stepping;
            
        if (intervals[n] != NULL)
        for (k = 0; k < dimension; k++) { // here the interval is converted to the nessesary form
            start->array[k] = intervals[n]->start[k] / h_l;
            lang->array[k] = round(intervals[n]->end[k] / h_l) - start->array[k];  // lang = length (integer), see p. 68; round() is needed to correct error
        }

        folge = folge_new(start, lang);
        folge_del(function->hierarchie[n]->vektor[0]);
        function->hierarchie[n]->vektor[0] = folge;
        for (k = 1, grad_count = get_degree_count(function->hierarchie[n]->grad); k < grad_count; k++) {
            folge_del(function->hierarchie[n]->vektor[k]);
            function->hierarchie[n]->vektor[k] = folge_copy(folge);
        }

        #ifdef __DEBUG
        printf("set grid structure\n  start\n");
        print_vec(start);
        printf("  lang\n");
        print_vec(lang);
        #endif
    }
}

fepc_real_t sqr(fepc_real_t x) {
    return x*x;
}

fepc_real_t integrate_legendre(vec_p v, vec_p degree, int step, fepc_real_t stepping) {
    int n;
    
    fepc_real_t h_l = get_h_l(step, stepping), result = 1.;
    
    for (n = 0; n < v->dim; n++) {
        if (degree->array[n] == 0) {
            result *= sqrt(h_l);
        } else {
            result *= sqrt(h_l/(4*(2*degree->array[n] + 1)))*(legendre(degree->array[n] + 1, -1) - legendre(degree->array[n]-1, -1));
        }
        //result *= sqrt((2*degree->array[n] + 1)/h_l)*h_l;//(legendre(degree->array[n], sqr((v->array[n]+1)*h_l)/2.) - legendre(degree->array[n], sqr(v->array[n]*h_l)/2.));
    }
}


fepc_real_t func_integrate(func_p function, fepc_real_t stepping) {
    int step, n, k, interval_element_count, position, grad_count;
    
    vec_p v, x;
    
    fepc_real_t result = 0.;
    
    for (step = 0; step <= function->maxlevel; step++) {
        grad_count = get_degree_count(function->hierarchie[step]->grad);
        for (k = 0; k < grad_count; k++) {
            for (n = 0, interval_element_count = get_interval_element_count(function->hierarchie[step]->vektor[k]); n < interval_element_count; n++) {
                v = generate_folgenvector(function->hierarchie[step]->vektor[k], n);
                position = entry_d2one(v, function->hierarchie[step]->vektor[k]->lang);
                vec_add2(v, function->hierarchie[step]->vektor[k]->start); // needed bc v starts from 0 
                if (step == function->maxlevel || !is_in_latter_interval(v, function->hierarchie[step+1]->vektor[k])) { // only use the vector if it is in the valid interval
                    x = entry_one2d_sloppy(k, function->hierarchie[step]->grad); // don't care about 0 as entries
                    result += function->hierarchie[step]->vektor[k]->glied[position]*integrate_legendre(v, x, step, stepping);                   
                    vec_del(x);
                }            
                vec_del(v);
            }
        }
    }
}

/**
 * Returns the number of discrete elmenets that fit in the given folge.
 *
 * @param folge the folge whose elements should be count
 *
 * @return #folge
 */
int get_interval_element_count(folge_p folge) {
    int n, result;
    
    for (n = 0, result = 1; n < folge->lang->dim; n++) {
        result *= folge->lang->array[n];
    }
    return result;
}

/**
 * This method generates the vector at the given position inside of the inverval that is defined in the given folge.
 *
 * @param folge            the folge that defines the interval
 * @param element_position position of the vector that is requested
 *
 * @return folge.interval(element_position)
 */
vec_p generate_folgenvector(folge_p folge, int element_position) {
    vec_p result = vec_new(folge->lang->dim);
    
    int n, length;
    
    for (n = 0; n < result->dim; n++) {
        length = folge->lang->array[n];
        result->array[n] = element_position % length;
        element_position = (element_position - result->array[n]) / length;
        //result->array[n] += folge->start->array[n]; that should be here but in that case the algorithm does not work --> use vec_add2 in add_folgenentries with the start-value
    }
    return result;
}

vec_real_p generate_x_for_integration(vec_p v, int calc_position, fepc_real_t h, fepc_real_t h_l, int stepcount) {
    vec_real_p result = vec_real_new(v->dim);
    
    int n, current;
       
    for (n = 0; n < v->dim; n++) {
        current = calc_position % stepcount;
        calc_position = (calc_position - current) / stepcount;
        result->array[n] = v->array[n]*h_l+current*h;
    }
    return result;
}

bool_t has_value_until(int dimension, int calc_position, int stepcount, int start, int end) {
    int n, current;
       
    for (n = 0; n < dimension; n++) {
        current = calc_position % stepcount;
        if (current >= start && current <= end) {
            return true;
        }
        calc_position = (calc_position - current) / stepcount;
    }
    return false;
}

int get_current_changing_dimension(int dimension, int calc_position, int stepcount) {
    int n, current;
    
    for (n = 0; n < dimension; n++) {
        current = calc_position % stepcount;
        if (current != calc_position -1 % stepcount) {
            return n;
        }
        calc_position = (calc_position - current) / stepcount;
    }
    return 0;

}

fepc_real_t integrate_coeff(Funcimpl function_impl, vec_p v, vec_p p, int step, fepc_real_t stepping) {
    int i, n, calc_count, position;
    
    n = INT_STEPS;
    
    fepc_real_t h, h_l, phi, f_x, result, k;
    
    bool_t second_sum = false;    
    h_l = get_h_l(step, stepping);
    
    
    h = h_l/n; // step size
    
    result = 0.;
    
    
    // interval length is always h_l! (by design)
    
    calc_count = pow(n+1, v->dim);
    
    vec_real_p x; 
    #ifdef __DEBUG
    printf(" calc_stepcount = %i\n", calc_count);
    #endif

    for (i = 0; i < calc_count; i++) {
        x = generate_x_for_integration(v, i, h, h_l, n+1);
        f_x = function_impl(x);
        phi = phi_l(step, v, p, x, stepping);
        if (has_value_until(v->dim, i, n+1, 0, 0)) { // first component
            result += 0.5*f_x*phi;
        }
        if (has_value_until(v->dim, i, n+1, n, n)) { // last component
            result += 0.5*f_x*phi;
        }
        if (has_value_until(v->dim, i, n+1, 1, n-1)) { // first sum
            result += f_x*phi;
            second_sum = true;
        }
        if (second_sum || has_value_until(v->dim, i, n+1, 1, n)) { // second sum
            position = i % v->dim;
            k = (x->array[position] - v->array[n]*h_l)/h;
            x->array[position] = (x->array[position] + v->array[n]*h_l+(k-1)*h)/2.;// only one coordinate has to be changed
            result += 2*function_impl(x)*phi_l(step, v, p, x, stepping);
            vec_real_del(x);
        }
        vec_real_del(x);
        second_sum = false;        
    }
    #ifdef __DEBUG
    printf("sum = %f, result = %f", result, pow(h_l, v->dim)*result / (2.*calc_count));
    #endif
    return pow(h/3, v->dim)*result;
}

/**
 * Generates the real-valued vector at the calc_position.
 *
 * @param v vector      that descripes the interval that is currently used
 * @param calc_position position in the interval of which the vector should be returned
 * @param h             step size
 * @param h_l           h_l
 * @param stepcount     total number of interval elements
 *
 * @return vector at calc_position
 */
vec_real_p generate_x_for_integration_st(vec_p v, int calc_position, fepc_real_t h, fepc_real_t h_l, int stepcount) {
    vec_real_p result = vec_real_new(v->dim);
    
    int n, current;
       
    for (n = 0; n < v->dim; n++) {
        current = calc_position % stepcount;
        calc_position = (calc_position - current) / stepcount;
        result->array[n] = v->array[n]*h_l+(2.0*current +1)*h/2.0;
    }
    return result;
}

/**
 * Returns the number of start or end values of the given vector in the calculated interval.
 * @param v vector      that descripes the interval that is currently used
 * @param calc_position position in the interval of which the vector should be returned
 * @param h             step size
 * @param h_l           h_l
 * @param stepcount     total number of interval elements
 *
 * @return startvalues + endvalues
 */
int has_start_or_end_value(vec_p v, int calc_position, fepc_real_t h, fepc_real_t h_l, int stepcount) { 
    int n, current, result;
    
    result = 0;
    for (n = 0; n < v->dim; n++) {
        current = calc_position % stepcount;
        
        calc_position = (calc_position - current) / stepcount;
        if (current == 0 || current == stepcount-1) {
            //result = 1
            result++;
        }
    }
    //printf("%i\n", result);
    return result;
}

/**
 * Integrates the given function via trapezoidal rule. Uses multiple cores where available!
 *
 * @param function_impl function to integrate
 * @param v             vector that defines the integration-bounds (v[n]*h_l till (v[n]+1)*h_l)
 * @param p             current polynomial degree
 * @param step          step
 *
 * @return integration result
 */
fepc_real_t integrate_coeff_st(Funcimpl function_impl, vec_p v, vec_p p, int step, fepc_real_t stepping) {
    // uses Sehnen-Trapez
    int n, calc_count;
    
    fepc_real_t h, result, h_l, temp;
    h_l = get_h_l(step, stepping);
    
    h = h_l/(INT_STEPS); // step size
    result = 0.0;
    // interval length is always h_l! (by design)
    
    calc_count = pow(INT_STEPS, v->dim);
    
    vec_real_p x; 
    
    //printf(" calc_stepcount = %i\n", calc_count);
    
    
    #pragma omp parallel for private(temp) reduction(+:result) schedule(static,1)
    for (n = 0; n < calc_count; n++) {
        x = generate_x_for_integration_st(v, n, h, h_l, INT_STEPS);
        //printf("%f\n", phi_l(step, v, p, x, stepping));
        temp = function_impl(x)*phi_l(step, v, p, x, stepping)*pow(.5,has_start_or_end_value(v, n, h, h_l, INT_STEPS));
        result += temp;
        vec_real_del(x);
    }
    
    return pow(h_l, v->dim)*result / pow(INT_STEPS-1, v->dim);
}


/**
 * Returns true if the given vector is in the latter interval that is defined by the given folge.
 *
 * @param v     vector whose position is going to be checked
 * @param folge defines the interval
 *
 * @return v in folge.interval 
 */
bool_t is_in_latter_interval(vec_p v, folge_p folge) {
    int n;
    
    for (n = 0; n < v->dim; n++) {
        if (v->array[n] < folge->start->array[n] / 2. || v->array[n] >= (folge->start->array[n]+folge->lang->array[n]) /2.) {
            return false;
        }
    }
    
    return true;
}

/**
 * Returns true if the given vector is in the latter interval that is defined by the given folge.
 *
 * @param v     vector whose position is going to be checked
 * @param folge defines the interval
 *
 * @return v in folge.interval 
 *
 * ACTS DIFFERENTLY TO is_in_latter_interval !
 */
bool_t is_in_latter_interval_real(vec_real_p v, folge_p folge, int step, fepc_real_t stepping) {
    int n;
    
    fepc_real_t h_l = get_h_l(step, stepping);
    for (n = 0; n < v->dim; n++) {//printf("Start = %i Laenge = %i\n",folge->start->array[n], folge->lang->array[n]);
        if (v->array[n] < h_l*folge->start->array[n]  || v->array[n] >= h_l*(folge->start->array[n]+folge->lang->array[n])) {
            return false;
        }
    }
    return true;
}

/**
 * Returns the number of degree-vectors that are possible in [0, degree).
 *
 * @param degree vector that contains the component-degrees
 *
 * @return |[0, degree)|
 */
int get_degree_count(vec_p degree) {
    int result, n;
    
    result = 1;
    for (n = 0; n < degree->dim; n++) {
        result *= degree->array[n]+1;    
    }
    return result;
}

/**
 * Sets all entries of the hp structure of the given function.
 *
 * @param function the function whose entries are goin to be set
 * @param function_impl R^n->R function from the convolution; maybe NULL if coeff_function is specified
 * @param coeff_function function that calculate the hp-structure as a closed expression without numerical integration; maybe NULL if function_impl is specified; preferred 
 */
void add_folgenentries(func_p function, Funcimpl function_impl, Funcimpl_step coeff_function, fepc_real_t stepping) {
    ASSERT(function_impl != NULL || coeff_function != NULL); // we need one of the two functions to generate the folge
    
    int step, n, k, interval_element_count, position, grad_count;
    
    vec_p v, p;
    
    for (step = 0; step <= function->maxlevel; step++) {
        grad_count = get_degree_count(function->hierarchie[step]->grad);
        for (k = 0; k < grad_count; k++) {
            for (n = 0, interval_element_count = get_interval_element_count(function->hierarchie[step]->vektor[k]); n < interval_element_count; n++) {
                v = generate_folgenvector(function->hierarchie[step]->vektor[k], n);
                position = entry_d2one(v, function->hierarchie[step]->vektor[k]->lang);
                vec_add2(v, function->hierarchie[step]->vektor[k]->start); // needed bc v starts from 0 
                if (step == function->maxlevel || !is_in_latter_interval(v, function->hierarchie[step+1]->vektor[k])) { // only set the vector if it is in the valid interval
                    p = entry_one2d_sloppy(k, function->hierarchie[step]->grad); // don't care about 0 as entries
                    // we prefer the coeff-function due to speed and acuracy:
                    function->hierarchie[step]->vektor[k]->glied[position] = coeff_function != NULL ? coeff_function(v, p, step, stepping) : integrate_coeff_st(function_impl, v, p, step, stepping);
                    vec_del(p);
                }            
                vec_del(v);
            }
        }
    }
}

/**
 * Sets the degree to every element in the given function.
 *
 * @param function the function where to set the degree
 * @param degree   degree
 */
void set_degree(func_p function, int degree) {
    int n, dimension;
    
    dimension = function->dim;
    vec_p p = vec_new(dimension);
    
    for (n = 0; n < dimension; n++) {
        p->array[n] = degree;
    }
    for (n = 0; n <= function->maxlevel; n++) {
        folgen_vektor_del(function->hierarchie[n]);
        function->hierarchie[n] = folgen_vektor_new(vec_copy(p));
    }
    vec_del(p);
}



/**
 * Sets up the needed structure from human readable intervals.
 *
 * @param func           function that should be set up
 * @param function       R^n -> R function
 * @param intervals      human readable inverals
 * @param interval_count number of intervals
 * @param degree         degree of the function (0 for piecewise constant)
 * @param stepping       stepping
 *
 */
void
setup_fepc_structure(func_t * func, Funcimpl function, interval_p* intervals, int interval_count, int degree, fepc_real_t stepping) {
    ASSERT(func->dim == intervals[interval_count-1]->dimension);
    ASSERT(func->maxlevel == interval_count-1);
	if (degree > 0) {
        set_degree(func, degree);
    }
    
    set_gridstructure(func, intervals, stepping);
    if (function != NULL) {
        add_folgenentries(func, function, NULL, stepping);
    }
}

/**
 * Creates a functions and sets up the correct fepc structure.
 *
 * @param function       R^n -> R function
 * @param intervals      human readable inverals
 * @param interval_count number of intervals
 * @param degree         degree of the function (0 for piecewise constant)
 * @param stepping       stepping
 *
 * @return setted up function
 */

func_t *
create_fepc_structure(Funcimpl function, interval_p* intervals, int interval_count, int degree, fepc_real_t stepping) {
	func_t * result = func_new(interval_count-1, intervals[interval_count-1]->dimension);

	setup_fepc_structure(result, function, intervals, interval_count, degree, stepping);
	return result;
}

fepc_real_t max(fepc_real_t a, fepc_real_t b) {
    return a > b ? a : b;
}

bool_t isSuccessor(vec_p vector1, vec_p vector2, int direction, int delta) {
    int n;
    
    for (n = 0; n < vector1->dim; n++) {
        if (n == direction) {
            if (vector2->array[n] - vector1->array[n] != delta) {
                return false;
            }
        } else {
            if (vector2->array[n] != vector1->array[n]) {
                return false;
            }
        }
    }
    return true;
}

vec_p getSuccessor(vec_p v, int position, folge_p folge, int direction, int element_count) {
    int n;
    
    vec_p result;
    
    for (n = position+1; n < element_count; n++) {
        result = generate_folgenvector(folge, n);
        if (isSuccessor(v, result, direction, 1)) {
            return result;
        } else {
            vec_del(result);
        }
    }
    return NULL;
}

vec_p getPredeccessor(vec_p v, int position, folge_p folge, int direction) {
    int n;
    
    vec_p result;
    
    for (n = position-1; n >-1; n--) {
        result = generate_folgenvector(folge, n);
        if (isSuccessor(v, result, direction, -1)) {
            return result;
        } else {
            vec_del(result);
        }
    }
    return NULL;
}

/**
 * Derives a function in the given direction.
 *
 * @param function  function to derive
 * @param direction direction of the derivative
 * @param stepping  stepping
 *
 * @return derived function in the direction
 */
func_p func_derive(func_p function, int direction, fepc_real_t stepping) {
    ASSERT(function->dim > direction);

    func_p result = func_new(function->maxlevel, function->dim);
    
    vec_p v, v2;
        
    int step, n, n2, i, grad_count, position, position2, interval_element_count, position3;
    
    for (step = 0; step <= function->maxlevel; step++) {
        result->hierarchie[step]->vektor[0] = folge_new(vec_copy(function->hierarchie[step]->vektor[0]->start), vec_copy(function->hierarchie[step]->vektor[0]->lang));
        for (n = 0, interval_element_count = get_interval_element_count(function->hierarchie[step]->vektor[0]); n < interval_element_count; n++) {
            v = generate_folgenvector(function->hierarchie[step]->vektor[0], n);
            
            v2 = getSuccessor(v, n, function->hierarchie[step]->vektor[0], direction, interval_element_count);
            position = entry_d2one(v, function->hierarchie[step]->vektor[0]->lang);
            if (v2 != NULL) {
                position2 = entry_d2one(v2, function->hierarchie[step]->vektor[0]->lang);
                result->hierarchie[step]->vektor[0]->glied[position] = (function->hierarchie[step]->vektor[0]->glied[position2] - function->hierarchie[step]->vektor[0]->glied[position])/get_h_l(0, stepping);
                vec_del(v2);
                vec_del(v);
                //printf("%f\n", result->hierarchie[step]->vektor[0]->glied[position]);
            } else {
                
                v2 = getPredeccessor(v, n, function->hierarchie[step]->vektor[0], direction);
                vec_del(v);
                if (v2 != NULL) {
                    v = getPredeccessor(v2, n, function->hierarchie[step]->vektor[0], direction);
                    position2 = entry_d2one(v2, function->hierarchie[step]->vektor[0]->lang);
                    position3 = entry_d2one(v, function->hierarchie[step]->vektor[0]->lang);
                    result->hierarchie[step]->vektor[0]->glied[position] = 2*result->hierarchie[step]->vektor[0]->glied[position2] - result->hierarchie[step]->vektor[0]->glied[position3];
                    vec_del(v2);
                    vec_del(v);
                } else { // this happens only, if stepping == real_interval_length[direction]
                    result->hierarchie[step]->vektor[0]->glied[position] = function->hierarchie[step]->vektor[0]->glied[position];
                }
                //vec_del(v2);
                    
            }
            
            
            
        }
        
    }
    
    return result;
    
}

/**
 * Computes the laplacian of a given function.
 *
 * @param function function that will be 'laplaced'
 * @param stepping stepping
 *
 * @return laplacian of the given function
 */
func_p func_laplace(func_p function, fepc_real_t stepping) {
    func_p result, temp1, temp2;
    
    int n;
    
    for (n = 0; n < function->dim; n++) {
        temp1 = func_derive(function, n, stepping);
        temp2 = func_derive(temp1, n, stepping);
        func_del(temp1);
        if (n == 0) {
            result = temp2;
        } else {
            func_add2(result, temp2);
            func_del(temp2);
        }
    }
    return result;
}

/**
 * This method simplifies the fecp convolution algorithm in a way that one does not need to take care of the structural expections the fepc algorithm has.
 *
 * @param function1       first R^n -> R function of the convolution
 * @param intervals1      refined grid of the first function where the i-th grid has to be finer than the (i-1)-th grid
 * @param interval_count1 number of the intervals 
 * @param degree_func1    degree of every component variable of the basis function that are associated to the first function
 * @param function2       second R^n -> R function of the convolution
 * @param intervals2      refined grid of the second function where the i-th grid has to be finer than the (i-1)-th grid
 * @param interval_count2 number of the intervals 
 * @param degree_func2    degree of every component variable of the basis function that are associated to the second function
 * @param intervals3      refined grid of the convolution result where the i-th grid has to be finer than the (i-1)-th grid
 * @param interval_count3 number of the intervals 
 * @param stepping        stepping of the first level refinement
 *
 * @return function1 * function2 in S(M)
 */
func_p fepc_easy(Funcimpl function1, interval_p* intervals1, int interval_count1, int degree_func1, Funcimpl function2, interval_p* intervals2, int interval_count2, int degree_func2, interval_p* intervals3, int interval_count3, fepc_real_t stepping) {   
    // check for coherence of the parameters
    ASSERT(interval_count1 > 0 && interval_count2 > 0 && interval_count3 > 0 && degree_func1 > -1 && degree_func2 > -1); 
    ASSERT(function1 != NULL && intervals1 != NULL && function2 != NULL && intervals2 != NULL && intervals3 != NULL);

    func_p f, g, w, fepc;
       
    f = create_fepc_structure(function1, intervals1, interval_count1, degree_func1, stepping);
    g = create_fepc_structure(function2, intervals2, interval_count2, degree_func2, stepping);
    w = create_fepc_structure(NULL, intervals3, interval_count3, 0, stepping);
    
    fepc = faltung_fepc(f, g, w, stepping);
    func_del(f);
    func_del(g);
    func_del(w);
    return fepc;
}


