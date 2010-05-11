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
 
#ifndef __FEPCEASY
#define __FEPCEASY

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "faltung.h"
#include "interval.h"

#ifndef INT_STEPS
#define INT_STEPS 20
#endif

#define SQRT_3 1.7320508075688772935274463415058723669428
#define DIV_1_6 0.16666666666666666666666666666666666666666666

typedef fepc_real_t (*Funcimpl_step) (vec_p, vec_p, int, fepc_real_t);

typedef vec_real_p (*Funcimpl_vec_step) (vec_p, vec_p, int, fepc_real_t);

typedef fepc_real_t (*Funcimpl) (vec_real_p);

typedef struct {
    int count;
    vec_real_p* elements;
    fepc_real_t* element_results;
} vec_real_set_t; 

typedef vec_real_set_t * vec_real_set_p;

/**
 * Converts the given intervals from Omega to A form.
 */
interval_p * convert_interval(int steps, interval_p * intervals, fepc_real_t stepping);

/**
 * Sets the gridstructure in the needed form out of a common list of intervals.
 */
void set_gridstructure(func_p function, interval_p* intervals, fepc_real_t stepping);

/**
 * Creates a new set with allocated space for elements.
 */
vec_real_set_p vec_real_set_new(int count);

/**
 * Prints out all elements of a set.
 */
void print_vec_real_set(vec_real_set_p set);

/**
 * Returns the value of the function at the given point.
 */
vec_real_set_p get_value(func_p function, Funcimpl_vec_step generate_points, fepc_real_t stepping);

void vec_real_set_del(vec_real_set_p set);

/**
 * Returns the h_l value corrensponding to the step.
 */
fepc_real_t get_h_l(int step, fepc_real_t stepping);

/**
 * Generates the vector of the folge at the given element_position.
 */
vec_p generate_folgenvector(folge_p folge, int element_position);

/**
 * Returns the number of elements of the given folge.
 */
int get_interval_element_count(folge_p folge);

/**
 * Integrates via trapezoidal rule.
 */
fepc_real_t integrate_coeff_st(Funcimpl function_impl, vec_p v, vec_p p, int step, fepc_real_t stepping);

/**
 * Integrates via Gauss-Legendre.
 */
fepc_real_t integrate_coeff(Funcimpl function_impl, vec_p v, vec_p p, int step, fepc_real_t stepping);

/**
 * Generates the x for integration based on the calc_position and the stepcount.
 */
vec_real_p generate_x_for_integration_st(vec_p v, int calc_position, fepc_real_t h, fepc_real_t h_l, int stepcount);

/**
 * Returns true if the vector is in the latter interval of the folge.
 */
bool_t is_in_latter_interval(vec_p v, folge_p folge);

bool_t is_in_latter_interval_real(vec_real_p v, folge_p folge, int step, fepc_real_t stepping);

/**
 * Returns the value for phi_l (basis function) for the given vector, point x and step.
 */
fepc_real_t phi_l(int step, vec_p v, vec_p p, vec_real_p x, fepc_real_t stepping);

/**
 * Sets the entries of the folge. function_impl xor coeff_function may be NULL. coeff_function is preffered if both are given.
 */
void add_folgenentries(func_p function, Funcimpl function_impl, Funcimpl_step coeff_function, fepc_real_t stepping);

func_p setup_fepc_structure(Funcimpl function, interval_p* intervals, int interval_count, int degree, fepc_real_t stepping);

vec_real_set_p get_value(func_p function, Funcimpl_vec_step generate_points, fepc_real_t stepping);

fepc_real_t get_value_at_step(func_p function, vec_real_p x, int step, fepc_real_t stepping);

/**
 * Finds the best interval.
 */
fepc_real_t get_value_at(func_p function, vec_real_p x, fepc_real_t stepping);

/**
 * Performs a convolution of the given functions in the given intervals with the given degree basis-functions.
 */
func_p fepc_easy(Funcimpl function1, interval_p* intervals1, int interval_count1, int degree_func1, Funcimpl function2, interval_p* intervals2, int interval_count2, int degree_func2, interval_p* intervals3, int interval_count3, fepc_real_t stepping);


fepc_real_t
func_integrate(func_p function, fepc_real_t stepping);

func_p 
func_derive(func_p function, int direction, fepc_real_t stepping);

func_p
func_laplace(func_p function, fepc_real_t stepping);

fepc_real_t max(fepc_real_t a, fepc_real_t b);

int get_degree_count(vec_p degree);

void set_degree(func_p function, int degree);

#endif

