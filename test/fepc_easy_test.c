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




#define __DEBUG

#define A 0.15
#define B 4.
#define C -3.5

#include "fepc_easy.h"
#include "fepc_easy_helper.h"
#include "seconds.h"
#include <math.h>

#define STEPPING 0.5


fepc_real_t const_1(vec_real_p x) {
    return 1.;
}

fepc_real_t f_coeff_peter(vec_p v, vec_p grad, int step) {
    fepc_real_t h_l = get_h_l(step, STEPPING);
    
    return log((v->array[0] +1)*h_l + A) - log(v->array[0]*h_l + A);
}

fepc_real_t g_coeff_peter(vec_p v, vec_p grad, int step) {
    fepc_real_t h_l = get_h_l(step, STEPPING);
    
    return B*h_l*h_l*(v->array[0] + .5) + C *h_l;
}

fepc_real_t f_peter(vec_real_p vector) {
    return 1./(vector->array[0]+A);
}

fepc_real_t g_peter(vec_real_p vector) {
    return B*vector->array[0] + C;
}

void set_result_peter_real(vec_real_set_p set) {
    int n;
    
    for (n = 0; n < set->count; n++) {
        set->element_results[n] -= set->elements[n]->array[1]*((B*(set->elements[n]->array[0]+A)+C)*(log(set->elements[n]->array[0]+A)-log(A))-B*set->elements[n]->array[0]);
    }
}

void test_example_peter() {
    func_p fepc;
    
    int dimension, steps_f, steps_g, steps_w, grad_f, grad_g;
    
    dimension = 2;
    
    steps_f = 4;
    steps_g = 1;
    steps_w = 3;
    
    grad_f = 0;
    grad_g = 0;
    
    interval_p * intervals_f = (interval_p*) malloc(sizeof(interval_p)*steps_f);
    
    
    intervals_f[0] = interval_new(dimension);
    
    intervals_f[0]->start[0] = 0;
    intervals_f[0]->start[1] = 0;
    intervals_f[0]->end[0] = 1;
    intervals_f[0]->end[1] = 1;
    
    
    intervals_f[1] = interval_new(dimension);
    
    intervals_f[1]->start[0] = 0.;
    intervals_f[1]->start[1] = 0.;
    intervals_f[1]->end[0] = 0.2;
    intervals_f[1]->end[1] = 1.;
    
    intervals_f[2] = interval_new(dimension);
    
    intervals_f[2]->start[0] = 0.;
    intervals_f[2]->start[1] = 0.;
    intervals_f[2]->end[0] = 0.2;
    intervals_f[2]->end[1] = 1.0;
    
    intervals_f[3] = interval_new(dimension);
    
    intervals_f[3]->start[0] = 0.;
    intervals_f[3]->start[1] = 0.;
    intervals_f[3]->end[0] = 0.2;
    intervals_f[3]->end[1] = 1.0;
    
    
    interval_p * intervals_g = (interval_p*) malloc(sizeof(interval_p)*steps_g);
    
    intervals_g[0] = interval_new(dimension);
    
    intervals_g[0]->start[0] = 0;
    intervals_g[0]->start[1] = 0;
    intervals_g[0]->end[0] = 1;
    intervals_g[0]->end[1] = 1;
    
    interval_p * intervals_w = (interval_p*) malloc(sizeof(interval_p)*steps_w);
    
    intervals_w[0] = interval_new(dimension);
    
    intervals_w[0]->start[0] = 0;
    intervals_w[0]->start[1] = 0;
    intervals_w[0]->end[0] = 1;
    intervals_w[0]->end[1] = 1;
    
    intervals_w[1] = interval_new(dimension);
    
    intervals_w[1]->start[0] = 0.2;
    intervals_w[1]->start[1] = 0.4;
    intervals_w[1]->end[0] = 0.8;
    intervals_w[1]->end[1] = 1;
    
    intervals_w[2] = interval_new(dimension);
    
    intervals_w[2]->start[0] = 0.4;
    intervals_w[2]->start[1] = 0.8;
    intervals_w[2]->end[0] = 0.6;
    intervals_w[2]->end[1] = 1;    
    
    fepc = fepc_easy(f_peter, intervals_f, steps_f, grad_f, g_peter, intervals_g, steps_g, grad_g, intervals_w, steps_w, STEPPING); // this will be our result
    
    #ifdef __DEBUG
    
    printf("==========================================\n");
    func_print(fepc, 3);
    #endif
    vec_real_set_p points = get_value(fepc, get_mean_points, STEPPING);
    print_vec_real_set(points);
    set_result_peter_real(points);
    print_vec_real_set(points);
}

void test_example_peter_orig() {
    func_p  f, g, w, fepc, ref;
	fepc_real_t  norm;
	fepc_real_t  a,b,c, mesh, temp1, temp2;
	int  dim, i,j , l, pos;
	vec_p  grad, start, lang, r, s;
	vec_real_p  h;
	char dateiname1[20];



	dim = 2;
	a = 0.15;
	b = 4.;
	c = -3.5;
	mesh = 0.2;


	h = vec_real_new(5);
	for(l=0;l<5;l++) {
		h->array[l] = mesh*pow(2,-l); // 5-stufiges h-array anlegen
	}





/*Die Funktion f erstellen*/
	f = func_new(2,dim);

	// hp-Funktion erstellen
	for(l=0;l<=f->maxlevel;l++) {
		grad = vec_new(dim);
		grad->array[0] = 0;
		grad->array[1] = 0;
		f->hierarchie[l] = folgen_vektor_new( grad );
	}

	// Gitterstruktur Stufe 0
	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 0;
	start->array[1] = 0;
	lang->array[0] = 5;
	lang->array[1] = 5;
	f->hierarchie[0]->vektor[0] = folge_new(start,lang);

	// Gitterstruktur Stufe 1 und 2
	for(l=1;l<=f->maxlevel;l++) {
		start = vec_new(dim);
		lang = vec_new(dim);
		start->array[0] = 0;
		start->array[1] = 0;
		lang->array[0] = 4;
		lang->array[1] = pow(2,l)*5;
		f->hierarchie[l]->vektor[0] = folge_new(start,lang);
	}

	
	l=0;
	for(i=2;i<5;i++) {
		for(j=0;j<5;j++) {
			r = vec_new(dim);
			r->array[0] = i;
			r->array[1] = j;
			pos = entry_d2one(r,f->hierarchie[l]->vektor[0]->lang);
			vec_del(r);
			f->hierarchie[l]->vektor[0]->glied[pos] = log((i+1)*h->array[l]+a) - log(i*h->array[l]+a);
		}
	}


	for(l=1;l<f->maxlevel;l++) {
		for(i=2;i<4;i++) {
			for(j=0;j<(5*pow(2,l));j++) {
				r = vec_new(dim);
				r->array[0] = i;
				r->array[1] = j;
				pos = entry_d2one(r,f->hierarchie[l]->vektor[0]->lang);
				vec_del(r);
				f->hierarchie[l]->vektor[0]->glied[pos] = log((i+1)*h->array[l]+a) - log(i*h->array[l]+a);
			}
		}
	}


	l=f->maxlevel;
	for(i=0;i<4;i++) {
		for(j=0;j<(5*pow(2,l));j++) {
			r = vec_new(dim);
			r->array[0] = i;
			r->array[1] = j;
			pos = entry_d2one(r,f->hierarchie[l]->vektor[0]->lang);
			vec_del(r);
			f->hierarchie[l]->vektor[0]->glied[pos] = log((i+1)*h->array[l]+a) - log(i*h->array[l]+a);
		}
	}




	/*Die Funktion g erstellen*/

	g = func_new(0,dim);

	for(l=0;l<=g->maxlevel;l++) {
		grad = vec_new(dim);
		grad->array[0] = 0;
		grad->array[1] = 0;
		g->hierarchie[l] = folgen_vektor_new( grad );
	}

	for(l=0;l<=g->maxlevel;l++) {
		start = vec_new(dim);
		lang = vec_new(dim);
		start->array[0] = 0;
		lang->array[0] = 5*pow(2,l);
		start->array[1] = 0;
		lang->array[1] = 5*pow(2,l);
		g->hierarchie[l]->vektor[0] = folge_new(start,lang);
	}

	l=g->maxlevel;
	for(i=0;i<(5*pow(2,l));i++) {
		for(j=0;j<(5*pow(2,l));j++) {
			r = vec_new(dim);
			r->array[0] = i;
			r->array[1] = j;
			pos = entry_d2one(r,g->hierarchie[l]->vektor[0]->lang);
			vec_del(r);
			temp1 = b/2.*pow(h->array[l],2);
			temp2 = 2*i + 1;
			g->hierarchie[l]->vektor[0]->glied[pos] = temp1*temp2 + c*h->array[l];
		}
	}


	/*Die Funktion w erstellen*/
	w = func_new(2,dim);

	for(l=0;l<=w->maxlevel;l++) {
		grad = vec_new(dim);
		grad->array[0] = 0;
		grad->array[1] = 0;
		w->hierarchie[l] = folgen_vektor_new( grad );
	}

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 0;
	lang->array[0] = 5;
	start->array[1] = 0;
	lang->array[1] = 5;
	w->hierarchie[0]->vektor[0] = folge_new(start,lang);

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 2;
	lang->array[0] = 6;
	start->array[1] = 4;
	lang->array[1] = 6;
	w->hierarchie[1]->vektor[0] = folge_new(start,lang);

	start = vec_new(dim);
	lang = vec_new(dim);
	start->array[0] = 8;
	lang->array[0] = 4;
	start->array[1] = 16;
	lang->array[1] = 4;
	w->hierarchie[2]->vektor[0] = folge_new(start,lang);
printf("==========================================\n");printf("==========================================\n");
    #ifdef __DEBUG
    //print_func(f);
    //print_func(g);
    //print_func(w);
    #endif


	/*Durchführen der Berechnungen*/
	fepc = faltung_fepc(f,g,w,mesh);
	ref = faltung_ref(f,g,w,mesh);
     #ifdef __DEBUG
    printf("==========================================\n");
    //print_func(fepc);
    #endif
    vec_real_set_p points = get_value(fepc, get_mean_points, STEPPING);
    norm = faltung_hilfe_norm(fepc,ref);
	printf("\n Norm %lf \n",norm);
    print_vec_real_set(points);

}


void some_other_test() {
    func_p f, g, w;
    
    int dimension = 2;
    
    int steps_f = 3;
    
    interval_p * intervals = (interval_p*) malloc(sizeof(interval_p)*steps_f);
    
    intervals[0] = interval_new(dimension);
    
    intervals[0]->start[0] = 0;
    intervals[0]->start[1] = 0;
    intervals[0]->end[0] = 2;
    intervals[0]->end[1] = 2;
    
    intervals[1] = interval_new(dimension);
    
    intervals[1]->start[0] = 0.8;
    intervals[1]->start[1] = 0.8;
    intervals[1]->end[0] = 2;
    intervals[1]->end[1] = 2;
    
    intervals[2] = interval_new(dimension);
    
    intervals[2]->start[0] = 1.2;
    intervals[2]->start[1] = 1.2;
    intervals[2]->end[0] = 2;
    intervals[2]->end[1] = 2;
    
    f = func_new(steps_f-1, dimension);
    
    #ifdef __DEBUG
    print_intervals(intervals, steps_f);
    printf("\nwill be converted to\n");
    print_intervals(convert_interval(steps_f, intervals, STEPPING), steps_f);
    #endif
    //print_intervals(intervals, 3);
    set_gridstructure(f, intervals, STEPPING); 
    
    int n;
    printf("Hier2\n");
    for (n = 0; n < steps_f; n++) {// add the hp-structure of f
        //add_folgenentries(f, n);
    }
    
    
    
    
    
    
    
    g = func_new(0, dimension);
    
    
    w = func_new(2, dimension);
    
    // generate the grid structure for w
    
    //func_p fepc = faltung_fepc(f, g, w, STEPPING); // this will be our result
    

}

void test_integration() {
    int step = 2;
    
    
    
    vec_p v = vec_new(2);
    v->array[0] = 0.;
    v->array[1] = 0.;
    
    vec_p p = vec_new(2);
    
    fepc_real_t result, h_l, check_result, result_2;
    
    h_l = get_h_l(step, STEPPING);
    
    
    check_result = g_coeff_peter(v, p, step);
    result = integrate_coeff(g_peter, v, p, step, STEPPING);
    result_2 = integrate_coeff_st(g_peter, v, p, step, STEPPING);
    double seconds_used = seconds();
    printf("\nIntegration %.10f\n", integrate_coeff(f_peter, v, p, step, STEPPING));
    printf("seconds: %f\n", seconds()-seconds_used);
    seconds_used = seconds();
    printf("\nIntegration_st %.10lf\n", integrate_coeff_st(f_peter, v, p, step, STEPPING));
    printf("seconds: %f\n", seconds()-seconds_used);
    printf("real result \n======\n %.10lf\n", f_coeff_peter(v, p, step));
    
    printf("\nIntegrate %.10f\nExact %.10lf\nIntegrate_st %.10lf\n", result, check_result, result_2);
    vec_del(v);
    vec_del(p);
}

void test_multi() {
    func_p f, g, multi;
    
    int dimension, steps_f, steps_g, grad_f, grad_g;
    
    dimension = 2;
    steps_f = 2;
    steps_g = 2;
    
    grad_f = 1;
    grad_g = 0;
    
    interval_p * intervals = (interval_p*) malloc(sizeof(interval_p)*steps_f);
    
    intervals[0] = interval_new(dimension);
    
    intervals[0]->start[0] = 0;
    intervals[0]->start[1] = 0;
    intervals[0]->end[0] = 1;
    intervals[0]->end[1] = 1;
    
    intervals[0] = interval_new(dimension);
    
    intervals[0]->start[0] = 0;
    intervals[0]->start[1] = 0;
    intervals[0]->end[0] = 1;
    intervals[0]->end[1] = 1;
   
    intervals[1] = interval_new(dimension);
    
    intervals[1]->start[0] = 0.;
    intervals[1]->start[1] = 0.;
    intervals[1]->end[0] = 0.4;
    intervals[1]->end[1] = 1.;

    f = setup_fepc_structure(f_peter, intervals, steps_f, grad_f, STEPPING);
    g = setup_fepc_structure(g_peter, intervals, steps_g, grad_g, STEPPING);
    multi = func_multi(f, g, STEPPING);
    
    vec_real_set_p points = get_value(multi, get_mean_points, STEPPING);
    print_vec_real_set(points);
    points = get_value(f, get_mean_points, STEPPING);
    print_vec_real_set(points);
    points = get_value(g, get_mean_points, STEPPING);
    print_vec_real_set(points);
}

void test_norm() {
    int dimension = 2;
    
    interval_p * intervals = (interval_p*) malloc(sizeof(interval_p)*1);

    
    intervals[0] = interval_new(dimension);
           
    intervals[0]->start[0] = 0.;
    intervals[0]->start[1] = 0.;
    intervals[0]->end[0] = 1.;
    intervals[0]->end[1] = 1.;

    func_p function = setup_fepc_structure(const_1, intervals, 1, 0, STEPPING);
    
    //func_print(function, 3);
    
    fepc_real_t integral = func_integrate(function, STEPPING);
    
    //vec_real_set_p set = get_value(function, get_mean_points, STEPPING);
    
    //print_vec_real_set(set);
    printf("Integral = %.10lf\n", integral);
    
}

int main(void) {
    test_norm();
    test_integration();
    //test_example_peter_orig();
    //test_example_peter();
    //test_multi();
    return 0;
}




