/*
 * FEPC
 * Copyright (C) 2009 Peter Gerds (gerds@mis.mpg.de)
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

#include "funktion.h"


/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/





/*******************************************************
 *
 * globale Funktionen
 *
 *******************************************************/


func_p
func_new(int maxlevel, int dim) {
	func_p  back;
	folgen_vektor_p  *hierarchie;
	int  k;
	vec_p  grad;

	ASSERT(maxlevel >= 0);
	ASSERT(dim > 0);
	back = (func_p) malloc(sizeof(func_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_vektor_p*) malloc(sizeof(folgen_vektor_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);


	for(k=0;k<=maxlevel;k++) {
		grad = vec_new(dim);
		hierarchie[k] = folgen_vektor_new(grad);
	}
	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;
	back->dim = dim;
	return back;
}



void
func_del(func_p f) {
    if (f && f != NULL) {
        int  k;

	    for(k=0;k<=f->maxlevel;k++) {
		    folgen_vektor_del(f->hierarchie[k]);
	    }
	    free(f->hierarchie);
	    free(f);
    }
}






func2_p
func2_new(int maxlevel, int dim) {
	func2_p  back;
	folgen_matrix_p  *hierarchie;
	int  k;
	vec_p  grad1, grad2;

	ASSERT(maxlevel >= 0);
	ASSERT(dim > 0);
	back = (func2_p) malloc(sizeof(func2_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_matrix_p*) malloc(sizeof(folgen_matrix_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		grad1 = vec_new(dim);
		grad2 = vec_new(dim);
		hierarchie[k] = folgen_matrix_new(grad1,grad2);
	}
	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;
	back->dim = dim;

	return back;
}


void
func2_del(func2_p f) {
	int  k;

	for(k=0;k<=f->maxlevel;k++) {
		folgen_matrix_del(f->hierarchie[k]);
	}
	free(f->hierarchie);
	free(f);
}



func_p
func_projekt(func_p f,func_p g) {
	func_p  back;
	int  l;

	back = func_new( g->maxlevel, g->dim );
	for(l=0;l<=g->maxlevel;l++) {
		if (l<=f->maxlevel) {
			folgen_vektor_del( back->hierarchie[l]);
			back->hierarchie[l] = folgen_vektor_projekt( f->hierarchie[l], g->hierarchie[l] );
		}
		else {
			folgen_vektor_del( back->hierarchie[l] );
			back->hierarchie[l] = folgen_vektor_copy_structure( g->hierarchie[l] );
		}
	}
	return back;
}

func_p
func_clone(func_p f) {
    return func_projekt(f, f);
}



void
func_print(func_p f, int info) {
	int  k;

	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	printf("\n#Ausgabe einer Funktion der Dimension %d und des Maxlevel %d ", f->dim, f->maxlevel);
	if( (info == 1) || (info == 2) || (info == 3) ) {
		info = info -1;
		for(k=0;k<=f->maxlevel;k++) {
			folgen_vektor_print( f->hierarchie[k], info );
		}
	}
	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n\n");
	return;
}


func_p
func_build( int maxlevel, int dim , int grad, int a, int n, int mod, bool_t random) {
	func_p  back;
	int  k, max;
	vec_p  p;

	ASSERT( maxlevel >= 0 );
	ASSERT( dim >= 0 );
	ASSERT( grad >= 0 );
	if( random == true) {
		max = rand()%(maxlevel+1);
	}
	else {
		max = maxlevel;
	}
	back = func_new( max, dim );
	p = vec_new(dim);
	for(k=0;k<dim;k++) {
		p->array[k] = grad;
	}
	for(k=0;k<=max;k++) {
		folgen_vektor_del( back->hierarchie[k] );
		back->hierarchie[k] = folgen_vektor_build( p, a, n, mod, random);
	}
	vec_del( p );
	return back;
}


func_p
func_add(func_p f, func_p g) {
	func_p  back;
	int  k, maxlevel, dim;
	folgen_vektor_p  temp;

	ASSERT( f->dim == g->dim );
	dim = f->dim;
	if (f->maxlevel < g->maxlevel) {
		maxlevel = g->maxlevel;
	}
	else {
		maxlevel = f->maxlevel;
	}


	back = func_new( maxlevel, dim );
	for (k=0;k<=maxlevel;k++) {
		if (k <= f->maxlevel) {
			temp = folgen_vektor_add( back->hierarchie[k], f->hierarchie[k] );
			folgen_vektor_del( back->hierarchie[k] );
			back->hierarchie[k] = temp;
		}
		if (k <= g->maxlevel) {
			temp = folgen_vektor_add( back->hierarchie[k], g->hierarchie[k] );
			folgen_vektor_del( back->hierarchie[k] );
			back->hierarchie[k] = temp;
		}
	}

	return back;
}


int
func_count( func_p f, func_p g, func_p w ) {
	int  back;
	int  l, k, size;
	vec_p  vec_1, n, lang;

	vec_1 = vec_one( f->dim );

	back = 0;

	for(l=0;l<=f->maxlevel;l++) {
		n = vec_add( f->hierarchie[l]->grad, vec_1 );
		size = vec_size( n );
		for(k=0;k<size;k++) {
			lang =  f->hierarchie[l]->vektor[k]->lang;
			back = back + vec_size( lang );
		}
		vec_del( n );
	}

	for(l=0;l<=g->maxlevel;l++) {
		n = vec_add( g->hierarchie[l]->grad, vec_1 );
		size = vec_size( n );
		for(k=0;k<size;k++) {
			lang =  g->hierarchie[l]->vektor[k]->lang;
			back = back + vec_size( lang );
		}
		vec_del( n );
	}

	for(l=0;l<=w->maxlevel;l++) {
		n = vec_add( w->hierarchie[l]->grad, vec_1 );
		size = vec_size( n );
		for(k=0;k<size;k++) {
			lang =  w->hierarchie[l]->vektor[k]->lang;
			back = back + vec_size( lang );
		}
		vec_del( n );
	}


	vec_del( vec_1 );
	return back;
}


int
func_modell_count( func_p f, func_p g, func_p w ){
	int  back;
	int  l, min_f_g, min_f_w, min_g_w;
	vec_p  temp;


	if (f->maxlevel < g->maxlevel) {
		min_f_g = f->maxlevel;
	}
	else {
		min_f_g = g->maxlevel;
	}

	if (f->maxlevel < w->maxlevel) {
		min_f_w = f->maxlevel;
	}
	else {
		min_f_w = w->maxlevel;
	}

	if (g->maxlevel < w->maxlevel) {
		min_g_w = g->maxlevel;
	}
	else {
		min_g_w = w->maxlevel;
	}


	back = 0;

	for(l=0;l<=min_f_g;l++) {
		temp = vec_add( f->hierarchie[l]->vektor[0]->lang, g->hierarchie[l]->vektor[0]->lang );
		back = back + vec_size( temp );
		vec_del( temp );
	}

	for(l=0;l<=min_f_w;l++) {
		temp = vec_add( f->hierarchie[l]->vektor[0]->lang, w->hierarchie[l]->vektor[0]->lang );
		back = back + vec_size( temp );
		vec_del( temp );
	}

	for(l=0;l<=min_g_w;l++) {
		temp = vec_add( g->hierarchie[l]->vektor[0]->lang, w->hierarchie[l]->vektor[0]->lang );
		back = back + vec_size( temp );
		vec_del( temp );
	}

	for(l=min_g_w+1;l<=f->maxlevel;l++) {
		temp = f->hierarchie[l]->vektor[0]->lang;
		back = back + vec_size( temp );
	}

	for(l=min_f_w+1;l<=g->maxlevel;l++) {
		temp = g->hierarchie[l]->vektor[0]->lang;
		back = back + vec_size( temp );
	}

	for(l=min_f_g+1;l<=w->maxlevel;l++) {
		temp = w->hierarchie[l]->vektor[0]->lang;
		back = back + vec_size( temp );
	}

	return back;
}


void
func_grid_zero(func_p f) {
	func_p  back;
	int  l, i, d, k;
	vec_p  n, grad, start, lang, vec_1, r, s;
	int  test, size, size_grad, a, b;

	vec_1 = vec_one(f->dim);

	for(l=0;l<f->maxlevel;l++) {
		grad =  f->hierarchie[l]->grad;
		n = vec_add(grad, vec_1);
		size_grad = vec_size(n);
		vec_del(n);

		lang = f->hierarchie[l]->vektor[0]->lang;
		start = f->hierarchie[l]->vektor[0]->start;
		size = vec_size(lang);

		for(i=0;i<size;i++) {
			/*Teste ob Folgenglied durch Folgenglieder hoeherer Level repraesentiert wird*/
			r = entry_one2d(i,lang);
			s = vec_add(r,start);
			vec_del(r);
			r = vec_multi(2, s);
			vec_del(s);
			test = 0;
			for(d=0;d<f->dim;d++) {
				a = f->hierarchie[l+1]->vektor[0]->start->array[d];
				b = f->hierarchie[l+1]->vektor[0]->lang->array[d] + a - 1;
				if( (r->array[d]<a) || (r->array[d]>b) ) {
					test = test + 1;
				}
			}

			vec_del(r);
			if(test==0) {
				for(k=0;k<size_grad;k++) {
					f->hierarchie[l]->vektor[k]->glied[i] = 0;
				}
			}
		}
	}
	vec_del(vec_1);
}

func_p 
func_factor_multi(func_p function, fepc_real_t factor) {
    func_p  back = func_new(function->maxlevel, function->dim);
    
    folgen_vektor_p  temp;
    
    int k;
    
	for (k=0;k<=function->maxlevel;k++) {
		temp = folgen_vektor_factor_multi( function->hierarchie[k], factor );
		folgen_vektor_del( back->hierarchie[k] );
		back->hierarchie[k] = temp;
	}

	return back;
}

void
funcs_del(func_p * array, int length) {
	int n;

	for (n = 0; n < length; n++) {
		func_del(array[n]);
	}
	free(array);
}
