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

#include "folge.h"


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

folge_p
folge_new(vec_p start,vec_p lang) {
	folge_p back;
	fepc_real_t *glied;
	int k,size;

	/*Testen auf Konsistenz*/
	ASSERT(start->dim == lang->dim);
	for(k=0;k<lang->dim;k++) {
		ASSERT( lang->array[k] >= 0 );
	}
	back = (folge_p) malloc(sizeof(folge_t));
	ASSERT(back != NULL);

	size = vec_size( lang );
	glied = (fepc_real_t*) malloc(sizeof(fepc_real_t)*size);
	ASSERT(glied != NULL);

	for (k=0;k<size;k++) {
		glied[k] = 0.0;
	}

	back->glied = glied;
	back->lang = lang;
	back->start = start;

	return back;
}


void
folge_del(folge_p f) {
	free(f->glied);
	vec_del(f->start);
	vec_del(f->lang);
	free(f);
}



folge_p
folge_slow_faltung(folge_p f,folge_p g) {
	folge_p back;
	int  k, l, i, test, wert;
	int  size_f, size_g, size_w, dim;
	vec_p  n_f, n_g, n_w, w_start;
	vec_p  r, s, temp;
	fepc_real_t  *f_glied, *g_glied, *w_glied;
	fepc_real_t  sum;

	/*Testen auf Konsistenz*/
	ASSERT(f->start->dim == f->lang->dim);
	ASSERT(g->start->dim == g->lang->dim);
	ASSERT(g->start->dim == f->start->dim);

	dim = g->start->dim;

	n_f = f->lang;
	size_f = vec_size( n_f );

	n_g = g->lang;
	size_g = vec_size( n_g );

	for(k=0;k<dim;k++) {
		ASSERT( f->lang->array[k] >= 0 );
		ASSERT( g->lang->array[k] >= 0 );
	}


	/*Rueckgabe einer leeren Folge, falls f oder g leere Folgen sind*/
	if ( (size_f == 0) || (size_g == 0) ) {
		r = vec_new(dim);
		s = vec_new(dim);
		back = folge_new( r, s );
		return back;
	}

	/*Berechnungen*/
	n_w = vec_new(dim);
	for(k=0;k<dim;k++) {
		n_w->array[k] = n_f->array[k]+n_g->array[k]-1;
	}


	f_glied = f->glied;
	g_glied = g->glied;

	/*Initialisierung der Faltung*/
	w_start = vec_add( f->start, g->start );
	back = folge_new(w_start , n_w);
	w_glied = back->glied;
	size_w = vec_size( n_w );

	/*Berechnung der Faltung*/
	for(k=0;k<size_w;k++) {
		sum = 0;
		r = entry_one2d(k, n_w);
		for(l=0;l<size_f;l++) {
			s = entry_one2d(l,n_f);
			temp = vec_op( 1, r, -1, s);
			test = 0;
			for(i = 0;i<dim;i++) {
				if ( (temp->array[i]<0) || (temp->array[i]>=n_g->array[i]) ) {
					test = test + 1;
				}
			}
			if (test == 0) {
				wert = entry_d2one( temp, n_g );
				sum = sum + f_glied[l] * g_glied[wert];
			}
			vec_del(s);
			vec_del(temp);
		}
	w_glied[k] = sum;
	vec_del(r);
	}

	return back;
}




folge_p
folge_faltung(folge_p f,folge_p g) {
	folge_p  back;
	int  k;
	int  size_f, size_g, dim;
	vec_p  n_f, n_g, n_w, w_start;
	vec_p  r, s;
	fepc_real_t  *f_glied, *g_glied, *w_glied;


	/*Testen auf Konsistenz*/
	ASSERT(f->start->dim == f->lang->dim);
	ASSERT(g->start->dim == g->lang->dim);
	ASSERT(g->start->dim == f->start->dim);

	dim = g->start->dim;

	n_f = f->lang;
	size_f = vec_size( n_f );

	n_g = g->lang;
	size_g = vec_size( n_g );


	for(k=0;k<dim;k++) {
		ASSERT( f->lang->array[k] >= 0 );
		ASSERT( g->lang->array[k] >= 0 );
	}



	/*Rueckgabe einer leeren Folge, falls f oder g leere Folgen sind*/
	if ( (size_f == 0) || (size_g == 0) ) {
		r = vec_new(dim);
		s = vec_new(dim);
		back = folge_new( r, s );
		return back;
	}

	/*Berechnungen*/
	n_w = vec_new(dim);
	for(k=0;k<dim;k++) {
		n_w->array[k] = n_f->array[k]+n_g->array[k]-1;
	}


	f_glied = f->glied;
	g_glied = g->glied;


	/*Berechnung der Faltung (Theorie hierzu ist in der Dokumentation zu finden)*/
	w_glied = fft_faltung( f_glied, n_f, g_glied, n_g);


	/*Initialisierung des Ergebnisses*/
	w_start = vec_add( f->start, g->start );
	back = folge_new(w_start , n_w);
	free( back->glied );
	back->glied = w_glied;

	return back;
}



void
folge_print(folge_p f, int info) {
	int  k, dim, size;

	dim = f->lang->dim;
	printf("\n------------------------------------------------------------");
	printf("\n#Ausgabe einer Folge der Dimension %d", dim);
	printf("\n\t#:start");
	for(k=0;k<dim;k++){
		printf("\t %d",f->start->array[k]);
	}
	printf("\n\t#:lang");
	for(k=0;k<dim;k++){
		printf("\t %d",f->lang->array[k]);
	}

	size = vec_size( f->lang );
	if(info == 1) {
		printf("\n\t#:glied");
		for(k=0;k<size;k++) {
			printf("\t%.5lf",f->glied[k]);
		}
	}
	printf("\n------------------------------------------------------------\n");

	return;
}


folge_p
folge_build(int dim, int a, int n, int mod, bool_t random) {
	vec_p  start, lang;
	int  k, size;
	fepc_real_t *glied;
	folge_p back;

	ASSERT(a>=0);
	ASSERT(n>0);
	start = vec_new(dim);
	lang = vec_new(dim);

	for(k=0;k<dim;k++) {
		if( random == true ) {
			start->array[k] = rand()%(2*a +1 ) - a;
			lang->array[k] = rand()%(n) +1;
		}
		else {
			start->array[k] = a;
			lang->array[k] = n;
		}
	}
	back = folge_new(start,lang);

	size = vec_size( lang );
	glied = back->glied;
	for (k=0;k<size;k++) {
		glied[k] = rand()%(2*mod+1) - mod;
	}

	return back;
}


fepc_real_t
folge_glied( vec_p r, folge_p f ) {
	int  k, dim;
	vec_p  anfang, ende, temp, vec_1;
	fepc_real_t  back;

	dim = f->start->dim;
	vec_1 = vec_one(dim);

	anfang = f->start;
	temp = vec_add( anfang, f->lang );
	ende = vec_op( 1, temp, -1, vec_1 );
	vec_del(temp);
	vec_del(vec_1);

	/*Falls der Vektor r nicht im Traeger von f enthalten ist
	Rueckgabe von 0.0*/
	for(k=0;k<dim;k++) {
		if (( r->array[k] > ende->array[k] ) || ( r->array[k] < anfang->array[k] )) {
			back = 0.0;
			vec_del(ende);
			return back;
		}
	}
	vec_del(ende);

	/*Vektor r ist im Traeger von f enthalten, ermitteln des entsprechenden Folgenwertes*/
	temp = vec_op( 1, r, -1, f->start);
	k = entry_d2one( temp, f->lang);
	vec_del(temp);

	back = f->glied[k];
	return back;
}




fepc_real_t
folge_norm(folge_p f, folge_p g) {
	vec_p  temp, temp1, temp2, min, max, lang, vec_1;
	int  k, dim, size;
	fepc_real_t  norm, diff;

	ASSERT(f->start->dim == g->start->dim);
	dim = f->start->dim;
	vec_1 = vec_one( dim );

	min = vec_min( f->start, g->start );
	temp1 = vec_add( f->start, f->lang );
	temp2 = vec_add( g->start, g->lang );
	temp = vec_max( temp1, temp2 );
	vec_del( temp1 );
	vec_del( temp2 );
	lang = vec_op( 1, temp, -1, min );
	vec_del( temp );

	size = vec_size( lang );
	norm = 0.0;
	for(k=0;k<size;k++) {
		temp = entry_one2d( k, lang );
		temp1 = vec_add( temp, min );
		diff = folge_glied( temp1, f ) - folge_glied( temp1, g );
		norm = norm + ( diff * diff );
		vec_del( temp );
		vec_del( temp1 );
	}
	vec_del( vec_1 );
	vec_del( min );
	vec_del( lang );

	norm = sqrt(norm);
	return norm;
}




folge_p
folge_add(folge_p f, folge_p g) {
	folge_p  back;
	int  size, k, size_g, size_f;
	fepc_real_t  x, y;
	vec_p  temp1, temp2, max, lang, min, r;


	size_f = vec_size( f->lang );
	size_g = vec_size( g->lang );
	if(size_f == 0) {
		back = folge_copy( g );
		return back;
	}

	if(size_g == 0) {
		back = folge_copy( f );
		return back;
	}

	if( (size_g!=0) && (size_f!=0) ) {
		min = vec_min( f->start, g->start );
		temp1 = vec_add( f->start, f->lang );
		temp2 = vec_add( g->start, g->lang );
		max = vec_max( temp1, temp2 );
		vec_del( temp1 );
		vec_del( temp2 );
		lang = vec_op( 1, max, -1, min);
		vec_del( max );
		back = folge_new( min, lang );
		size = vec_size( lang );
		for(k=0;k<size;k++) {
			temp1 = entry_one2d( k, lang );
			r = vec_add( min, temp1 );
			vec_del( temp1 );
			x = folge_glied( r, f );
			y = folge_glied( r, g );
			vec_del( r );
			back->glied[k] = x + y;
		}
		return back;
	}
}

folge_p
folge_subtract(folge_p f, folge_p g) {
	folge_p  back;
	int  size, k, size_g, size_f;
	fepc_real_t  x, y;
	vec_p  temp1, temp2, max, lang, min, r;


	size_f = vec_size( f->lang );
	size_g = vec_size( g->lang );
	if(size_f == 0) {
		back = folge_copy( g );
		folge_multi_factor(back, -1);
		return back;
	}

	if(size_g == 0) {
		back = folge_copy( f );
		return back;
	}

	if( (size_g!=0) && (size_f!=0) ) {
		min = vec_min( f->start, g->start );
		temp1 = vec_add( f->start, f->lang );
		temp2 = vec_add( g->start, g->lang );
		max = vec_max( temp1, temp2 );
		vec_del( temp1 );
		vec_del( temp2 );
		lang = vec_op( 1, max, -1, min);
		vec_del( max );
		back = folge_new( min, lang );
		size = vec_size( lang );
		for(k=0;k<size;k++) {
			temp1 = entry_one2d( k, lang );
			r = vec_add( min, temp1 );
			vec_del( temp1 );
			x = folge_glied( r, f );
			y = folge_glied( r, g );
			vec_del( r );
			back->glied[k] = x - y;
		}
		return back;
	}
}


folge_p
folge_copy( folge_p f) {
	folge_p  back;
	int  k, size;
	vec_p  start, lang;

	start = vec_copy( f->start );
	lang = vec_copy( f->lang );
	size = vec_size( lang );
	back = folge_new( start, lang );
	for(k=0;k<size;k++) {
		back->glied[k] = f->glied[k];
	}

	return back;
}



folge_p
folge_projekt(folge_p f, folge_p g) {
	folge_p  back;
	int  k, size;
	vec_p  r, start, lang, temp;

	ASSERT( f->start->dim == g->start->dim );
	start = vec_copy( g->start );
	lang = vec_copy( g->lang );

	back = folge_new( start, lang );
	size = vec_size( lang );
	for(k=0;k<size;k++) {
		temp = entry_one2d( k, lang );
		r = vec_add( temp, start );
		vec_del( temp );
		back->glied[k] = folge_glied( r, f );
		vec_del( r );
	}

	return back;
}

folge_p
folge_multi_factor(folge_p folge, fepc_real_t factor) {
    folge_p  back;
    vec_p lang, start;
	int  k, size;

    start = vec_copy( folge->start );
	lang = vec_copy( folge->lang );

	back = folge_new( start, lang );
	size = vec_size( lang );
	
	for(k=0;k<size;k++) {

		back->glied[k] = folge->glied[k]*factor;
	}

	return back;
}




