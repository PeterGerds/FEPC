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

#if defined(HAS_FFTW3)
#include <fftw3.h>

#endif

#include "fft_faltung.h"


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
fepc_real_t*
fft_faltung(fepc_real_t* a, vec_p n_a, fepc_real_t* b, vec_p n_b) {
#if defined(HAS_FFTW3)

	int  size_a, size_b, size_c, dim;
	int  k, i, wert, test;
	int  *n;
	vec_p  temp, n_c;
	fepc_real_t  *c;
	fftw_complex  *in, *out,*in_a, *out_a, *in_b, *out_b;
	fftw_plan  p;



	/*Auf Testen von Konsistenz wird verzichtet, da Input bereits auf Konsistenz getestet*/

	/*Faltung ueber Fouriertrafo (Theorie ist in Dokumentation zu finden)*/
	dim = n_a->dim;
	n_c = vec_new(dim);
	for(k=0;k<dim;k++) {
		n_c->array[k] = n_a->array[k]+n_b->array[k]-1;
	}
	n = n_c->array;
	size_a = vec_size( n_a );
	size_b = vec_size( n_b );
	size_c = vec_size( n_c );


	/*Initialisieren des Ergebnis Array*/
	c = (fepc_real_t*) malloc(sizeof(fepc_real_t) * size_c);


	/*Berechnen der Fouriertrafo von in_a*/
	in_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	out_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	for (k=0;k<size_c;k++) {
		temp = entry_one2d(k,n_c);
		test = 0;
		for(i=0;i<dim;i++) {
			if ((temp->array[i] <0)||(temp->array[i]>=n_a->array[i])) {
				test = test + 1;
			}
		}
		if (test == 0) {
			wert = entry_d2one(temp,n_a);
			in_a[k][0] = a[wert];
			in_a[k][1] = 0;
		}
		else {
			in_a[k][0] = 0;
			in_a[k][1] = 0;
		}
		vec_del(temp);
	}
	p = fftw_plan_dft(dim,n,in_a,out_a,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);


	/*Berechnen der Fouriertrafo von in_b*/
	in_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	out_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	for (k=0;k<size_c;k++) {
		temp = entry_one2d(k,n_c);
		test = 0;
		for(i=0;i<dim;i++) {
			if ((temp->array[i] <0)||(temp->array[i]>=n_b->array[i])) {
				test = test + 1;
			}
		}
		if (test == 0) {
			wert = entry_d2one(temp,n_b);
			in_b[k][0] = b[wert];
			in_b[k][1] = 0;
		}
		else {
			in_b[k][0] = 0;
			in_b[k][1] = 0;
		}
		vec_del(temp);
	}

	p = fftw_plan_dft(dim,n,in_b,out_b,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);


	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);

	for (k=0;k<size_c;k++) {
		in[k][0] = out_a[k][0]*out_b[k][0] - out_a[k][1]*out_b[k][1];
		in[k][1] = out_a[k][1]*out_b[k][0] + out_a[k][0]*out_b[k][1];
	}

	/*Berechnung der Inversen Fouriertrafo von in*/
	p = fftw_plan_dft(dim,n,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	for (k=0;k<size_c;k++) {
		c[k] = (fepc_real_t) out[k][0]/size_c;
	}

	vec_del(n_c);
	fftw_free(in);
	fftw_free(in_a);
	fftw_free(in_b);
	fftw_free(out);
	fftw_free(out_a);
	fftw_free(out_b);
	return c;



#else
    printf( "\n (fft_faltung) FEHLER : keine FFT Bibliothek verfuegbar\n" );

    exit( 1 );
#endif
}
