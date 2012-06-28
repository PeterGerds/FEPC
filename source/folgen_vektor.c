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

#include "folgen_vektor.h"


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


folgen_vektor_p
folgen_vektor_new(vec_p grad) {
	folgen_vektor_p  back;
	folge_p  *vektor;
	int  k, size, dim;
	vec_p  r, s;

	/*Testen auf Konsistenz*/
	dim = grad->dim;
	size = 1;
	for(k=0;k<dim;k++) {
		ASSERT( grad->array[k] >= 0 );
		size = size * (grad->array[k] + 1);
	}



	back = (folgen_vektor_p) malloc(sizeof(folgen_vektor_t));
	ASSERT(back != NULL);


	vektor = (folge_p*) malloc(sizeof(folge_p) * size);
	ASSERT(vektor != NULL);

	for(k=0;k<size;k++) {
		r = vec_new(dim);
		s = vec_new(dim);
		vektor[k] = folge_new( r, s);
	}
	back->vektor = vektor;
	back->grad = grad;

	return back;
}



void
folgen_vektor_del(folgen_vektor_p f) {
	int  k, size, dim;
	int  *array;

	dim = f->grad->dim;
	array = f->grad->array;
	size = 1;
	for(k=0;k<dim;k++) {
		size = size * (array[k] + 1);
	}

	for(k=0;k<size;k++) {
		folge_del(f->vektor[k]);
	}
	free(f->vektor);
	vec_del(f->grad);
	free(f);
}



folgen_matrix_p
folgen_matrix_new(vec_p grad1, vec_p grad2) {
	folgen_matrix_p  back;
	folge_p  **matrix;
	int  i, j, size1, size2, dim;
	vec_p  r, s;

	/*Testen auf Konsistenz*/
	ASSERT(grad1->dim == grad2->dim);
	dim = grad1->dim;
	size1 = 1;
	for(i=0;i<dim;i++) {
		ASSERT( grad1->array[i] >= 0 );
		size1 = size1 * (grad1->array[i] + 1);
	}
	size2 = 1;
	for(i=0;i<dim;i++) {
		ASSERT( grad2->array[i] >= 0 );
		size2 = size2 * (grad2->array[i] + 1);
	}


	back = (folgen_matrix_p) malloc(sizeof(folgen_matrix_t));
	ASSERT(back != NULL);
	matrix = (folge_p**) malloc(sizeof(folge_p*) * size1);
	ASSERT(matrix != NULL);

	matrix[0] = (folge_p*) malloc(sizeof(folge_p) * size1 * size2);
	for(i=0;i<size1;i++) {
		matrix[i] = matrix[0] + i * size2;
	}

	for(i=0;i<size1;i++) {
		for(j=0;j<size2;j++) {
			r = vec_new(dim);
			s = vec_new(dim);
			matrix[i][j] = folge_new( r, s);
		}
	}
	back->matrix = matrix;
	back->grad1 = grad1;
	back->grad2 = grad2;

	return back;
}




void
folgen_matrix_del(folgen_matrix_p f) {
	int  i, j, size1, size2, dim;
	int  *array1, *array2;

	dim = f->grad1->dim;
	array1 = f->grad1->array;
	array2 = f->grad2->array;
	size1 = 1;
	for(i=0;i<dim;i++) {
		size1 = size1 * (array1[i] + 1);
	}
	size2 = 1;
	for(i=0;i<dim;i++) {
		size2 = size2 * (array2[i] + 1);
	}



	for(i=0;i<size1;i++) {
		for(j=0;j<size2;j++) {
			folge_del(f->matrix[i][j]);
		}
	}
	free(*(f->matrix));
	free(f->matrix);
	vec_del(f->grad1);
	vec_del(f->grad2);
	free(f);
}


void
folgen_vektor_print(folgen_vektor_p f, int info) {
	int  k, dim, size;

	dim = f->grad->dim;
	printf("\n======================================================================");
	printf("\n#Ausgabe eines Folgenvektor der Dimension %d", dim);
	printf("\n\t#grad");
	for(k=0;k<dim;k++){
		printf("\t %d",f->grad->array[k]);
	}
	size = 1;
	for(k=0;k<dim;k++){
		size = size * (f->grad->array[k] + 1);
	}


	if(info==1) {
		for(k=0;k<size;k++) {
			folge_print( f->vektor[k], 0);
		}
	}

	if(info==2) {
		for(k=0;k<size;k++) {
			folge_print( f->vektor[k], 1);
		}
	}
	printf("\n======================================================================\n\n");
	return;
}



folgen_vektor_p
folgen_vektor_build(vec_p p, int a, int n, int mod, bool_t random) {
	vec_p  grad;
	int  dim, k, size;
	folgen_vektor_p  back;

	
	dim = p->dim;
	for(k=0;k<dim;k++) {
		ASSERT( p->array[k] >= 0 );
	}

	grad = vec_new(dim);
	for(k=0;k<dim;k++) {
		if(random == true) {
			grad->array[k] = rand()%(p->array[k] + 1);
		}
		else {
			grad->array[k] = p->array[k];
		}
	}

	back = folgen_vektor_new(grad);
	size = 1;
	for(k=0;k<dim;k++){
		size = size * (grad->array[k] + 1);
	}

	for(k=0;k<size;k++) {
		folge_del(back->vektor[k]);
		back->vektor[k] = folge_build( dim, a, n, mod, random );
	}

	return back;
}



fepc_real_t
folgen_vektor_norm(folgen_vektor_p f, folgen_vektor_p g) {
	vec_p  max, temp, temp1, temp2, r, vec_1;
	int  k, size, test_f, test_g;
	int  i, j, d, dim;
	folge_p  temp_folge;
	fepc_real_t  wert, norm;

	ASSERT( f->grad->dim == g->grad->dim );
	dim = f->grad->dim;
	vec_1 = vec_one(dim);
	temp = vec_max( f->grad, g->grad );
	max = vec_add( temp, vec_1 );
	vec_del(temp);
	size = vec_size( max );
	norm = 0.0;
	for(k=0;k<size;k++) {
		r = entry_one2d( k, max );

		/*Test ob Folge von f den Grad r besitzt*/
		test_f = 0;
		for(d=0;d<dim;d++) {
			if (( r->array[d] ) > ( f->grad->array[d] )) {
				test_f = test_f + 1;
			}
		}
		/*Test ob Folge von g den Grad r besitzt*/
		test_g = 0;
		for(d=0;d<dim;d++) {
			if (( r->array[d] ) > ( g->grad->array[d] )) {
				test_g = test_g + 1;
			}
		}

		/*Berechnung der Norm*/
		if( (test_f==0) && (test_g==0) ) {
			temp = vec_add(f->grad, vec_1);
			i = entry_d2one( r, temp);
			vec_del(temp);

			temp = vec_add(g->grad, vec_1);
			j = entry_d2one( r, temp);
			vec_del(temp);
			wert = folge_norm(f->vektor[i],g->vektor[j]);
			norm = norm + wert;
		}
		if( (test_f==0) && (test_g>0) ) {
			temp = vec_add(f->grad, vec_1);
			i = entry_d2one( r, temp);
			vec_del(temp);

			temp1 = vec_new(dim);
			temp2 = vec_new(dim);
			temp_folge = folge_new( temp1, temp2 );
			wert = folge_norm(f->vektor[i],temp_folge);
			norm = norm + wert;
			folge_del(temp_folge);
		}
		if( (test_f>0) && (test_g==0) ) {
			temp = vec_add(g->grad, vec_1);
			j = entry_d2one( r, temp);
			vec_del(temp);

			temp1 = vec_new(dim);
			temp2 = vec_new(dim);
			temp_folge = folge_new( temp1, temp2 );
			wert = folge_norm(g->vektor[j],temp_folge);
			norm = norm + wert;
			folge_del(temp_folge);
		}

		vec_del( r );
	}

	vec_del(vec_1);
	vec_del(max);

	return norm;
}



folgen_vektor_p
folgen_vektor_faltung(folgen_vektor_p f, folgen_matrix_p Gamma) {
	folgen_vektor_p  back;
	vec_p  vec_1, n_1, n_2, grad_1, grad_2, r, s;
	int  k, j, a, b;
	int  size_1, size_2, dim;
	folge_p  temp, temp1, sum;

	ASSERT( f->grad->dim == Gamma->grad1->dim );
	dim = f->grad->dim;
	vec_1 = vec_one(dim);

	grad_2 = vec_min( f->grad, Gamma->grad2 );
	n_2 = vec_add( grad_2, vec_1 );
	size_2 = vec_size( n_2 );

	grad_1 = vec_copy( Gamma->grad1 );
	n_1 = vec_add( grad_1, vec_1 );
	size_1 = vec_size( n_1 );

	back = folgen_vektor_new( grad_1 );
	for(k=0;k<size_1;k++) {
		sum = folge_new( vec_new(dim), vec_new(dim) );
		for(j=0;j<size_2;j++) {
			r = entry_one2d( j, n_2 );

			s = vec_add( f->grad, vec_1 );
			a = entry_d2one( r, s );
			vec_del( s );

			s = vec_add( Gamma->grad2, vec_1 );
			b = entry_d2one( r, s );
			vec_del( s );
			vec_del( r );

			temp = folge_faltung( f->vektor[a], Gamma->matrix[k][b] );
			temp1 = folge_add( sum, temp );
			folge_del( temp );
			folge_del( sum );

			sum = temp1;
		}
		folge_del( back->vektor[k] );
		back->vektor[k] = sum;
	}

	vec_del( n_1 );
	vec_del( n_2 );
	vec_del( grad_2 );
	vec_del( vec_1 );

	return back;
}



folgen_vektor_p
folgen_vektor_copy_structure(folgen_vektor_p w) {
	folgen_vektor_p  back;
	int  k, dim, size;
	vec_p  grad, start, lang;

	grad = vec_copy( w->grad );
	back = folgen_vektor_new( grad );
	dim = grad->dim;
	size = 1;
	for(k=0;k<dim;k++) {
		size = size * (grad->array[k] + 1);
	}
	for(k=0;k<size;k++) {
		start = vec_copy( w->vektor[k]->start );
		lang = vec_copy( w->vektor[k]->lang );
		folge_del( back->vektor[k] );
		back->vektor[k] = folge_new( start, lang );
	}
	return back;
}



folgen_vektor_p
folgen_vektor_projekt(folgen_vektor_p f,folgen_vektor_p g) {
	folgen_vektor_p  back;
	int  k, size_g, test, dim, d, i;
	vec_p  grad, r, n_f, n_g, start, lang, vec_1;

	ASSERT( f->grad->dim == g->grad->dim );
	dim = f->grad->dim;
	vec_1 = vec_one( dim );

	n_g = vec_add( g->grad, vec_1 );
	n_f = vec_add( f->grad, vec_1 );

	size_g = vec_size( n_g );
	grad = vec_copy( g->grad );
	back = folgen_vektor_new( grad );
	for(k=0;k<size_g;k++) {
		r = entry_one2d( k, n_g );
		test = 0;
		for(d=0;d<dim;d++) {
			if( r->array[d] > f->grad->array[d] ) {
				test = test + 1;
			}
		}
		if(test == 0) {
			i = entry_d2one( r, n_f );
			folge_del( back->vektor[k] );
			back->vektor[k] = folge_projekt( f->vektor[i], g->vektor[k] );
		}
		else {
			folge_del( back->vektor[k] );
			start = vec_copy( g->vektor[k]->start );
			lang = vec_copy( g->vektor[k]->lang );
			back->vektor[k] = folge_new( start, lang );
		}
		vec_del( r );
	}
	vec_del( vec_1 );
	vec_del( n_g );
	vec_del( n_f );

	return back;
}



folgen_vektor_p
folgen_vektor_add(folgen_vektor_p f, folgen_vektor_p g) {
	folgen_vektor_p  back;
	int  k, d, size, dim, a, b, test_f, test_g;
	vec_p  max, vec_1, n, n_f, n_g, r;

	ASSERT( f->grad->dim == g->grad->dim );
	dim = f->grad->dim;
	max = vec_max( f->grad, g->grad );


	vec_1 = vec_one( dim );
	n = vec_add( vec_1, max );
	n_f = vec_add( vec_1, f->grad );
	n_g = vec_add( vec_1, g->grad );
	size = vec_size( n );

	back = folgen_vektor_new( max );
	for(k=0;k<size;k++) {
		r = entry_one2d( k, n );
		test_f = 0;
		test_g = 0;
		for(d=0;d<dim;d++) {
			if( r->array[d] > f->grad->array[d] ) {
				test_f = test_f + 1;
			}
			if( r->array[d] > g->grad->array[d] ) {
				test_g = test_g + 1;
			}
		}
		if( (test_f == 0) && (test_g == 0) ) {
			a = entry_d2one( r, n_f );
			b = entry_d2one( r, n_g );
			folge_del( back->vektor[k] );
			back->vektor[k] = folge_add( f->vektor[a], g->vektor[b] );
		}
		if( (test_f != 0) && (test_g == 0) ) {
			b = entry_d2one( r, n_g );
			folge_del( back->vektor[k] );
			back->vektor[k] = folge_copy( g->vektor[b] );
		}
		if( (test_f == 0) && (test_g != 0) ) {
			a = entry_d2one( r, n_f );
			folge_del( back->vektor[k] );
			back->vektor[k] = folge_copy( f->vektor[a] );
		}
		vec_del( r );
	}
	vec_del( n );
	vec_del( n_f );
	vec_del( n_g );
	vec_del( vec_1 );

	return back;
}

folgen_vektor_p
folgen_vektor_subtract(folgen_vektor_p f, folgen_vektor_p g) {
	folgen_vektor_p  back;
	int  k, d, size, dim, a, b, test_f, test_g;
	vec_p  max, vec_1, n, n_f, n_g, r;

	ASSERT( f->grad->dim == g->grad->dim );
	dim = f->grad->dim;
	max = vec_max( f->grad, g->grad );


	vec_1 = vec_one( dim );
	n = vec_add( vec_1, max );
	n_f = vec_add( vec_1, f->grad );
	n_g = vec_add( vec_1, g->grad );
	size = vec_size( n );

	back = folgen_vektor_new( max );
	
	folge_p temp;
	
	for(k=0;k<size;k++) {
		r = entry_one2d( k, n );
		test_f = 0;
		test_g = 0;
		for(d=0;d<dim;d++) {
			if( r->array[d] > f->grad->array[d] ) {
				test_f = test_f + 1;
			}
			if( r->array[d] > g->grad->array[d] ) {
				test_g = test_g + 1;
			}
		}
		if( (test_f == 0) && (test_g == 0) ) {
			a = entry_d2one( r, n_f );
			b = entry_d2one( r, n_g );
			folge_del( back->vektor[k] );
			back->vektor[k] = folge_subtract( f->vektor[a], g->vektor[b] );
		}
		if( (test_f != 0) && (test_g == 0) ) {
			b = entry_d2one( r, n_g );
			folge_del( back->vektor[k] );
			temp = folge_copy( g->vektor[b] );
			back->vektor[k] =  folge_multi_factor(temp, -1.);
			folge_del(temp); 
		}
		if( (test_f == 0) && (test_g != 0) ) {
			a = entry_d2one( r, n_f );
			folge_del( back->vektor[k] );
			back->vektor[k] = folge_copy( f->vektor[a] );
		}
		vec_del( r );
	}
	vec_del( n );
	vec_del( n_f );
	vec_del( n_g );
	vec_del( vec_1 );

	return back;
}

folgen_vektor_p
folgen_vektor_factor_multi(folgen_vektor_p f, fepc_real_t factor) {
    folgen_vektor_p  back;
	int  k, size;
	vec_p vec_1, n;
	
	vec_1 = vec_one( f->grad->dim );
	n = vec_add( vec_1, f->grad );
	vec_del(vec_1);
	size = vec_size( n );
    vec_del(n);
	back = folgen_vektor_new( f->grad );
	
	
	for(k=0;k<size;k++) {
		back->vektor[k] = folge_multi_factor(f->vektor[k], factor);
	}
	return back;

}



folgen_matrix_p
folgen_matrix_add(folgen_matrix_p f, folgen_matrix_p g) {
	folgen_matrix_p  back;
	int  k, j, d, size1, size2, dim;
	int  a1, a2, b1, b2;
	int  test_f, test_g, test_f1, test_g1, test_f2, test_g2;
	vec_p  grad1, grad2, vec_1, n1, n2;
	vec_p  n_f1, n_g1, n_f2, n_g2, r1, r2;


	ASSERT( f->grad1->dim == g->grad1->dim );
	dim = f->grad1->dim;
	grad1 = vec_max( f->grad1, g->grad1 );
	grad2 = vec_max( f->grad2, g->grad2 );



	vec_1 = vec_one( dim );
	n1 = vec_add( vec_1, grad1 );
	n2 = vec_add( vec_1, grad2 );
	n_f1 = vec_add( vec_1, f->grad1 );
	n_g1 = vec_add( vec_1, g->grad1 );
	n_f2 = vec_add( vec_1, f->grad2 );
	n_g2 = vec_add( vec_1, g->grad2 );
	size1 = vec_size( n1 );
	size2 = vec_size( n2 );

	back = folgen_matrix_new( grad1, grad2 );

	for(k=0;k<size1;k++) {
		r1 = entry_one2d( k, n1 );
		test_f1 = 0;
		test_g1 = 0;
		for(d=0;d<dim;d++) {
			if( r1->array[d] > f->grad1->array[d] ) {
				test_f1 = test_f1 + 1;
			}
			if( r1->array[d] > g->grad1->array[d] ) {
				test_g1 = test_g1 + 1;
			}
		}
		for(j=0;j<size2;j++) {
			r2 = entry_one2d( j, n2 );
			test_f2 = 0;
			test_g2 = 0;
			for(d=0;d<dim;d++) {
				if( r2->array[d] > f->grad2->array[d] ) {
					test_f2 = test_f2 + 1;
				}
				if( r2->array[d] > g->grad2->array[d] ) {
					test_g2 = test_g2 + 1;
				}
			}
			test_g = test_g1 + test_g2;
			test_f = test_f1 + test_f2;
			if( (test_f == 0) && (test_g == 0) ) {
				a1 = entry_d2one( r1, n_f1 );
				a2 = entry_d2one( r2, n_f2 );
				b1 = entry_d2one( r1, n_g1 );
				b2 = entry_d2one( r2, n_g2 );
				folge_del( back->matrix[k][j] );
				back->matrix[k][j] = folge_add( f->matrix[a1][a2], g->matrix[b1][b2] );
			}
			if( (test_f != 0) && (test_g == 0) ) {
				b1 = entry_d2one( r1, n_g1 );
				b2 = entry_d2one( r2, n_g2 );
				folge_del( back->matrix[k][j] );
				back->matrix[k][j] = folge_copy( g->matrix[b1][b2] );
			}
			if( (test_f == 0) && (test_g != 0) ) {
				a1 = entry_d2one( r1, n_f1 );
				a2 = entry_d2one( r2, n_f2 );
				folge_del( back->matrix[k][j] );
				back->matrix[k][j] = folge_copy( g->matrix[a1][a2] );
			}
			vec_del( r2 );
		}
		vec_del( r1 );
	}

	vec_del( n1 );
	vec_del( n2 );
	vec_del( n_f1 );
	vec_del( n_g1 );
	vec_del( n_f2 );
	vec_del( n_g2 );
	vec_del( vec_1 );

	return back;
}























