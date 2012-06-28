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

#include "faltung_hilfe.h"


/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/


static folge_p
lambda(folgen_vektor_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p alpha, vec_p mu) {
	folge_p  back;
	folge_p  *vektor;
	int  test, size, size_back, size_2;
	vec_p  min, max , anfang, ende, start, lang;
	vec_p  vec_1, vec_2, n, temp, temp1, temp2;
	vec_p  nu, r, kappa;
	int  k, j, t, d, dim;
	fepc_real_t  x, sum, wert;


	dim = g->grad->dim;
	vec_1 = vec_one(dim);
	n = vec_add( g->grad, vec_1 );
	size = vec_size( n );

	vektor = g->vektor;
	/*Bestimmen des Traegers der Folge*/
	min = vec_copy( vektor[0]->start );
	temp = vec_add( min, vektor[0]->lang );
	max = vec_op(1, temp, -1, vec_1);
	vec_del(temp);

	for(k=0;k<size;k++) {
		anfang = vektor[k]->start;
		temp = vec_add( anfang, vektor[k]->lang );
		ende = vec_op(1, temp, -1, vec_1);
		vec_del(temp);

		temp1 = vec_min( min, anfang );
		vec_del(min);
		temp2 = vec_max( max, ende );
		vec_del(max);
		vec_del(ende);

		min = temp1;
		max = temp2;
	}

	start = min;
	temp = vec_op( 1, max, -1, min);
	vec_del(max);
	lang = vec_op( 1, temp, 2, vec_1 );
	vec_del(temp);

	/*Berechnen der Eintraege der Folge*/
	back = folge_new( start, lang );
	size_back = vec_size( lang );

	vec_2 = vec_multi( 2, vec_1 );
	size_2 = vec_size( vec_2 );
	for(j=0;j<size_back;j++) {
		sum = 0.0;
		temp = entry_one2d( j, lang );
		nu = vec_add( start, temp );
		vec_del( temp );
		for(k=0;k<size;k++) {
			kappa = entry_one2d( k, n );

			/*Test ob zugehoeriger Gamma Koeffizient Null ist*/
			test = 0;
			for(d=0;d<dim;d++) {
				if( alpha->array[d] - mu->array[d] -1 > kappa->array[d] ) {
					test = test + 1;
				}
				if( alpha->array[d] + mu->array[d] +1 < kappa->array[d] ) {
					test = test + 1;
				}
				if( mu->array[d] - alpha->array[d] -1 > kappa->array[d] ) {
					test = test + 1;
				}
			}
			if( test==0 ) {
				for(t=0;t<size_2;t++) {
					temp = entry_one2d( t, vec_2 );
					r = vec_multi( -1, temp );
					vec_del( temp );
					temp = vec_add( nu, r );
					wert = folge_glied( temp, vektor[k] );
					vec_del( temp );
					if(wert != 0.0) {
						x = koeffizienten_gamma(level, r, mesh, alpha, mu, kappa, gamma_koef);
						wert = wert * x;
						sum = sum + wert;
					}
					vec_del(r);
				}
			}
			vec_del(kappa);
		}
		vec_del(nu);
		back->glied[j] = sum;
	}

	vec_del( n );
	vec_del( vec_1 );
	vec_del( vec_2 );

	return back;
}






static folge_p
LAMBDA( folgen_matrix_p Gamma, matrix_p xi_koef, vec_p alpha, vec_p mu ) {
	folge_p  back;
	folge_p  **matrix;
	vec_p  n_1, n_2, n_mu, n_alpha, vec_1, vec_2, temp, temp1, temp2;
	int  size_alpha, size_mu, size_2, size_back;
	vec_p  min, max , anfang, ende, start, lang;
	int  i, j, k, a, b, j2, k2;
	vec_p  nu, tau, sigma, r, t;
	fepc_real_t  sum, sub_sum, wert, x, faktor;
	int  skalar_prod1, skalar_prod2, dim;


	dim = mu->dim;
	vec_1 = vec_one(dim);
	n_alpha = vec_add( alpha, vec_1 );
	size_alpha = vec_size( n_alpha );
	n_mu = vec_add( mu, vec_1 );
	size_mu = vec_size( n_mu );
	n_1 = vec_add( Gamma->grad1, vec_1 );
	n_2 = vec_add( Gamma->grad2, vec_1 );

	matrix = Gamma->matrix;
	/*Bestimmen des Traegers der Folge*/
	min = vec_copy( matrix[0][0]->start );
	max = vec_add( min, matrix[0][0]->lang );


	for(k=0;k<size_alpha;k++) {
		temp = entry_one2d( k, n_alpha );
		k2 = entry_d2one( temp, n_1 );
		vec_del( temp );
		for(j=0;j<size_mu;j++) {
			temp = entry_one2d( j, n_mu );
			j2 = entry_d2one( temp, n_2 );
			vec_del( temp );
			anfang = matrix[k2][j2]->start;
			ende = vec_add( anfang, matrix[k2][j2]->lang );

			temp1 = vec_min( min, anfang );
			vec_del(min);
			temp2 = vec_max( max, ende );
			vec_del(max);
			vec_del(ende);

			min = temp1;
			max = temp2;
		}
	}


	temp = vec_op( 1, min, -1, vec_1 );
	vec_del( min );
	start = vec_div( 2, temp );
	vec_del( temp );

	temp1 = vec_div( 2, max );
	vec_del( max );
	temp2 = vec_op( 1, temp1, -1, start );
	vec_del( temp1 );
	lang = vec_add( temp2, vec_1 );
	vec_del( temp2 );

	/*Berechnen der Eintraege der Folge*/
	back = folge_new( start, lang );
	size_back = vec_size( lang );

	vec_2 = vec_multi( 2, vec_1 );
	size_2 = vec_size( vec_2 );
	for(i=0;i<size_back;i++) {
		sum = 0.0;
		temp = entry_one2d( i, lang );
		nu = vec_add( start, temp );
		vec_del( temp );
		for(k=0;k<size_alpha;k++) {
			tau = entry_one2d( k, n_alpha );
			k2 = entry_d2one( tau, n_1 );
			for(j=0;j<size_mu;j++) {
				sigma = entry_one2d( j, n_mu );
				j2 = entry_d2one( sigma, n_2 );

				x = koeffizienten_xi( alpha, tau , xi_koef ) * koeffizienten_xi( mu, sigma, xi_koef );
				sub_sum = 0.0;
				for(a=0;a<size_2;a++) {
					t = entry_one2d( a, vec_2 );
					for(b=0;b<size_2;b++) {
						r = entry_one2d( b, vec_2 );

						temp1 = vec_op( 1, t, -1, r );
						temp2 = vec_op( 2, nu, 1, temp1 );
						vec_del( temp1 );
						wert = folge_glied( temp2, matrix[k2][j2] );
						vec_del( temp2 );
						if( wert != 0.0 ) {
							temp1 = vec_add( tau, alpha);
							temp2 = vec_op( 1, vec_1, -1, t );
							skalar_prod1 = vec_skalar_prod( temp1, temp2 );
							vec_del( temp1 );
							vec_del( temp2 );

							temp1 = vec_add( sigma, mu);
							temp2 = vec_op( 1, vec_1, -1, r );
							skalar_prod2 = vec_skalar_prod( temp1, temp2 );
							vec_del( temp1 );
							vec_del( temp2 );

							faktor = pow( -1, skalar_prod1 + skalar_prod2 );
							sub_sum = sub_sum + (wert * faktor);
						}
						vec_del( r );
					}
					vec_del( t );
				}
				sum = sum + (x * sub_sum);
				vec_del( sigma );
			}
			vec_del( tau );
		}
		vec_del( nu );
		back->glied[i] = sum;
	}
	vec_del( n_1 );
	vec_del( n_2 );
	vec_del( n_alpha );
	vec_del( n_mu );
	vec_del( vec_1 );
	vec_del( vec_2 );

	return back;
}





/*******************************************************
 *
 * globale Funktionen
 *
 *******************************************************/




folgen_vektor_p
faltung_hilfe_prolong(folgen_vektor_p f,matrix_p xi_koef) {
	folgen_vektor_p  back;
	folge_p  *vektor;
	folge_p  g;
	vec_p  grad, n, vec_1;
	vec_p  min, max , anfang, ende, start, lang;
	vec_p  temp, temp1, temp2;
	vec_p  rest;
	int  size, size_g;
	vec_p  alpha, q, nu, r;
	vec_p  *array;
	fepc_real_t  sum, prod , wert, x;
	int  k, i, j, c, dim, a;


	grad = vec_copy( f->grad );
	dim = grad->dim;
	vec_1 = vec_one( dim );

	n = vec_add( grad , vec_1 );
	size = vec_size( n );

	vektor = f->vektor;
	back = folgen_vektor_new( grad );

	/*Berechnungen*/
	for(k=0;k<size;k++) {

		/*Bestimmen des Traegers der Folge*/
		min = vec_copy(vektor[k]->start);
		temp = vec_add( min, vektor[k]->lang);
		max = vec_op(1, temp, -1, vec_1);
		vec_del(temp);

		q = entry_one2d( k, n);
		rest = vec_r_s_n( q, n );

		for(i=0;i<rest->dim;i++) {
			c = rest->array[i];
			anfang = vektor[c]->start;
			temp = vec_add( anfang, vektor[c]->lang);
			ende = vec_op(1, temp, -1, vec_1);
			vec_del(temp);

			temp1 = vec_min( min, anfang );
			vec_del(min);
			temp2 = vec_max( max, ende );
			vec_del(max);
			vec_del(ende);

			min = temp1;
			max = temp2;
		}

		start = vec_multi( 2, min);
		vec_del(min);
		temp = vec_op( 2 , max , -1 ,start );
		vec_del(max);
		lang = vec_op( 1, temp , 2, vec_1 );
		vec_del(temp);
		g = folge_new( start, lang );


		/*Berechnen der Eintraege der Folge*/
		size_g = vec_size( lang );
		for(j=0;j<size_g;j++) {
			sum = 0.0;
			temp = entry_one2d( j, lang);
			temp1 = vec_add( temp, start );
			vec_del(temp);
			array = vec_zerlegung( temp1 );
			vec_del(temp1);
			nu = array[0];
			r = array[1];

			for(i=0;i<rest->dim;i++) {
				c = rest->array[i];
				alpha = entry_one2d( c, n );

				/*entsprechendes Folgenglied berechnen*/
				wert = folge_glied( nu , vektor[c] );
				if( wert != 0.0 ) {
					prod = wert;

					/*Skalarprodukt berechnen*/
					temp1 = vec_add( q , alpha );
					temp2 = vec_op( 1, vec_1, -1, r );
					a = vec_skalar_prod( temp1, temp2 );
					vec_del(temp1);
					vec_del(temp2);
					prod = prod * pow( -1, a );

					/*Xi-Koeffizienten*/
					x = koeffizienten_xi( alpha, q, xi_koef );
					prod = prod * x;

					sum = sum + prod;
				}
				vec_del(alpha);
			}
			vec_del(r);
			vec_del(nu);
			free(array);
			g->glied[j] = sum;
		}
		vec_del(q);
		vec_del(rest);
		folge_del( back->vektor[k] );
		back->vektor[k] = g;
	}
	vec_del(vec_1);
	vec_del(n);
	return back;
}





folgen_vektor_p
faltung_hilfe_restrict(folgen_vektor_p f, matrix_p xi_koef) {
	folgen_vektor_p  back;
	folge_p  *vektor;
	folge_p  g;
	vec_p  grad, n, vec_1, vec_2;
	vec_p  min, max , anfang, ende, start, lang;
	vec_p  temp, temp1, temp2;
	vec_p  rest;
	int  size, size_g, size_2;
	vec_p  alpha, q, nu, r, s;
	vec_p  *array;
	fepc_real_t  sum, prod , wert, x;
	int  t, k, i, j, c, dim, a;



	grad = vec_copy( f->grad );
	dim = grad->dim;
	vec_1 = vec_one( dim );
	vec_2 = vec_multi( 2, vec_1 );
	size_2 = vec_size(vec_2);
	
	n = vec_add( grad , vec_1 );
	size = vec_size( n );

	vektor = f->vektor;
	back = folgen_vektor_new( grad );



	/*Berechnungen*/
	for(k=0;k<size;k++) {

		/*Bestimmen des Traegers der Folge*/
		min = vec_copy(vektor[k]->start);
		temp = vec_add( min, vektor[k]->lang);
		max = vec_op(1, temp, -1, vec_1);
		vec_del(temp);

		q = entry_one2d( k, n);
		rest = vec_0_s_r( q, n );

		for(i=0;i<rest->dim;i++) {
			c = rest->array[i];
			anfang = vektor[c]->start;
			temp = vec_add( anfang, vektor[c]->lang);
			ende = vec_op(1, temp, -1, vec_1);
			vec_del(temp);

			temp1 = vec_min( min, anfang );
			vec_del(min);
			temp2 = vec_max( max, ende );
			vec_del(max);
			vec_del(ende);

			min = temp1;
			max = temp2;
		}

		temp = vec_op( 1, min, -1, vec_1 );
		vec_del(min);
		start = vec_div( 2 , temp );
		vec_del(temp);
		temp = vec_div( 2 , max );
		vec_del(max);
		temp1 = vec_op( 1, temp, -1 , start);
		vec_del(temp);
		lang = vec_add( temp1, vec_1 );
		vec_del(temp1);
		g = folge_new( start, lang );

		/*Berechnen der Eintraege der Folge*/
		size_g = vec_size( lang );
		for(j=0;j<size_g;j++) {
			sum = 0.0;
			temp = entry_one2d( j, lang);
			nu = vec_add( temp, start );
			vec_del(temp);
			for(i=0;i<rest->dim;i++) {
				c = rest->array[i];
				alpha = entry_one2d( c, n );

				for(t=0;t<size_2;t++) {
					r = entry_one2d( t, vec_2 );
					/*entsprechendes Folgenglied berechnen*/
					s = vec_op( 2, nu, 1, r);
					wert = folge_glied( s , vektor[c] );
					if( wert != 0.0 ) {
						prod = wert;

						/*Skalarprodukt berechnen*/
						temp1 = vec_add( q , alpha );
						temp2 = vec_op( 1, vec_1, -1, r );
						a = vec_skalar_prod( temp1, temp2 );
						vec_del(temp1);
						vec_del(temp2);
						prod = prod * pow( -1, a );

						/*Xi-Koeffizienten*/
						x = koeffizienten_xi( q, alpha, xi_koef );
						prod = prod * x;

						sum = sum + prod;
					}
					vec_del(r);
					vec_del(s);
				}
				vec_del(alpha);
			}
			vec_del(nu);
			g->glied[j] = sum;
		}
		vec_del(q);
		vec_del(rest);
		folge_del( back->vektor[k] );
		back->vektor[k] = g;
	}
	vec_del(vec_1);
	vec_del(vec_2);
	vec_del(n);
	return back;
}




folgen_matrix_p
faltung_hilfe_lambda(folgen_vektor_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p maxgrad) {
	folgen_matrix_p  back;
	int  j, k, size, max, temp, dim;
	vec_p  vec_1, n, alpha, mu, grad1, grad2;

	ASSERT( maxgrad->dim == g->grad->dim);
	dim = g->grad->dim;
	ASSERT( dim == mesh->dim );

	max = 0;
	for(k=0;k<dim;k++) {
		temp = maxgrad->array[k];
		if( temp > max ) {
			max = temp;
		}
		temp = g->grad->array[k];
		if( temp > max ) {
			max = temp;
		}
	}

	ASSERT( max >= 0 );
	ASSERT( gamma_koef->d1 > max );

	vec_1 = vec_one( dim );
	n = vec_add( maxgrad, vec_1 );
	for(k=0;k<dim;k++) {
		ASSERT( n->array[k] > 0 );
	}
	size = vec_size( n );


	grad1 = vec_copy( maxgrad );
	grad2 = vec_copy( maxgrad );
	back = folgen_matrix_new( grad1, grad2 );

	for(k=0;k<size;k++) {
		alpha = entry_one2d( k, n );
		for(j=0;j<size;j++) {
			mu = entry_one2d( j, n );
			folge_del( back->matrix[k][j] );
			back->matrix[k][j] = lambda( g, level, mesh, gamma_koef, alpha, mu );
			vec_del( mu );
		}
		vec_del( alpha );
	}

	vec_del( n );
	vec_del( vec_1 );
	return back;
}



folgen_matrix_p
faltung_hilfe_LAMBDA(folgen_matrix_p Gamma, matrix_p xi_koef) {
	folgen_matrix_p  back;
	vec_p  grad1, grad2, vec_1, n, alpha, mu;
	int  dim, k, j, temp, max, size;


	grad1 = vec_copy( Gamma->grad1 );
	grad2 = vec_copy( Gamma->grad2 );

	dim = grad1->dim;
	for(k=0;k<dim;k++) {
		ASSERT( grad1->array[k] == grad2->array[k] );
	}

	max = 0;
	for(k=0;k<dim;k++) {
		temp = grad1->array[k];
		if( temp > max ) {
			max = temp;
		}
	}
	ASSERT( xi_koef->zeilen > max );



	vec_1 = vec_one(dim);
	n = vec_add( grad1, vec_1 );
	size = vec_size( n );


	back = folgen_matrix_new( grad1, grad2 );

	for(k=0;k<size;k++) {
		alpha = entry_one2d( k, n );
		for(j=0;j<size;j++) {
			mu = entry_one2d( j, n );
			folge_del( back->matrix[k][j] );
			back->matrix[k][j] = LAMBDA( Gamma, xi_koef, alpha, mu );
			vec_del( mu );
		}
		vec_del( alpha );
	}

	vec_del( n );
	vec_del( vec_1 );
	return back;
}



fepc_real_t
faltung_hilfe_norm(func_p f, func_p g) {
	func_p  f_dach, g_dach;
	int  maxlevel, l, k, maxgrad_1dim, dim;
	vec_p  maxgrad, p;
	folgen_vektor_p  temp, temp1;
	matrix_p  xi_koef;
	fepc_real_t  norm;


	ASSERT( f->dim == g->dim );
	dim = f->dim;

	/*Berechnen der maximal vorkommenden Grade*/
	maxgrad = vec_new( dim );
	for(k=0;k<=f->maxlevel;k++) {
		p = vec_max( f->hierarchie[k]->grad, maxgrad );
		vec_del( maxgrad );
		maxgrad = p;
	}
	for(k=0;k<=g->maxlevel;k++) {
		p = vec_max( g->hierarchie[k]->grad, maxgrad );
		vec_del( maxgrad );
		maxgrad = p;
	}

	maxgrad_1dim = 0;
	for(k=0;k<dim;k++) {
		if( maxgrad_1dim < maxgrad->array[k] ) {
			maxgrad_1dim = maxgrad->array[k];
		}
	}

	/*Berechnung der eindimensionalen Xi-Koeffizienten*/
	xi_koef = koeffizienten_xi_1dim( maxgrad_1dim );

	/*Berechnen des maximal vorkommenden Level*/
	if(f->maxlevel > g->maxlevel) {
		maxlevel = f->maxlevel;
	}
	else {
		maxlevel = g->maxlevel;
	}



	/*Verfeinerung der gesamten Funktion f und g auf den hoechsten Level maxlevel derart, dass
	die Funktionen f = f_dach->hierarchie[maxlevel] und g = g_dach->hierarchie[maxlevel] entsprechen*/
	f_dach = func_new( maxlevel, dim );
	g_dach = func_new( maxlevel, dim );


	for(l=0;l<=maxlevel;l++) {
		if (l>0) {
			folgen_vektor_del(f_dach->hierarchie[l]);
			f_dach->hierarchie[l] = faltung_hilfe_prolong( f_dach->hierarchie[l-1], xi_koef );
		}
		if (l<=f->maxlevel) {
			temp = f_dach->hierarchie[l];
			temp1 = folgen_vektor_add(temp,f->hierarchie[l]);

			folgen_vektor_del(temp);

			f_dach->hierarchie[l] = temp1;
		}
	}
	for(l=0;l<=maxlevel;l++) {
		if (l>0) {
			folgen_vektor_del(g_dach->hierarchie[l]);
			g_dach->hierarchie[l] = faltung_hilfe_prolong( g_dach->hierarchie[l-1], xi_koef );
		}
		if (l<=g->maxlevel) {
			temp = g_dach->hierarchie[l];
			temp1 = folgen_vektor_add(temp,g->hierarchie[l]);

			folgen_vektor_del(temp);

			g_dach->hierarchie[l] = temp1;
		}
	}



	/*Berechnung der l2_Norm zwischen f_dach->hierarchie[maxlevel] und
	g_dach->hierarchie[maxlevel]*/
	temp1 = f_dach->hierarchie[maxlevel];
	temp = g_dach->hierarchie[maxlevel];

	norm = folgen_vektor_norm( temp1, temp );

	vec_del( maxgrad );
	func_del( f_dach );
	func_del( g_dach );
	matrix_del( xi_koef );

	return norm;
}





func2_p
faltung_hilfe_Gamma_bauen_1( func_p g, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef ) {
	func2_p  back;
	int  l, maxlevel;
	folgen_matrix_p  temp, temp1, temp2;

	maxlevel = g->maxlevel;
	back = func2_new( maxlevel, g->dim );
	for(l=maxlevel;l>=0;l--) {
		if ( l == maxlevel ) {
			temp = faltung_hilfe_lambda( g->hierarchie[l], l, mesh, gamma_koef, maxgrad );
			folgen_matrix_del( back->hierarchie[l] );
			back->hierarchie[l] = temp;
		}
		else {
			temp1 = faltung_hilfe_LAMBDA( back->hierarchie[l+1], xi_koef );
			temp2 = faltung_hilfe_lambda( g->hierarchie[l], l, mesh, gamma_koef, maxgrad );
			temp = folgen_matrix_add(temp1,temp2);

			folgen_matrix_del( temp1 );
			folgen_matrix_del( temp2 );

			folgen_matrix_del( back->hierarchie[l] );
			back->hierarchie[l] = temp;
		}
	}

	return back;
}






func2_p
faltung_hilfe_Gamma_bauen_2( func_p g, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef ) {
	func2_p  back;
	int  l, maxlevel;
	folgen_matrix_p  temp, temp1, temp2, temp3;

	maxlevel = g->maxlevel - 1;

	/*Dieser Spezialfall muss an dieser Stelle abgefangen werden,da das Ergebnis
	in diesem Fall weiterverwendet wird (siehe Funktion faltung_sum2)*/
	if(maxlevel == -1) {
		back = func2_new( 0, g->dim );
		return back;
	}
	
	back = func2_new( maxlevel, g->dim );
	for(l=maxlevel;l>=0;l--) {
		if ( l == maxlevel ) {
			temp1 = faltung_hilfe_lambda( g->hierarchie[l+1], l+1, mesh, gamma_koef, maxgrad );
			temp2 = faltung_hilfe_LAMBDA( temp1, xi_koef );
			folgen_matrix_del( temp1 );
			folgen_matrix_del( back->hierarchie[l] );
			back->hierarchie[l] = temp2;
		}
		else {
			temp1 = back->hierarchie[l+1];
			temp2 = faltung_hilfe_lambda( g->hierarchie[l+1], l+1, mesh, gamma_koef, maxgrad );
			temp3 = folgen_matrix_add( temp1, temp2 );
			temp = faltung_hilfe_LAMBDA( temp3, xi_koef );

			folgen_matrix_del(temp2);
			folgen_matrix_del(temp3);

			folgen_matrix_del( back->hierarchie[l] );
			back->hierarchie[l] = temp;
		}
	}

	return back;
}



folgen_matrix_p
faltung_hilfe_lambda_c(folgen_vektor_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p maxgrad) {
	folgen_matrix_p  back;
	int  j, k, size_1, size_2, max, temp, dim;
	vec_p  vec_1, n_1, n_2, alpha, mu, grad1, grad2;

	ASSERT( maxgrad->dim == g->grad->dim);
	dim = g->grad->dim;
	ASSERT( dim == mesh->dim );

	max = 0;
	for(k=0;k<dim;k++) {
		temp = maxgrad->array[k];
		if( temp > max ) {
			max = temp;
		}
		temp = g->grad->array[k];
		if( temp > max ) {
			max = temp;
		}
	}

	ASSERT( max >= 0 );
	ASSERT( gamma_koef->d1 > 2*max+1 );
	vec_1 = vec_one(dim);

	grad2 = vec_copy( maxgrad );
	n_2 = vec_add( grad2, vec_1 );
	size_2 = vec_size( n_2 );
	for(k=0;k<dim;k++) {
		ASSERT( n_2->array[k] > 0 );
	}

	grad1 = vec_op( 2, maxgrad, 1, vec_1 );
	n_1 = vec_add( grad1, vec_1 );
	size_1 = vec_size( n_1 );

	back = folgen_matrix_new( grad1, grad2 );
	for(k=0;k<size_1;k++) {
		alpha = entry_one2d( k, n_1 );
		for(j=0;j<size_2;j++) {
			mu = entry_one2d( j, n_2 );
			folge_del( back->matrix[k][j] );
			back->matrix[k][j] = lambda( g, level, mesh, gamma_koef, alpha, mu );
			vec_del( mu );
		}
		vec_del( alpha );
	}

	vec_del( n_1 );
	vec_del( n_2 );
	vec_del( vec_1 );
	return back;
}



/*******************************************************
 *
 * alternative Funktionen
 *
 *******************************************************/


/*
static folge_p
lambda_sum( folge_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p alpha, vec_p mu, vec_p kappa) {
	folge_p  back;
	int  k, dim, t, size, size_2;
	vec_p  start, lang, vec_1, vec_2;
	vec_p  temp, r, nu;
	fepc_real_t  x, sum, wert;


	dim = g->lang->dim;
	vec_1 = vec_one(dim);
	start = vec_copy( g->start );
	lang = vec_add( g->lang, vec_1 );


	back = folge_new( start, lang );
	size = vec_size( lang );
	vec_2 = vec_multi( 2, vec_1 );
	size_2 = vec_size( vec_2 );

	for(k=0;k<size;k++) {
		sum = 0.0;
		temp = entry_one2d( k, lang );
		nu = vec_add( start, temp );
		vec_del( temp );
		for(t=0;t<size_2;t++) {
			temp = entry_one2d( t, vec_2 );
			r = vec_multi( -1, temp );
			vec_del( temp );
			temp = vec_add( nu, r );
			wert = folge_glied( temp, g );
			vec_del( temp );
			if( wert != 0.0 ) {
				x = koeffizienten_gamma( level, r, mesh, alpha, mu, kappa, gamma_koef );
				wert = wert * x;
				sum = sum + wert;
			}
			vec_del(r);
		}
		vec_del( nu );
		back->glied[k] = sum;
	}

	vec_del( vec_1 );
	vec_del( vec_2 );

	return back;
}

static folge_p
lambda(folgen_vektor_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p alpha, vec_p mu) {
	folge_p  back, temp1, temp2;
	vec_p  n, vec_1, kappa;
	int  k, size, test, d, dim;


	dim = g->grad->dim;
	vec_1 = vec_one(dim);
	n = vec_add( g->grad, vec_1 );
	size = vec_size( n );

	back = folge_new( vec_new(dim), vec_new(dim) );
	for(k=0;k<size;k++) {
		kappa = entry_one2d( k, n );
		test = 0;
		for(d=0;d<dim;d++) {
			if( alpha->array[d] - mu->array[d] -1 > kappa->array[d] ) {
				test = test + 1;
			}
			if( alpha->array[d] + mu->array[d] +1 < kappa->array[d] ) {
				test = test + 1;
			}
			if( mu->array[d] - alpha->array[d] -1 > kappa->array[d] ) {
				test = test + 1;
			}
		}
		if( test == 0 ) {
			temp1 = lambda_sum( g->vektor[k], level, mesh, gamma_koef, alpha, mu, kappa );
			temp2 = folge_add( back, temp1 );
			folge_del( back );
			folge_del( temp1 );
			back = temp2;
		}
		vec_del( kappa );
	}

	vec_del( vec_1 );
	vec_del( n );
	return back;
}
*/

