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

#include "faltung.h"



/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/


/*Faltung der Funktion f und g mit Projektion auf die Gitterstruktur von w. Algorithmus A, Fall l'<=l*/
static func_p
faltung_a1( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p Gamma ) {
	func_p  back, omega;
	folgen_vektor_p  vektor, temp1, temp2;
	folgen_matrix_p  matrix;
	int  l, min, dim;


	dim = f->dim;

	/*Omega bauen fuer den Fall A1*/
	if ( f->maxlevel < g->maxlevel ) {
		min = f->maxlevel;
	}
	else {
		min = g->maxlevel;
	}

	omega = func_new( min, dim );
	for(l=min;l>=0;l--) {
		if (l==min) {
			vektor = f->hierarchie[l];
			matrix = Gamma->hierarchie[l];
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = folgen_vektor_faltung( vektor , matrix );
		}
		else {
			vektor = f->hierarchie[l];
			matrix = Gamma->hierarchie[l];
			temp1 = folgen_vektor_faltung( vektor, matrix );
			temp2 = faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef );
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = folgen_vektor_add( temp1, temp2 );

			folgen_vektor_del(temp1);
			folgen_vektor_del(temp2);
		}
	}

	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	func_del( omega );
	return back;
}




/*Faltung der Funktion f und g mit Projektion auf die Gitterstruktur von w. Algorithmus A, Fall l'>l*/
static func_p
faltung_a2( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p Gamma ) {
	func_p  back, omega;
	folgen_vektor_p  vektor, temp1, temp2;
	folgen_matrix_p  matrix;
	int  l, min, dim;


	dim = f->dim;
	/*Omega bauen fuer den Fall A2*/
	if (f->maxlevel < (g->maxlevel-1)) {
		min = f->maxlevel;
	}
	else {
		min = g->maxlevel-1;
	}

	if (min < 0) {
		back = func_new( 0, dim );
		return back;
	}

	omega = func_new( min, dim );
	for(l=min;l>=0;l--) {
		if (l==min) {
			vektor = f->hierarchie[l];
			matrix = Gamma->hierarchie[l];
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = folgen_vektor_faltung( vektor , matrix );
		}
		else {
			vektor = f->hierarchie[l];
			matrix = Gamma->hierarchie[l];
			temp1 = folgen_vektor_faltung( vektor, matrix );
			temp2 = faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef );
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = folgen_vektor_add( temp1, temp2 );

			folgen_vektor_del(temp1);
			folgen_vektor_del(temp2);
		}
	}

	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	func_del( omega );
	return back;
}




/*Faltung der Funktion f und g mit Projektion auf die Gitterstruktur von w. Algorithmus B, Fall l'<=l*/
static func_p
faltung_b1( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p Gamma ) {
	func_p  back, omega, F;
	folgen_vektor_p  vektor;
	folgen_matrix_p  matrix;
	int  l, min, dim;


	dim = f->dim;

	/*F bauen fuer den Fall B1 und zusaetzliches Einkuerzen des Traegers (siehe Dokumentation)*/
	F = kuerzen_F_bauen_b( f, g, w, Gamma, xi_koef );

	/*Omega bauen fuer den Fall B1*/
	if ( w->maxlevel < g->maxlevel ) {
		min = w->maxlevel;
	}
	else {
		min = g->maxlevel;
	}

	omega = func_new( min, dim );
	for(l=1;l<=min;l++) {
		vektor = F->hierarchie[l];
		matrix = Gamma->hierarchie[l];
		folgen_vektor_del( omega->hierarchie[l] );
		omega->hierarchie[l] = folgen_vektor_faltung( vektor, matrix );
	}

	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	func_del( omega );
	func_del( F );

	return back;
}







/*Faltung der Funktion f und g mit Projektion auf die Gitterstruktur von w. Algorithmus B, Fall l'>l*/
static func_p
faltung_b2( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p Gamma ) {
	func_p  back, omega, F;
	func2_p  Gamma2;
	folgen_vektor_p  vektor;
	folgen_matrix_p  matrix, temp1, temp2;
	int  l, min, dim;


	dim = f->dim;
	/*Gamma Daten des Algorithmus A2 anpassen an den Fall B2*/
	Gamma2 = func2_new( g->maxlevel, g->dim);
	for(l=0;l<=g->maxlevel;l++) {
		if( l<g->maxlevel ) {
			temp1 = faltung_hilfe_lambda( g->hierarchie[l], l, mesh, gamma_koef, maxgrad );
			temp2 = folgen_matrix_add( Gamma->hierarchie[l], temp1 );
			folgen_matrix_del( temp1 );
			folgen_matrix_del( Gamma2->hierarchie[l] );
			Gamma2->hierarchie[l] = temp2;
		}
		else {
			temp1 = faltung_hilfe_lambda( g->hierarchie[l], l, mesh, gamma_koef, maxgrad );
			folgen_matrix_del( Gamma2->hierarchie[l] );
			Gamma2->hierarchie[l] = temp1;
		}
	}


	/*F bauen fuer den Fall B2 und zusaetzliches Einkuerzen des Traegers (siehe Dokumentation)*/
	F = kuerzen_F_bauen_b( f, g, w, Gamma2, xi_koef );
	/*Omega bauen fuer den Fall B2*/
	if ( w->maxlevel < g->maxlevel ) {
		min = w->maxlevel;
	}
	else {
		min = g->maxlevel;
	}

	omega = func_new( min, dim );
	for(l=1;l<=min;l++) {
		vektor = F->hierarchie[l];
		matrix = Gamma2->hierarchie[l];
		folgen_vektor_del( omega->hierarchie[l] );
		omega->hierarchie[l] = folgen_vektor_faltung( vektor, matrix );
	}

	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	func2_del( Gamma2 );
	func_del( omega );
	func_del( F );

	return back;
}


/*Faltung der Funktion f und g mit Projektion auf die Gitterstruktur von w. Algorithmus C, Fall l'<=l*/
static func_p
faltung_c1( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef ) {
	func_p  back, omega, F;
	folgen_vektor_p  vektor, temp1, temp2, temp3;
	folgen_matrix_p  matrix;
	support_p  sup;
	func2_p  Gamma;
	int  min, dim, l;


	if ( (w->maxlevel-1) < g->maxlevel) {
		min = w->maxlevel-1;
	}
	else {
		min = g->maxlevel;
	}
	dim = f->dim;
	if ( min < 0 ) {
		back = func_new( 0, dim );
		return back;
	}


	/*Gamma bauen fuer den Fall C1*/
	Gamma = func2_new( min, dim );
	for(l=0;l<=min;l++) {
		temp1 = g->hierarchie[l];
		matrix = faltung_hilfe_lambda_c( temp1, l, mesh, gamma_koef, maxgrad );
		folgen_matrix_del( Gamma->hierarchie[l] );
		Gamma->hierarchie[l] = matrix;
	}

	/*F bauen und Einkuerzen auf den fuer Faltung F*Gamma notwendigen Traeger*/
	F = kuerzen_F_bauen_c1( f, g, w, Gamma, xi_koef );

	/*Berechnen des Definitionsbereiches von Omega, der bei der
	Konstruktion von Omega tatsaechlich benoetigt wird*/
	sup = support_omega( w );

	/*Omega bauen fuer den Fall C1*/
	omega = func_new( w->maxlevel, dim );
	for(l=1;l<=w->maxlevel;l++) {
		if ( l<= (g->maxlevel+1) ) {
			vektor = F->hierarchie[l-1];
			matrix = Gamma->hierarchie[l-1];
			temp1 = folgen_vektor_faltung( vektor, matrix );
			temp2 = omega->hierarchie[l-1];
			temp3 = folgen_vektor_add( temp1, temp2 );
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = faltung_hilfe_prolong( temp3, xi_koef );
			folgen_vektor_del(temp1);
			folgen_vektor_del(temp3);
		}
		else {
			temp1 = omega->hierarchie[l-1];
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = faltung_hilfe_prolong( temp1, xi_koef );
		}

		/*Einkuerzen von Omega auf den Definitionsbereich, der bei der
		Konstruktion von Omega tatsaechlich benoetigt wird*/
		temp1 = omega->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup->a[l], sup->b[l] );
		folgen_vektor_del(temp1);
		omega->hierarchie[l] = temp2;
	}

	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	func2_del( Gamma );
	func_del( omega );
	func_del( F );
	support_del( sup );

	return back;
}



/*Faltung der Funktion f und g mit Projektion auf die Gitterstruktur von w. Algorithmus C, Fall l'>l*/
static func_p
faltung_c2( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef ) {
	func_p  back, omega, F;
	folgen_vektor_p  vektor, temp1, temp2, temp3;
	folgen_matrix_p  matrix;
	support_p  sup;
	func2_p  Gamma;
	int  min, dim, l;


	if ( (w->maxlevel-1) < g->maxlevel) {
		min = w->maxlevel-1;
	}
	else {
		min = g->maxlevel;
	}
	dim = f->dim;
	if ( min < 0 ) {
		back = func_new( 0, dim );
		return back;
	}


	/*Gamma bauen fuer den Fall C2*/
	Gamma = func2_new( min, dim );
	for(l=0;l<=min;l++) {
		temp1 = g->hierarchie[l];
		matrix = faltung_hilfe_lambda_c( temp1, l, mesh, gamma_koef, maxgrad );
		folgen_matrix_del( Gamma->hierarchie[l] );
		Gamma->hierarchie[l] = matrix;
	}

	/*F bauen und Einkuerzen auf den fuer Faltung F*Gamma notwendigen Traeger*/
	F = kuerzen_F_bauen_c2( f, g, w, Gamma, xi_koef );

	/*Berechnen des Definitionsbereiches von Omega, der bei der
	Konstruktion von Omega tatsaechlich benoetigt wird*/
	sup = support_omega( w );

	/*Omega bauen fuer den Fall C1*/
	omega = func_new( w->maxlevel, dim );
	for(l=2;l<=w->maxlevel;l++) {
		if ( l<= (g->maxlevel+1) ) {
			vektor = F->hierarchie[l-1];
			matrix = Gamma->hierarchie[l-1];
			temp1 = folgen_vektor_faltung( vektor, matrix );
			temp2 = omega->hierarchie[l-1];
			temp3 = folgen_vektor_add( temp1, temp2 );
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = faltung_hilfe_prolong( temp3, xi_koef );
			folgen_vektor_del(temp1);
			folgen_vektor_del(temp3);
		}
		else {
			temp1 = omega->hierarchie[l-1];
			folgen_vektor_del( omega->hierarchie[l] );
			omega->hierarchie[l] = faltung_hilfe_prolong( temp1, xi_koef );
		}

		/*Einkuerzen von Omega auf den Definitionsbereich, der bei der
		Konstruktion von Omega tatsaechlich benoetigt wird*/
		temp1 = omega->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup->a[l], sup->b[l] );
		folgen_vektor_del(temp1);
		omega->hierarchie[l] = temp2;
	}

	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	func2_del( Gamma );
	func_del( omega );
	func_del( F );
	support_del( sup );

	return back;
}



/* Berechnung der Algorithmen A1, B1, C1 (Fall l'<=l) */
static func_p
faltung_sum1( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef ) {
	func_p  back, sum, temp1, temp2;
	func2_p  Gamma;

	/*Gamma bauen fuer den Fall A1 und B1*/
	Gamma = faltung_hilfe_Gamma_bauen_1( g, mesh, maxgrad, gamma_koef, xi_koef );

	/*Berechnung eines Teil der Faltung (Fall l'<=l)*/
	sum = faltung_a1( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma );

	temp1 = faltung_b1( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma );
	temp2 = func_add( sum, temp1 );
	func_del( sum );
	func_del( temp1 );
	sum = temp2;

	temp1 = faltung_c1( f, g, w, mesh, maxgrad, gamma_koef, xi_koef );
	temp2 = func_add( sum, temp1 );
	func_del( sum );
	func_del( temp1 );
	sum = temp2;

	func2_del( Gamma );
	back = sum;

	return back;
}



/* Berechnung der Algorithmen A2, B2, C2 (Fall l'>l) */
static func_p
faltung_sum2( func_p f, func_p g, func_p w, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef ) {
	func_p  back, sum, temp1, temp2;
	func2_p  Gamma;

	/*die Funktionen f und g muessen vertauscht werden*/
	temp1 = f;
	f = g;
	g = temp1;

	/*Gamma bauen fuer den Fall A2*/
	Gamma = faltung_hilfe_Gamma_bauen_2( g, mesh, maxgrad, gamma_koef, xi_koef );

	/*Berechnung eines Teil der Faltung (Fall l'<=l)*/
	sum = faltung_a2( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma );

	temp1 = faltung_b2( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma );
	temp2 = func_add( sum, temp1 );
	func_del( sum );
	func_del( temp1 );
	sum = temp2;

	temp1 = faltung_c2( f, g, w, mesh, maxgrad, gamma_koef, xi_koef );
	temp2 = func_add( sum, temp1 );
	func_del( sum );
	func_del( temp1 );
	sum = temp2;

	func2_del( Gamma );
	back = sum;

	return back;
}





/*******************************************************
 *
 * globale Funktionen
 *
 *******************************************************/



func_p
faltung_ref(func_p f, func_p g, func_p w, fepc_real_t h) {
	func_p  back, F, G, omega;
	int  maxlevel, l, k, maxgrad_1dim, dim;
	vec_p  maxgrad, p;
	folgen_matrix_p  Gamma;
	folgen_vektor_p  temp, temp1;
	matrix_p xi_koef;
	matrix3_p gamma_koef;
	vec_real_p  mesh;

	ASSERT( f->dim == g->dim );
	ASSERT( f->dim == w->dim );
	dim = f->dim;

	mesh = vec_real_new( dim );
	for(k=0;k<dim;k++) {
		mesh->array[k]=h;
	}


	/*Berechnen der maximal vorkommenden Grade*/
	maxgrad = vec_new( dim );
	for(k=0;k<=w->maxlevel;k++) {
		p = vec_max( w->hierarchie[k]->grad, maxgrad );
		vec_del( maxgrad );
		maxgrad = p;
	}
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
	/*Berechnung der eindimensionalen Gamma-Koeffizienten*/
	gamma_koef = koeffizienten_gamma_1dim( maxgrad_1dim );


	/*Berechnen des maximal vorkommenden Level*/
	if(f->maxlevel > g->maxlevel) {
		maxlevel = f->maxlevel;
	}
	else {
		maxlevel = g->maxlevel;
	}
	if(w->maxlevel > maxlevel) {
		maxlevel = w->maxlevel;
	}


	/*Verfeinerung der gesamten Funktion f und g auf den hoechsten Level maxlevel derart, dass
	die Funktionen f = F->hierarchie[maxlevel] und g = G->hierarchie[maxlevel] entsprechen*/
	F = func_new( maxlevel, dim );
	G = func_new( maxlevel, dim );


	for(l=0;l<=maxlevel;l++) {
		if (l>0) {
			folgen_vektor_del( F->hierarchie[l] );
			F->hierarchie[l] = faltung_hilfe_prolong( F->hierarchie[l-1], xi_koef );
		}
		if (l<=f->maxlevel) {
			temp = F->hierarchie[l];
			temp1 = folgen_vektor_add(temp,f->hierarchie[l]);

			folgen_vektor_del(temp);

			F->hierarchie[l] = temp1;
		}
	}
	for(l=0;l<=maxlevel;l++) {
		if (l>0) {
			folgen_vektor_del( G->hierarchie[l] );
			G->hierarchie[l] = faltung_hilfe_prolong( G->hierarchie[l-1], xi_koef );
		}
		if (l<=g->maxlevel) {
			temp = G->hierarchie[l];
			temp1 = folgen_vektor_add(temp,g->hierarchie[l]);

			folgen_vektor_del(temp);

			G->hierarchie[l] = temp1;
		}
	}


	/*Implementation einer abgespeckten Version des Teilalgorithmus A, der hinreichend ist
	um die Projektion der Faltung f*g zu berechnen siehe Dokumentation*/
	Gamma = faltung_hilfe_lambda( G->hierarchie[maxlevel], maxlevel, mesh, gamma_koef, maxgrad);
	/*Berechnung der Omegakoeffizienten*/
	omega = func_new( maxlevel, dim );
	for(l=maxlevel;l>=0;l--) {
		if (l==maxlevel) {
			folgen_vektor_del(omega->hierarchie[l]);
			omega->hierarchie[l] = folgen_vektor_faltung( F->hierarchie[l], Gamma );
		}
		else {
			folgen_vektor_del(omega->hierarchie[l]);
			omega->hierarchie[l] = faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef );
		}
	}



	/*Einschraenken der Omegafunktion auf die tatsaechlich vorgegebenen Gitterstruktur von w*/
	back = func_projekt( omega, w );

	vec_del( maxgrad );
	func_del( omega );
	func_del( F );
	func_del( G );
	folgen_matrix_del( Gamma );
	matrix3_del( gamma_koef );
	matrix_del( xi_koef );
	vec_real_del( mesh );

	/*Folgenwerte, die sich auf Intervallen befinden, welche verfeinert werden, werden gleich Null gesetzt*/
	func_grid_zero( back );
	return back;
}



func_p
faltung_fepc(func_p f, func_p g, func_p w, fepc_real_t h) {
	func_p back = func_new(w->maxlevel, w->dim);

	faltung_fepc_overwrite(back, f, g, w, h);
	return back;
}

void
faltung_fepc_overwrite(func_t * result, func_p f, func_p g, func_p w, fepc_real_t h) {
	func_p temp1, temp2;

	int  k, maxgrad_1dim, dim;
	vec_p  maxgrad, p;
	matrix_p xi_koef;
	matrix3_p gamma_koef;
	vec_real_p  mesh;

	ASSERT( f->dim == g->dim );
	ASSERT( f->dim == w->dim );
	dim = f->dim;

	mesh = vec_real_new( dim );
	for(k=0;k<dim;k++) {
		mesh->array[k]=h;
	}

	/*Berechnen der maximal vorkommenden Grade*/
	maxgrad = vec_new( dim );
	for(k=0;k<=w->maxlevel;k++) {
		p = vec_max( w->hierarchie[k]->grad, maxgrad );
		vec_del( maxgrad );
		maxgrad = p;
	}
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
	xi_koef = koeffizienten_xi_1dim( 2 * maxgrad_1dim + 1 );
	/*Berechnung der eindimensionalen Gamma-Koeffizienten*/
	gamma_koef = koeffizienten_gamma_1dim( 2 * maxgrad_1dim + 1 );
	/*Berechnung der Faltung*/
	temp1 = faltung_sum1( f, g, w, mesh, maxgrad, gamma_koef, xi_koef );
	temp2 = faltung_sum2( f, g, w, mesh, maxgrad, gamma_koef, xi_koef );
	func_add_overwrite(result, temp1, temp2 );
	vec_del( maxgrad );
	func_del( temp1 );
	func_del( temp2 );
	matrix3_del( gamma_koef );
	matrix_del( xi_koef );
	vec_real_del( mesh );
	/*Folgenwerte, die sich auf Intervallen befinden, welche verfeinert werden, werden gleich Null gesetzt*/
	func_grid_zero( result );
}




