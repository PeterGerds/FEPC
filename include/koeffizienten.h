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

#ifndef __KOEFFIZIENTEN_H
#define __KOEFFIZIENTEN_H

#include "funktion.h"


typedef struct {
	int  zeilen;									/* Anzahl der Zeilen */
	int  spalten;									/* Anzahl der Spalten */
	fepc_real_t  **a;							/* Matrixglieder a[i][j]; i-te Zeile, j-te Spalte */
} matrix_t;

typedef matrix_t *  matrix_p;

/* Datenstruktur entspicht einer d1 x d2 x d3 Matrix mit den Eintraegen a[i][j][k]
mit 0<=i<d1 , 0<=j<d2 , 0<=k<d3 */
typedef struct {
	int  d1;
	int  d2;
	int  d3;
	fepc_real_t  ***a;				/* Eintraege a[i][j][k] der Matrix */
} matrix3_t;

typedef matrix3_t *  matrix3_p;


/* Initialisieren einer Matrix mit m-Zeilen und n-Spalten */
matrix_p
matrix_new(int m,int n);


void
matrix_del(matrix_p matrix);


/*Initialisieren einer d1 x d2 x d3 Matrix*/
matrix3_p
matrix3_new(int d1,int d2,int d3);


void
matrix3_del(matrix3_p matrix);


/* Algorithmus zur Berechnnung der eindimensionalen Xi-Koeffizienten */
matrix_p
koeffizienten_xi_1dim(int grad);


/* Berechnung der eindimensionalen Gammakoeffizienten mit h = 1.0 und l = 0 */
matrix3_p
koeffizienten_gamma_1dim(int grad);


/* Berechnung der mehrdimensionalen Xi-Koeffizienten bezueglich der Vektoren alpha und q, wobei 0 <= alpha, q. */
fepc_real_t
koeffizienten_xi(vec_p alpha, vec_p q, matrix_p xi_koef);


/* Berechnung der mehrdimensionalen Gamma-Koeffizienten bezueglich des Level level, der Vektoren r, alpha, mu, kappa
 und des Gittervektors mesh. Es gilt 0<=alpha, mu, kappa und -1<=r<=0. */
fepc_real_t
koeffizienten_gamma(int level, vec_p r, vec_real_p mesh, vec_p alpha, vec_p mu, vec_p kappa, matrix3_p gamma_koef);


#endif



