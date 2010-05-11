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

#ifndef __FALTUNG_HILFE_H
#define __FALTUNG_HILFE_H

#include "koeffizienten.h"


/* Berechnung der Prolongation siehe Dokumentation */
folgen_vektor_p
faltung_hilfe_prolong(folgen_vektor_p f,matrix_p xi_koef);

/* Berechnung der Restriktion siehe Dokumentation */
folgen_vektor_p
faltung_hilfe_restrict(folgen_vektor_p f, matrix_p xi_koef);

/* Funktion zur Berechnung von Gamma[l][l], siehe Dokumentation */
folgen_matrix_p
faltung_hilfe_lambda(folgen_vektor_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p maxgrad);

/* Funktion zur Berechnung von Gamma_dach[l][l] fuer den Fall C, siehe Dokumentation */
folgen_matrix_p
faltung_hilfe_lambda_c(folgen_vektor_p g, int level, vec_real_p mesh, matrix3_p gamma_koef, vec_p maxgrad);

/* Funktion zur Berechnung von Gamma[l], siehe Dokumentation */
folgen_matrix_p
faltung_hilfe_LAMBDA(folgen_matrix_p Gamma, matrix_p xi_koef);


/* Berechnet eine Norm zwischen 2 Funktionen. Dies geschieht, indem beide Funktionen auf deren hoechste Stufe
 verfeinert werden um anschliessend auf die entstandenen Zahlenfolgen die l2 Norm anzuwenden */
fepc_real_t
faltung_hilfe_norm(func_p f, func_p g);


/* Funktion zur Berechnung der Gamma Daten fuer den Algorithmus A1 und B1 (Fall l'<=l siehe Dokumentation) */
func2_p
faltung_hilfe_Gamma_bauen_1( func_p g, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef );

/* Funktion zur Berechnung der Gamma Daten fuer den Algorithmus A2 (Fall l'>l siehe Dokumentation) */
func2_p
faltung_hilfe_Gamma_bauen_2( func_p g, vec_real_p mesh, vec_p maxgrad, matrix3_p gamma_koef, matrix_p xi_koef );


#endif










