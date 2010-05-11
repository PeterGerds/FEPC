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

#ifndef __FOLGEN_VEKTOR_H
#define __FOLGEN_VEKTOR_H

#include "folge.h"

typedef struct {
	vec_p  grad;				/*Polynomgrad des Folgenvektors*/
	folge_p  *vektor;		/*Speicherung von Folgen durch multidimensionales Array dessen Dimension durch Vektor grad repraesentiert wird*/
} folgen_vektor_t;

typedef folgen_vektor_t *  folgen_vektor_p;


typedef struct {
	vec_p  grad1;						/*Polynomgrad der 1.Komponente*/
	vec_p  grad2;						/*Polynomgrad der 2.Komponente*/
	folge_p  **matrix;			/*Speicherung von Folgen durch multidimensionale Matrix dessen Dimension durch Vektor grad1 und grad2
	repraesentiert wird. Ausfuehrliche Erlaeuterung ist in Dokumentation zu finden*/
} folgen_matrix_t;

typedef folgen_matrix_t *  folgen_matrix_p;



/*folgen_vektor wird initialisiert. Beachte: alle Folgen des Vektors werden durch folge_new(r,r) initialisiert, wobei
r = vec_new(dim) ist. Das heisst alle Folgen des Vektors haben die Laenge 0 (siehe Implementierung von folge_new)*/
folgen_vektor_p
folgen_vektor_new(vec_p grad);

void
folgen_vektor_del(folgen_vektor_p f);


/*folgen_matrix wird initialisiert. Beachte: alle Folgen der Matrix werden durch folge_new(r,r) initialisiert, wobei
 r = vec_new(dim) ist. Das heisst alle Folgen der Matrix haben die Laenge 0 (siehe Implementierung von folge_new)*/
folgen_matrix_p
folgen_matrix_new(vec_p grad1, vec_p grad2);

void
folgen_matrix_del(folgen_matrix_p f);



/* Gibt die einzelnen Details eines Folgenvektor auf dem Bildschirm wieder. Je groesser der Wert info ist, umso mehr
 Informationen werden wiedergegeben. 0, 1, 2 sind moegliche Infowerte. 0 wenig Info, 2 viele Infos */
void
folgen_vektor_print(folgen_vektor_p f, int info);

/* Erstellen eines Folgenvektor fuer dessen grad gilt: 0<=grad<=p. Die Folgen des Vektors werden durch
 die Funktion folge_build(int dim, int a, int n, int mod, bool_t random) erstellt. Paramter random
 beeinflusst Zufaelligkeit des Folgenvektors. */
folgen_vektor_p
folgen_vektor_build(vec_p p, int a, int n, int mod, bool_t random);


/* Berechnet eine Norm zwischen zwei Folgenvektoren. Hierbei werden einfach die l2-Normen zwischen den
 einzelnen Folgen addiert */
fepc_real_t
folgen_vektor_norm(folgen_vektor_p f, folgen_vektor_p g);


/* Berechnung der erweiterten Version der diskreten Faltung. Faltung von Tupeln von Folgen siehe Dokumentation */
folgen_vektor_p
folgen_vektor_faltung(folgen_vektor_p f, folgen_matrix_p Gamma);


/* Es wird ein Vektor von Folgen erstellt, der die gleiche Struktur wie der Eingabevektor besitzt. Alle Folgen des
 Vektors haben als Folgenglied den Wert 0.0 */
folgen_vektor_p
folgen_vektor_copy_structure(folgen_vektor_p w);

/* Uebertragen der gesamten Folgenwerte des Folgenvektors f auf einen Folgenvektor mit
 der gleichen Gitterstruktur wie der Folgenvektor g */
folgen_vektor_p
folgen_vektor_projekt(folgen_vektor_p f,folgen_vektor_p g);

/* Addition von Folgenvektor f und g */
folgen_vektor_p
folgen_vektor_add(folgen_vektor_p f, folgen_vektor_p g);

/* Subtraktion von Folgenvektor f und g */
folgen_vektor_p
folgen_vektor_subtract(folgen_vektor_p f, folgen_vektor_p g);

/* Multiplikation von Folgenvektor f mit faktor a */
folgen_vektor_p
folgen_vektor_factor_multi(folgen_vektor_p f, fepc_real_t a);

/* Addition von Folgenmatritzen f und g */
folgen_matrix_p
folgen_matrix_add(folgen_matrix_p f, folgen_matrix_p g);

#endif


