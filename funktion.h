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

#ifndef __FUNKTION_H
#define __FUNKTION_H

#include "folgen_vektor.h"

typedef struct {
	int  dim;														/* beschreibt zugrundeliegende Dimension */
	int  maxlevel;											/* hoechste Level in Funktionendarstellung */
	folgen_vektor_p  *hierarchie;				/* Array von Zeigern auf die einzelnen Folgenvektoren */
} func_t;

typedef func_t *  func_p;


typedef struct {
	int  dim;														/* beschreibt zugrundeliegende Dimension */
	int  maxlevel;											/* hoechste Level in Funktionendarstellung */
	folgen_matrix_p  *hierarchie;   		/* Array von Zeigern auf die einzelnen Folgenmatritzen */
} func2_t;

typedef func2_t *  func2_p;



/* Funktion wird initialisiert. Beachte: alle Vektoren von Folgen der Funktion werden durch folgen_vektor_new( vec_new(dim) )
initialisiert. f = func_new(maxlevel) bedeutet:
f besteht aus maxlevel+1 Vektoren von Folgen f->hierarchie[k] fuer k=0,...,maxlevel.
Jeder dieser Vektoren von Folgen ist vom grad=vec_new(dim) ,dh. besteht aus einer trivialen Folge.
Somit ist die neu initialisierte Funktionen volldefiniert und vollwertig. */
func_p
func_new(int maxlevel, int dim);

void
func_del(func_p f);

/* Funktion2 wird initialisiert. Beachte: alle Matritzen von Folgen der Funktion werden durch
folgen_matrix_new( vec_new(dim), vec_new(dim) ) initialisiert. f = func2_new(maxlevel) bedeutet:
f besteht aus maxlevel+1 Matritzen von Folgen f->hierarchie[k] fuer k=0,...,maxlevel.
Jeder dieser Matritzen von Folgen ist vom grad1 = vec_new(dim) und grad2 = vec_new(dim) ,dh. besteht aus einer trivialen Folge.
Somit ist die neu initialisierte Funktionen volldefiniert und vollwertig. */
func2_p
func2_new(int maxlevel, int dim);

void
func2_del(func2_p f);


/* Uebertragen der gesamten Folgenwerte der Funtkion f auf eine Funktion mit
 der gleichen Gitterstruktur wie die Funktion g */
func_p
func_projekt(func_p f,func_p g);



/* Gibt die einzelnen Details einer Funktion auf dem Bildschirm wieder. Je groesser der Wert info ist, umso mehr
 Informationen werden wiedergegeben. 0, 1, 2, 3 sind moegliche Infowerte. 0 wenig Info, 3 viele Infos */
void
func_print(func_p f, int info);


/* Erstellt wird eine zufaellige Funktion f der Dimension dim mit folgender Gestalt:
 - 0 <= f->maxlevel <= maxlevel
 - jeder Folgenvektor der Funktion wird durch die Funktion
   folgen_vektor_build(vec_p p, int a, int n, int mod, bool_t random)
   erstellt, wobei fuer den Vektor p gilt: p_i = grad
 - der Paramter random beeinflusst die Zufaelligkeit der Funktion. */
func_p
func_build( int maxlevel, int dim , int grad, int a, int n, int mod, bool_t random);


/* Es werden zwei Funktionen addiert */
func_p
func_add(func_p f, func_p g);

/* Funktion gibt die Anzahl der Freiheitsgrade der Funktionen f, g und w wieder */
int
func_count( func_p f, func_p g, func_p w );


/* Funktion gibt die Anzahl der im Modell verwendeten Freiheitsgrade wieder. Die Definition dieser Zahl ist aus dem 
 zu diesem Programm zugehoerigen Skript zu entnehmen. Die Funktion macht nur dann Sinn, falls der Vektor der 
 zulaessigen Polynomgrade im gesamten Modell gleich bleibt. */
int
func_modell_count( func_p f, func_p g, func_p w );


/* Entsprechend der Gitterstruktur der hp-Funktion werden gegebenenfalls Folgenwerte gleich Null gesetzt. Falls das zu
 einem Folgenwert zugehoerige Intervall weiter verfeinert wird und somit durch Folgenwerte hoeherer
 Level repraesentiert wird, so wird dieser Folgenwert gleich Null gesetzt. */
void
func_grid_zero(func_p f);


#endif
