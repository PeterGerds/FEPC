/*
 * FEPC
 * Copyright (C) 2009 Peter Gerds (gerds@mis.mpg.de), 2011 Stefan Handschuh (handschu@mis.mpg.de)
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

#ifndef __FOLGE_H
#define __FOLGE_H

#include "basic.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	vec_p start;				/* Startvektor der Folge */
	vec_p lang;					/* Laenge der Folge */
	fepc_real_t *glied;	/* Speicherung d. Folgenwerte durch multidimensionales Array dessen Dimension durch Vektor size repraesentiert wird */
} folge_t;

typedef folge_t *  folge_p;

/* Folge wird initialisiert. Beachte: alle Gliedwerte werden gleich 0.0 gesetzt */
folge_p
folge_new(vec_p start,vec_p lang);

/* Folge wird geloescht */
void
folge_del(folge_p f);

/* Faltung der Folgen f und g wird berechnet nach der herkoemmlichen Art und Weise wie die Faltung ist */
folge_p
folge_slow_faltung(folge_p f,folge_p g);

/* Faltung der Folgen f und g wird berechnet ueber die Fouriertransformation */
folge_p
folge_faltung(folge_p f,folge_p g);

/* Gibt die einzelnen Details einer Folge auf dem Bildschirm wieder. Je groesser der Wert info ist,umso
 mehr Informationen werden wiedergegeben. 0, 1, sind moegliche Infowerte. 0 wenig Info, 1 viele Infos */
void
folge_print(folge_p f, int info);

/* Erstellen einer Folge der Dimension dim, deren Eintraege im Intervall [-mod,mod] liegen.
 Die Entraege des Vektors start liegen im Intervall [-a,a] und die des Vektors lang in [1,n].
 Paramter random beeinflusst Zufaelligkeit der Folge. */
folge_p
folge_build(int dim, int a, int n, int mod, bool_t random);

/* Ermittelt ob der Vektor r im Traeger der Folge f liegt. Falls ja, dann Rueckgabe des Folgengliedes
 von f an der Stelle r, falls nein Rueckgabe von 0.0 . */
fepc_real_t
folge_glied( vec_p r, folge_p f );

/* Berechnet die l2-Norm zwischen den Folgen f und g */
fepc_real_t
folge_norm(folge_p f, folge_p g);

/* Addition von Folgen f und g */
folge_p
folge_add(folge_p f, folge_p g);

/* Subtraktion von Folgen f und g */
folge_p
folge_subtract(folge_p f, folge_p g);

/* Rueckgabe einer einer Kopie der Folge f */
folge_p
folge_copy( folge_p f);

/* Uebertragen der gesamten Folgenwerte der Folge f auf eine Folge mit der gleichen Gitterstruktur wie
 die Folge g */
folge_p
folge_projekt(folge_p f, folge_p g);

folge_p
folge_multi_factor(folge_p folge, fepc_real_t factor);

#ifdef __cplusplus
}
#endif

#endif


