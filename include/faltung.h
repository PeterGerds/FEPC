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

#ifndef __FALTUNG_H
#define __FALTUNG_H

#include "kuerzen.h"

/* Berechnet die Projektion der Faltung der Funktionen f und g auf die Gitterstruktur von w.
 Von f und g sind die Gitterstrukturen sowie die dazugehoerigen Eintraege gegeben. Von w ist
 nur die Gitterstruktur gegeben.

 Eingabewerte: Funktionen f,g mit Gitterstruktur und Eintraegen (dh. zugehoerigen Folgen)
							Funktion w mit Gitterstruktur ohne Eintraege
 Ausgabewert:	Faltung von f und g, projeziert auf die Gitterstruktur von w */

/* Berechnung der Referenz-Faltung */
func_p
faltung_ref(func_p f, func_p g, func_p w, fepc_real_t h);

/* Berechnung der FEPC-Faltung */
func_p
faltung_fepc(func_p f, func_p g, func_p w, fepc_real_t h);

void
faltung_fepc_overwrite(func_t * result, func_p f, func_p g, func_p w, fepc_real_t h);

#endif
