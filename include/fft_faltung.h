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

#ifndef __FFT_FALTUNG_H
#define __FFT_FALTUNG_H

#include "basic.h"

#ifdef __cplusplus
extern "C" {
#endif

/*Berechnung der Faltung von zwei Arrays ueber die Fouriertransformation*/
fepc_real_t*
fft_faltung(fepc_real_t* a, vec_p n_a, fepc_real_t* b, vec_p n_b);

#ifdef __cplusplus
}
#endif

#endif
