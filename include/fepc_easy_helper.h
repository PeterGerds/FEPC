/*
 * FEPC
 * Copyright (C) 2010 Stefan Handschuh (handschu@mis.mpg.de)
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
 
#ifndef _FEPC_EASY_HELPER
#define _FEPC_EASY_HELPER

#include "fepc_easy.h"

func_p
func_multi(func_p f, func_p g, fepc_real_t stepping);

func_p
func_div(func_p f, func_p g, fepc_real_t stepping);

func_p
func_reflect(func_p f);

func_p
func_reflect_x(func_p f);

folgen_vektor_p
folgen_vektor_multi(folgen_vektor_p f, folgen_vektor_p g, int step, fepc_real_t stepping);

folgen_vektor_p
folgen_vektor_div(folgen_vektor_p f, folgen_vektor_p g, int step, fepc_real_t stepping);

fepc_real_t folge_norm2(folge_p folge, int step, fepc_real_t stepping);

folge_p
folge_multi(folge_p f, folge_p g, int step, fepc_real_t stepping);

folge_p
folge_div(folge_p f, folge_p g, int step, fepc_real_t stepping);

vec_real_p 
get_mean_points(vec_p v, vec_p grad, int step, fepc_real_t stepping);

void
func_add2(func_p f, func_p g);

void
func_add3(func_p f, fepc_real_t factor, func_p g);


func_p
func_subtract(func_p f, func_p g);

void
func_subtract2(func_p f, func_p g);

/**
 * Calculates the square of L2 norm of the given function.
 */
fepc_real_t
func_norm_l2_sqr(func_p function, fepc_real_t stepping);


#endif
