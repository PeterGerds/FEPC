/*
 * FEPC
 * Copyright (C) 2010, 2011 Stefan Handschuh (handschu@mis.mpg.de)
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

#ifndef __INTERVALS
#define __INTERVALS

#include "basic.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int dimension;
    fepc_real_t* start;
    fepc_real_t* end;
} interval_t; 

typedef interval_t * interval_p;

/**
 * Creates a new interval with the given dimension.
 */
interval_p 
interval_new(int dimension);

void
interval_init(interval_p interval, int dimension);

void
interval_del(interval_p interval);

void
interval_clear(interval_t * interval);

void
intervals_del(interval_p* intervals, int interval_count);

void
intervals_del_type(interval_t * intervals, int interval_count);

/**
 * Prints out the range of the given interval.
 */
void 
print_interval(interval_p interval);

/**
 * Prints out the ranges of all given intervals.
 */
void 
print_intervals(interval_p * intervals, int count);

interval_p 
interval_clone(interval_p interval);

interval_p *
intervals_clone(interval_p *interval, int count);

#ifdef __cplusplus
}
#endif

#endif

