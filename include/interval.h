#ifndef __INTERVALS
#define __INTERVALS

#include "basic.h"

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

#endif

