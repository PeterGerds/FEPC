/*
 * Copyright (C) 2010 Philipp Waehnert (waehnert@mis.mpg.de), 2011 Stefan Handschuh (handschu@mis.mpg.de)
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

#ifndef __gauss_h
#define __gauss_h

#ifdef __cplusplus
extern "C" {
#endif

typedef struct node_ node;
struct node_ { 
  double x;
  double w;
};

node * create_gauss_nodes(int n);
int node_compare(const void* a, const void* b);
void sort_gauss_nodes(node* nodes, unsigned int n);
void shift_gauss_nodes(node* nodes, unsigned int n, double a, double b);
void free_gauss_nodes(node* nodes);

#ifdef __cplusplus
}
#endif

#endif // __gauss_h
