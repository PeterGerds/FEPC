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

#ifndef __CONFIG_H
#define __CONFIG_H
/**
 **
 ** provides basic types and configuration settings
 **
 ** arch-tag: 3c370550-398d-496c-ac9e-b4b7cb79c45f
 **/

/*
 * included header
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "seconds.h"

/*
 * basic types
 */

/*
 * Adds c++ compatibility
 */
#if !defined(__cplusplus)
	typedef enum { false, true } bool_t;
#else
	typedef bool bool_t;
#endif

typedef double  fepc_real_t;

/*
 * consistency check
 */

#if !defined(NDEBUG)
#  define ASSERT( cond )   assert( cond )
#else
#  define ASSERT( cond )
#endif

/*
 * enables debugging/output
 */

#define DEBUG  if ( false )
#define PRINT  if ( false )

/*
 * inline function decl.
 */

#if defined(__USE_ISOC99) && ( defined(__GNUC__) || defined(__ICC) || defined(__ECC) )
#define _INLINE_  inline
#else
#define _INLINE_  static
#endif

#endif  /* __CONFIG_H */
