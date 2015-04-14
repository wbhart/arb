/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2015 Tommy Hofmann

******************************************************************************/

#define ACB_POLY_INLINES_C
#include "acb_poly.h"

long _acb_poly_length(const acb_poly_t poly)
{
    return poly->length;
}
long _acb_poly_degree(const acb_poly_t poly)
{
    return poly->length - 1;
}

acb_ptr _acb_poly_arr_get(acb_ptr vec, long i)
{
  acb_ptr a = vec+i;
  return a;
}

/* Following two functions are copied from the example */

int check_accuracy(acb_ptr vec, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (mag_cmp_2exp_si(arb_radref(acb_realref(vec + i)), -prec) >= 0
         || mag_cmp_2exp_si(arb_radref(acb_imagref(vec + i)), -prec) >= 0)
            return 0;
    }

    return 1;
}

acb_ptr
poly_roots(const fmpz_poly_t poly,
    long initial_prec,
    long target_prec)
{
    long i, prec, deg, isolated, maxiter;
    acb_poly_t cpoly;
    acb_ptr roots;

    deg = poly->length - 1;

    acb_poly_init(cpoly);
    roots = _acb_vec_init(deg);

    for (prec = initial_prec; ; prec *= 2)
    {
        acb_poly_set_fmpz_poly(cpoly, poly, prec);
        maxiter = FLINT_MIN(FLINT_MAX(deg, 32), prec);

        isolated = acb_poly_find_roots(roots, cpoly,
            prec == initial_prec ? NULL : roots, maxiter, prec);

        if (isolated == deg && check_accuracy(roots, deg, target_prec))
        {
            break;
        }
    }
    acb_poly_clear(cpoly);
    return roots;
}

