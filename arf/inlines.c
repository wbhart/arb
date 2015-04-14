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

#define ARF_INLINES_C
#include "arf.h"

arf_ptr arf_vec_init(long n)
{
    long i;
    arf_ptr v = (arf_ptr) flint_malloc(sizeof(arf_struct) * n);

    for (i = 0; i < n; i++)
        arf_init(v + i);

    return v;
}

void arf_vec_clear(arf_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arf_clear(v + i);
    flint_free(v);
}

void _arf_mul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
{
    ((rnd == FMPR_RND_DOWN) ? arf_mul_rnd_down(z, x, y, prec) : arf_mul_rnd_any(z, x, y, prec, rnd));
}
