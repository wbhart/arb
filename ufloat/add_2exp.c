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

******************************************************************************/

#include "ufloat.h"

void
ufloat_add_2exp(ufloat_t z, const ufloat_t x, long e)
{
    if (x->m == 0)
    {
        z->m = UFLOAT_ONE_HALF;
        z->e = e + 1;
    }
    else
    {
        long shift;
        shift = x->e - e;

        if (shift >= 0)
        {
            z->e = x->e;
            if (shift >= UFLOAT_BITS)
                z->m = x->m + 1;
            else
                z->m = x->m + (1UL << (UFLOAT_BITS - shift));
        }
        else
        {
            shift = -shift;
            z->e = e;
            if (shift >= UFLOAT_BITS)
                z->m = (1UL << UFLOAT_BITS) + 1;
            else
                z->m = (1UL << UFLOAT_BITS) + (x->m >> shift);
        }

        UFLOAT_ADJUST_ONE_TOO_LARGE(z->m, z->e)
    }
}

