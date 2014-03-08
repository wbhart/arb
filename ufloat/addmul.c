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
ufloat_addmul(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    if (z->m == 0)
    {
        ufloat_mul(z, x, y);
    }
    else if (x->m == 0 || y->m == 0)
    {
        return;
    }
    else
    {
        long shift, e;

        /* x*y < 2^e */
        e = x->e + y->e;
        shift = z->e - e;

        if (shift >= 0)
        {
            if (shift >= UFLOAT_BITS)
                z->m++;
            else
                z->m = z->m + (UFLOAT_FIXMUL(x->m, y->m) >> shift) + 1;
        }
        else
        {
            shift = -shift;
            z->e = e;

            if (shift >= UFLOAT_BITS)
                z->m = UFLOAT_FIXMUL(x->m, y->m) + 2;
            else
                z->m = UFLOAT_FIXMUL(x->m, y->m) + (z->m >> shift) + 2;
        }

        UFLOAT_ADJUST_ONE_TOO_LARGE(z->m, z->e)
    }
}

