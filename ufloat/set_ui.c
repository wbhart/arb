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

int
ufloat_set_ui(ufloat_t rr, ulong x)
{
    if (x == 0)
    {
        ufloat_zero(rr);
    }
    else
    {
        long bits, e, shift;
        mp_limb_t t = x;

        count_leading_zeros(bits, t);
        bits = FLINT_BITS - bits;
        e = bits;
        shift = bits - UFLOAT_BITS;

        if (shift <= 0)
            t = t << (-shift);
        else
            t = (t >> shift) + 1;

        UFLOAT_ADJUST_ONE_TOO_LARGE(t, e)

        rr->m = t;
        rr->e = e;
    }

    return 1;
}

