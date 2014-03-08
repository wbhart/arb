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

/* assumes xn >= 1, top limb of xp set */
int
ufloat_set_mpn_2exp_fmpz(ufloat_t rr, mp_srcptr xp,
    mp_size_t xn, const fmpz_t exp)
{
    mp_limb_t m;
    long e, bits, shift;

    if (*exp < UFLOAT_MIN_EXP || *exp > UFLOAT_MAX_EXP)
        return 0;

    m = xp[xn - 1];

    count_leading_zeros(bits, m);
    bits = FLINT_BITS - bits;

    e = *exp + (xn - 1) * FLINT_BITS + bits;
    shift = bits - UFLOAT_BITS;

    if (shift >= 0)
        m = (m >> shift) + 1;
    else if (xn == 1)
        m = m << (-shift);
    else
        m = ((m << (-shift)) | (xp[xn-2] >> (FLINT_BITS - (-shift)))) + 1;

    UFLOAT_ADJUST_ONE_TOO_LARGE(m, e)

    if (e < UFLOAT_MIN_EXP || e > UFLOAT_MAX_EXP)
        return 0;

    rr->m = m;
    rr->e = e;

    return 1;
}

