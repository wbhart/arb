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
ufloat_set_fmpr(ufloat_t rr, const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            ufloat_zero(rr);
            return 1;
        }
        else
            return 0;
    }
    else
    {
        mp_limb_t t;
        fmpz v = *fmpr_manref(x);

        if (!COEFF_IS_MPZ(v))
        {
            long e, bits, shift;
            fmpz exp = *fmpr_expref(x);

            t = FLINT_ABS(v);

            if (exp < UFLOAT_MIN_EXP || exp > UFLOAT_MAX_EXP)
                return 0;

            count_leading_zeros(bits, t);
            bits = FLINT_BITS - bits;
            e = exp + bits;
            shift = bits - UFLOAT_BITS;

            if (shift <= 0)
                t = t << (-shift);
            else
                t = (t >> shift) + 1;

            UFLOAT_ADJUST_ONE_TOO_LARGE(t, e)

            if (e < UFLOAT_MIN_EXP || e > UFLOAT_MAX_EXP)
                return 0;

            rr->m = t;
            rr->e = e;

            return 1;
        }
        else
        {
            __mpz_struct * z = COEFF_TO_PTR(v);
            return ufloat_set_mpn_2exp_fmpz(rr, z->_mp_d, FLINT_ABS(z->_mp_size), fmpr_expref(x));
        }
    }
}

