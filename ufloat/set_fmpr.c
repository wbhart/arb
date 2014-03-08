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
            t = FLINT_ABS(v);
            return ufloat_set_mpn_2exp_fmpz(rr, &t, 1, fmpr_expref(x));
        }
        else
        {
            __mpz_struct * z = COEFF_TO_PTR(v);
            return ufloat_set_mpn_2exp_fmpz(rr, z->_mp_d, FLINT_ABS(z->_mp_size), fmpr_expref(x));
        }
    }
}

