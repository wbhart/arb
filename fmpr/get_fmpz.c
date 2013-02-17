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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

void
fmpr_get_fmpz(fmpz_t z, const fmpr_t x, fmpr_rnd_t rnd)
{
    long exp;

    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            fmpz_zero(z);
        }
        else
        {
            printf("fmpr_get_fmpz: cannot convert infinity or nan to integer\n");
            abort();
        }
    }

    if (!COEFF_IS_MPZ(*fmpr_expref(x)))
    {
        exp = *fmpr_expref(x);
    }
    else
    {
        /* tiny */
        if (fmpz_sgn(fmpr_expref(x)) < 0)
        {
            int sign = fmpz_sgn(fmpr_manref(x));

            if (rnd == FMPR_RND_NEAR
                || rnd == FMPR_RND_DOWN
                || (rnd == FMPR_RND_FLOOR && sign > 0)
                || (rnd == FMPR_RND_CEIL && sign < 0))
                fmpz_zero(z);
            else
                fmpz_set_si(z, sign);
            return;
        }
        else
        {
            printf("fmpr_get_fmpz: number too large to convert to integer\n");
            abort();
        }
    }

    if (exp >= 0)
    {
        fmpz_mul_2exp(z, fmpr_manref(x), exp);
    }
    else
    {
        exp = -exp;

        if (rnd == FMPR_RND_NEAR)
        {
            int sign = fmpz_sgn(fmpr_manref(x));

            /* half-integer -> tie, so round to even */
            if (exp == 1)
            {
                fmpz_tdiv_q_2exp(z, fmpr_manref(x), 1);

                if (fmpz_is_odd(z))
                {
                    if (sign > 0)
                        fmpz_add_ui(z, z, 1);
                    else
                        fmpz_sub_ui(z, z, 1);
                }
            }
            else
            {
                /* add 1/2 and round to floor */
                fmpz_t t;
                fmpz_init(t);
                fmpz_one(t);
                fmpz_mul_2exp(t, t, exp - 1);
                fmpz_add(t, t, fmpr_manref(x));
                fmpz_fdiv_q_2exp(z, t, exp);
                fmpz_clear(t);
            }
        }
        else if (rnd == FMPR_RND_DOWN)
            fmpz_tdiv_q_2exp(z, fmpr_manref(x), exp);
        else if (rnd == FMPR_RND_UP)
            fmpz_adiv_q_2exp(z, fmpr_manref(x), exp);
        else if (rnd == FMPR_RND_FLOOR)
            fmpz_fdiv_q_2exp(z, fmpr_manref(x), exp);
        else if (rnd == FMPR_RND_CEIL)
            fmpz_cdiv_q_2exp(z, fmpr_manref(x), exp);
    }
}

