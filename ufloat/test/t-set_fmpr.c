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

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_fmpr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t x, y, z;
        ufloat_t u;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);

        fmpr_randtest_special(x, state, 2 + n_randint(state, 1000), 2 + n_randint(state, 100));
        fmpr_randtest_special(y, state, 2 + n_randint(state, 1000), 2 + n_randint(state, 100));

        if (ufloat_set_fmpr(u, x))
        {
            ufloat_get_fmpr(y, u);
            fmpr_mul_2exp_si(z, x, 1);

            /* check that |x| <= u <= 2|x| */
            if (!(fmpr_cmpabs(x, y) <= 0 && fmpr_cmpabs(y, z) <= 0))
            {
                printf("FAIL\n\n");
                printf("x = "); fmpr_print(x); printf("\n\n");
                printf("y = "); fmpr_print(y); printf("\n\n");
                printf("z = "); fmpr_print(z); printf("\n\n");
                abort();
            }

        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
