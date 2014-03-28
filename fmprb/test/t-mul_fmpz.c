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

    Copyright (C) 2012-2014 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

int main()
{
    long iter, iter2;
    flint_rand_t state;

    printf("mul_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* main test */
    for (iter = 0; iter < 1000; iter++)
    {
        fmprb_t x, y, z;
        fmpz_t c;
        long prec;

        fmprb_init(x);
        fmprb_init(y);
        fmprb_init(z);
        fmpz_init(c);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            fmprb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);
            fmprb_randtest_special(y, state, n_randint(state,2) ? 2000 : 200, 200);
            fmpz_randtest(c, state, 1 + n_randint(state, 2000));

            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 2))
            {
                fmprb_mul_fmpz(y, x, c, prec);
                fmprb_mul_fmpz_naive(z, x, c, prec);

                if (!fmprb_equal(y, z))
                {
                    printf("FAIL!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("y = "); fmprb_print(y); printf("\n\n");
                    printf("z = "); fmprb_print(z); printf("\n\n");
                    printf("c = "); fmpz_print(c); printf("\n\n");
                    abort();
                }
            }
            else
            {
                fmprb_mul_fmpz(y, x, c, prec);
                fmprb_mul_fmpz(x, x, c, prec);

                if (!fmprb_equal(x, y))
                {
                    printf("FAIL (aliasing)!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("y = "); fmprb_print(y); printf("\n\n");
                    printf("c = "); fmpz_print(c); printf("\n\n");
                    abort();
                }
                break;
            }
        }

        fmprb_clear(x);
        fmprb_clear(y);
        fmprb_clear(z);
        fmpz_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
