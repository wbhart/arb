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

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("submul_arf....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c, d;
        arf_t x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arf_init(x);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        arf_randtest_special(x, state, 1 + n_randint(state, 2000), 100);

        prec = 2 + n_randint(state, 2000);

        arb_set_arf(b, x);
        arb_set(d, c);
        arb_submul_arf(c, a, x, prec);
        arb_submul(d, a, b, prec);

        if (!arb_equal(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arf_clear(x);
    }

    /* aliasing */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b, c;
        arf_t x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arf_init(x);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        arf_randtest_special(x, state, 1 + n_randint(state, 2000), 100);

        prec = 2 + n_randint(state, 2000);

        arb_set_arf(b, x);
        arb_set(c, a);
        arb_submul_arf(c, a, x, prec);
        arb_submul_arf(a, a, x, prec);

        if (!arb_equal(a, c))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arf_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
