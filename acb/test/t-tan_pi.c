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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("tan_pi....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t x, y, a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(y);
        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        acb_tan_pi(a, x, prec1);
        acb_tan_pi(b, x, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            abort();
        }

        /* compare with tan */
        arb_const_pi(acb_realref(c), prec1);
        acb_mul_arb(y, x, acb_realref(c), prec1);
        acb_tan(c, y, prec1);

        if (!acb_overlaps(a, c))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("y = "); acb_print(y); flint_printf("\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            abort();
        }

        acb_tan_pi(x, x, prec1);

        if (!acb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

