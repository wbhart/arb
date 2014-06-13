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

#include <pthread.h>
#include "acb_poly.h"

typedef struct
{
    acb_ptr z;
    acb_srcptr s;
    acb_srcptr a;
    long n0;
    long n1;
    long d0;
    long len;
    long prec;
}
powsum_arg_t;

void *
_acb_zeta_powsum_evaluator(void * arg_ptr)
{
    powsum_arg_t arg = *((powsum_arg_t *) arg_ptr);
    long i, k;

    acb_t t, u, v;
    arb_t f;

    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_init(f);

    _acb_vec_zero(arg.z, arg.len);

    for (k = arg.n0; k < arg.n1; k++)
    {
        /* t = log(a+k) */
        acb_add_ui(t, arg.a, k, arg.prec);
        acb_log(t, t, arg.prec);

        /* u = (a+k)^(-s) */
        acb_mul(u, t, arg.s, arg.prec);
        acb_neg(u, u);
        acb_exp(u, u, arg.prec);

        /* forward: u *= (-1)^d * log(a+k)^d / d! */
        if (arg.d0 != 0)
        {
            acb_pow_ui(v, t, arg.d0, arg.prec);
            acb_mul(u, u, v, arg.prec);
            arb_fac_ui(f, arg.d0, arg.prec);
            arb_div(acb_realref(u), acb_realref(u), f, arg.prec);
            arb_div(acb_imagref(u), acb_imagref(u), f, arg.prec);
            if (arg.d0 % 2)
                acb_neg(u, u);
        }

        acb_add(arg.z, arg.z, u, arg.prec);

        for (i = 1; i < arg.len; i++)
        {
            acb_mul(u, u, t, arg.prec);
            acb_div_si(u, u, -(arg.d0 + i), arg.prec);
            acb_add(arg.z + i, arg.z + i, u, arg.prec);
        }
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    arb_clear(f);

    flint_cleanup();

    return NULL;
}

void
_acb_poly_powsum_series_naive_threaded(acb_ptr z,
    const acb_t s, const acb_t a, long n, long len, long prec)
{
    pthread_t * threads;
    powsum_arg_t * args;
    long i, num_threads;
    int split_each_term;

    num_threads = flint_get_num_threads();

    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(powsum_arg_t) * num_threads);

    split_each_term = (len > 1000);

    for (i = 0; i < num_threads; i++)
    {
        args[i].s = s;
        args[i].a = a;

        if (split_each_term)
        {
            long n0, n1;
            n0 = (len * i) / num_threads;
            n1 = (len * (i + 1)) / num_threads;
            args[i].z = z + n0;
            args[i].n0 = 0;
            args[i].n1 = n;
            args[i].d0 = n0;
            args[i].len = n1 - n0;
        }
        else
        {
            args[i].z = _acb_vec_init(len);
            args[i].n0 = (n * i) / num_threads;
            args[i].n1 = (n * (i + 1)) / num_threads;
            args[i].d0 = 0;
            args[i].len = len;
        }

        args[i].prec = prec;
        pthread_create(&threads[i], NULL, _acb_zeta_powsum_evaluator, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    if (!split_each_term)
    {
        _acb_vec_zero(z, len);
        for (i = 0; i < num_threads; i++)
        {
            _acb_vec_add(z, z, args[i].z, len, prec);
            _acb_vec_clear(args[i].z, len);
        }
    }

    flint_free(threads);
    flint_free(args);
}

