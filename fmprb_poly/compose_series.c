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

#include "fmprb_poly.h"

void
_fmprb_poly_compose_series(fmprb_ptr res, fmprb_srcptr poly1, long len1,
                            fmprb_srcptr poly2, long len2, long n, long prec)
{
    if (len2 == 2)  /* fast linear case (TODO: all monomials) */
    {
        long i;
        fmprb_t t;

        fmprb_init(t);

        fmprb_set(t, poly2 + 1);
        fmprb_set_round(res, poly1, prec);

        for (i = 1; i < n; i++)
        {
            fmprb_mul(res + i, poly1 + i, t, prec);
            if (i + 1 < n)
                fmprb_mul(t, t, poly2 + 1, prec);
        }

        fmprb_clear(t);
    }
    else if (len1 < 6 || n < 6)
    {
        _fmprb_poly_compose_series_horner(res, poly1, len1, poly2, len2, n, prec);
    }
    else
    {
        _fmprb_poly_compose_series_brent_kung(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
fmprb_poly_compose_series(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmprb_is_zero(poly2->coeffs))
    {
        printf("exception: compose_series: inner "
                "polynomial must have zero constant term\n");
        abort();
    }

    if (len1 == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmprb_poly_set_fmprb(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        fmprb_poly_fit_length(res, lenr);
        _fmprb_poly_compose_series(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _fmprb_poly_set_length(res, lenr);
        _fmprb_poly_normalise(res);
    }
    else
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, lenr);
        _fmprb_poly_compose_series(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _fmprb_poly_set_length(t, lenr);
        _fmprb_poly_normalise(t);
        fmprb_poly_swap(res, t);
        fmprb_poly_clear(t);
    }
}
