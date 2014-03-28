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
#include "ufloat.h"

void
fmprb_mul_ui_naive(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    fmprb_mul_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_mul_si_naive(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_si(t, y);
    fmprb_mul_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_mul_fmpz_naive(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    fmprb_mul_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

static __inline__ int
fmprb_is_special(const fmprb_t x)
{
    return fmpr_is_special(fmprb_midref(x)) || fmpr_is_special(fmprb_radref(x));
}

#define NORMALISE_LIMB(_exp,_y,_x) \
    do { \
        unsigned int __bc; \
        mp_limb_t __t; \
        __t = (_x); \
        count_trailing_zeros(__bc, __t); \
        (_y) = (__t) >> __bc; \
        (_exp) = __bc; \
    } while (0)

static __inline__ void
ufloat_add_r_get_error(fmpr_t rad, const fmpr_t mid, const ufloat_t err, long r)
{
    if (r == FMPR_RESULT_EXACT)
    {
        ufloat_get_fmpr(rad, err);
    }
    else
    {
        fmpz exp = *fmpr_expref(mid);

        if (exp < UFLOAT_MIN_EXP || exp > UFLOAT_MAX_EXP)
        {
            ufloat_get_fmpr(rad, err);
            fmpr_add_error_result(rad, rad, mid, r, FMPRB_RAD_PREC, FMPR_RND_UP);
        }
        else
        {
            ufloat_t err2;
            ufloat_add_2exp(err2, err, exp - r);
            ufloat_get_fmpr(rad, err2);
        }
    }
}

void
fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    ufloat_t zr, ar, yr;
    int xsign, ysign;
    mp_ptr xptr, yptr;
    mp_limb_t xtmp, ytmp;
    mp_size_t xn, yn;
    fmpz_t yexp;
    long r;

    fmprb_mul_fmpz_naive(z, x, y, prec);
    return;

    if (fmprb_is_exact(x))
    {
        r = fmpr_mul_fmpz(fmprb_midref(z),
            fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        return;
    }
    else if (fmpz_is_zero(y) && fmprb_is_finite(x))
    {
        fmprb_zero(z);
        return;
    }
    else if (fmpz_is_zero(y) || fmprb_is_special(x))
    {
        fmprb_mul_fmpz_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))
    FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, *y)
    *yexp = 0;

    if (!( ufloat_set_fmpr(ar, fmprb_radref(x))
        && ufloat_set_mpn_2exp_fmpz(yr, yptr, yn, yexp)))
    {
        fmprb_mul_fmpz_naive(z, x, y, prec);
        return;
    }

    ufloat_mul(zr, yr, ar);

    if (xn == 1 && yn == 1)
    {
        mp_limb_t t;
        NORMALISE_LIMB(*yexp, t, yptr[0]);
        r = _fmpr_mul_1x1(fmprb_midref(z),
            xptr[0], fmpr_expref(fmprb_midref(x)),
            t, yexp, xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else if (xn >= yn)
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            yptr, yn, yexp,
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            yptr, yn, yexp,
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            ysign ^ xsign, prec, FMPR_RND_DOWN);
    }

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}



void
fmprb_mul_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    ufloat_t zr, ar, yr;
    int xsign, ysign;
    mp_ptr xptr, yptr;
    mp_limb_t xtmp, ytmp;
    long xn, yn;
    long r;

    if (fmprb_is_exact(x))
    {
        r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        return;
    }
    else if (fmprb_is_special(x) || fmpr_is_special(y))
    {
        fmprb_mul_fmpr_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))
    FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, *fmpr_manref(y))

    if (!( ufloat_set_fmpr(ar, fmprb_radref(x))
        && ufloat_set_mpn_2exp_fmpz(yr, yptr, yn, fmpr_expref(y))))
    {
        fmprb_mul_fmpr_naive(z, x, y, prec);
        return;
    }

    ufloat_mul(zr, yr, ar);

    if (xn == 1 && yn == 1)
    {
        r = _fmpr_mul_1x1(fmprb_midref(z),
            xptr[0], fmpr_expref(fmprb_midref(x)),
            yptr[0], fmpr_expref(y),
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else if (xn >= yn)
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            yptr, yn, fmpr_expref(y),
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            yptr, yn, fmpr_expref(y),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            ysign ^ xsign, prec, FMPR_RND_DOWN);
    }

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}

void
fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    ufloat_t zr, xr, ar, yr, br;
    int xsign, ysign;
    mp_ptr xptr, yptr;
    mp_limb_t xtmp, ytmp;
    long xn, yn;
    long r;

    if (fmprb_is_exact(y))
    {
        fmprb_mul_fmpr(z, x, fmprb_midref(y), prec);
        return;
    }
    else if (fmprb_is_exact(x))
    {
        fmprb_mul_fmpr(z, y, fmprb_midref(x), prec);
        return;
    }
    else if (fmprb_is_special(x) || fmprb_is_special(y))
    {
        fmprb_mul_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))
    FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, *fmpr_manref(fmprb_midref(y)))

    if (!(ufloat_set_fmpr(ar, fmprb_radref(x))
        && ufloat_set_fmpr(br, fmprb_radref(y))
        && ufloat_set_mpn_2exp_fmpz(xr, xptr, xn, fmpr_expref(fmprb_midref(x)))
        && ufloat_set_mpn_2exp_fmpz(yr, yptr, yn, fmpr_expref(fmprb_midref(y)))))
    {
        fmprb_mul_naive(z, x, y, prec);
        return;
    }

    ufloat_mul(zr, xr, br);
    ufloat_addmul(zr, yr, ar);
    ufloat_addmul(zr, ar, br);

    if (xn == 1 && yn == 1)
    {
        r = _fmpr_mul_1x1(fmprb_midref(z),
            xptr[0], fmpr_expref(fmprb_midref(x)),
            yptr[0], fmpr_expref(fmprb_midref(y)),
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else if (xn >= yn)
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            yptr, yn, fmpr_expref(fmprb_midref(y)),
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            yptr, yn, fmpr_expref(fmprb_midref(y)),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            ysign ^ xsign, prec, FMPR_RND_DOWN);
    }

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}

void
fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    ufloat_t zr, ar, yr;
    int xsign, ysign;
    mp_ptr xptr, yptr;
    mp_limb_t xtmp, ytmp;
    mp_size_t xn, yn;
    fmpz_t yexp;
    long r;

    if (fmprb_is_exact(x))
    {
        r = fmpr_mul_ui(fmprb_midref(z),
            fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        return;
    }
    else if ((y == 0) && fmprb_is_finite(x))
    {
        fmprb_zero(z);
        return;
    }
    else if ((y == 0) || fmprb_is_special(x))
    {
        fmprb_mul_ui_naive(z, x, y, prec);
        return;
    }

    if (!ufloat_set_fmpr(ar, fmprb_radref(x)))
    {
        fmprb_mul_ui_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))

    ytmp = y;
    ysign = 0;
    yn = 1;
    yptr = &ytmp;
    NORMALISE_LIMB(*yexp, ytmp, ytmp);

    ufloat_set_ui(yr, ytmp); yr->e += *yexp;
    ufloat_mul(zr, yr, ar);

    if (xn == 1)
    {
        r = _fmpr_mul_1x1(fmprb_midref(z),
            xptr[0], fmpr_expref(fmprb_midref(x)),
            yptr[0], yexp, xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            yptr, yn, yexp,
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}

void
fmprb_mul_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    ufloat_t zr, ar, yr;
    int xsign, ysign;
    mp_ptr xptr, yptr;
    mp_limb_t xtmp, ytmp;
    mp_size_t xn, yn;
    fmpz_t yexp;
    long r;

    if (fmprb_is_exact(x))
    {
        r = fmpr_mul_si(fmprb_midref(z),
            fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        return;
    }
    else if ((y == 0) && fmprb_is_finite(x))
    {
        fmprb_zero(z);
        return;
    }
    else if ((y == 0) || fmprb_is_special(x))
    {
        fmprb_mul_si_naive(z, x, y, prec);
        return;
    }

    if (!ufloat_set_fmpr(ar, fmprb_radref(x)))
    {
        fmprb_mul_si_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))

    if (y > 0)
    {
        ytmp = y;
        ysign = 0;
    }
    else
    {
        ytmp = -y;
        ysign = 1;
    }

    NORMALISE_LIMB(*yexp, ytmp, ytmp);
    yn = 1;
    yptr = &ytmp;

    ufloat_set_ui(yr, ytmp); yr->e += *yexp;
    ufloat_mul(zr, yr, ar);

    if (xn == 1)
    {
        r = _fmpr_mul_1x1(fmprb_midref(z),
            xptr[0], fmpr_expref(fmprb_midref(x)),
            yptr[0], yexp, xsign ^ ysign, prec, FMPR_RND_DOWN);
    }
    else
    {
        r = _fmpr_mul_mpn(fmprb_midref(z),
            xptr, xn, fmpr_expref(fmprb_midref(x)),
            yptr, yn, yexp,
            xsign ^ ysign, prec, FMPR_RND_DOWN);
    }

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}

