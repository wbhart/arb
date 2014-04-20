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

#include "fmprb.h"

void
fmprb_addmul_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_ui(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_si(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_fmpz(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_fmpr(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

#include "ufloat.h"


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

#define FMPR_ADD_MPN(_r,_res,_tptr,_tn,_tsign,_texp,_uptr,_un,_usign,_uexp,_prec)      \
  do {                                                                      \
    long __shift = _fmpz_sub_small(_uexp, _texp);                             \
    if (__shift >= 0)                                                       \
    {                                                                       \
        if ((_tn == 1) && (_un == 1) && (__shift < FLINT_BITS))               \
            _r = _fmpr_add_1x1(_res, _tptr[0], _tsign, _texp,                    \
                _uptr[0], _usign, _uexp, __shift, _prec, FMPR_RND_DOWN);        \
        else                                                                \
            _r = _fmpr_add_mpn(_res, _tptr, _tn, _tsign, _texp,                   \
                _uptr, _un, _usign, _uexp, __shift, _prec, FMPR_RND_DOWN);       \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        __shift = -__shift;                                                     \
        if ((_tn == 1) && (_un == 1) && (__shift < FLINT_BITS))                 \
            _r = _fmpr_add_1x1(_res, _uptr[0], _usign, _uexp,                    \
                _tptr[0], _tsign, _texp, __shift, _prec, FMPR_RND_DOWN);        \
        else                                                                \
            _r = _fmpr_add_mpn(_res, _uptr, _un, _usign, _uexp,                   \
                _tptr, _tn, _tsign, _texp, __shift, _prec, FMPR_RND_DOWN);       \
    } \
  } while (0);

#define FMPR_MPN_MUL(_cptr,_cn,_aptr,_an,_bptr,_bn)                         \
  do {                                                                      \
    mp_limb_t __cy;                                                         \
    if (_an >= _bn)                                                         \
    {                                                                       \
        if (_bn == 1)                                                       \
            _cptr[_an] = __cy = mpn_mul_1(_cptr, _aptr, _an, _bptr[0]);     \
        else                                                                \
            __cy = mpn_mul(_cptr, _aptr, _an, _bptr, _bn);                  \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        if (_an == 1)                                                       \
            _cptr[_bn] = __cy = mpn_mul_1(_cptr, _bptr, _bn, _aptr[0]);     \
        else                                                                \
            __cy = mpn_mul(_cptr, _bptr, _bn, _aptr, _an);                  \
    }                                                                       \
    _cn = _an + _bn - (__cy == 0);                                          \
  } while (0);

extern TLS_PREFIX mp_ptr __mul_tmp;
extern TLS_PREFIX long __mul_alloc;

void _mul_tmp_cleanup(void);

#define MUL_STACK_ALLOC 40
#define MUL_TLS_ALLOC 1000

#define MUL_TMP_ALLOC(tmp, alloc) \
    if (alloc <= MUL_STACK_ALLOC) \
    { \
        tmp = tmp_stack; \
    } \
    else if (alloc <= MUL_TLS_ALLOC) \
    { \
        if (__mul_alloc < alloc) \
        { \
            if (__mul_alloc == 0) \
            { \
                flint_register_cleanup_function(_mul_tmp_cleanup); \
            } \
            __mul_tmp = flint_realloc(__mul_tmp, sizeof(mp_limb_t) * alloc); \
            __mul_alloc = alloc; \
        } \
        tmp = __mul_tmp; \
    } \
    else \
    { \
        tmp = flint_malloc(sizeof(mp_limb_t) * alloc); \
    }

#define MUL_TMP_FREE(tmp, alloc) \
    if (alloc > MUL_TLS_ALLOC) \
        flint_free(tmp);


void
fmprb_addmul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    ufloat_t xr, ar, yr, br, zr;
    int xsign, ysign, tsign, zsign;
    mp_ptr xptr, yptr, tptr, zptr;
    mp_limb_t xtmp, ytmp, ztmp;
    mp_limb_t tmp_stack[MUL_STACK_ALLOC];
    long xn, yn, zn, tn, alloc, r;
    fmpz_t texp;

    if (fmprb_is_exact(y))
    {
        fmprb_addmul_fmpr(z, x, fmprb_midref(y), prec);
        return;
    }
    else if (fmprb_is_exact(x))
    {
        fmprb_addmul_fmpr(z, y, fmprb_midref(x), prec);
        return;
    }
    else if (fmprb_is_zero(z))
    {
        fmprb_mul(z, x, y, prec);
        return;
    }
    else if (fmprb_is_special(x) || fmprb_is_special(y) || fmprb_is_special(z))
    {
        fmprb_addmul_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))
    FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, *fmpr_manref(fmprb_midref(y)))
    FMPZ_GET_MPN_READONLY(zsign, zn, zptr, ztmp, *fmpr_manref(fmprb_midref(z)))

    if (!( ufloat_set_fmpr(ar, fmprb_radref(x))
        && ufloat_set_fmpr(br, fmprb_radref(y))
        && ufloat_set_fmpr(zr, fmprb_radref(z))
        && ufloat_set_mpn_2exp_fmpz(xr, xptr, xn, fmpr_expref(fmprb_midref(x)))
        && ufloat_set_mpn_2exp_fmpz(yr, yptr, yn, fmpr_expref(fmprb_midref(y)))))
    {
        fmprb_addmul_naive(z, x, y, prec);
        return;
    }

    ufloat_addmul(zr, xr, br);
    ufloat_addmul(zr, yr, ar);
    ufloat_addmul(zr, ar, br);

    fmpz_init(texp);

    alloc = xn + yn;

    MUL_TMP_ALLOC(tptr, alloc)

    fmpz_add_inline(texp, fmpr_expref(fmprb_midref(x)), fmpr_expref(fmprb_midref(y)));
    tsign = xsign ^ ysign;
    FMPR_MPN_MUL(tptr, tn, xptr, xn, yptr, yn)
    FMPR_ADD_MPN(r, fmprb_midref(z), zptr, zn, zsign, fmpr_expref(fmprb_midref(z)), tptr, tn, tsign, texp, prec)

    MUL_TMP_FREE(tptr, alloc)

    fmpz_clear(texp);

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}

void
fmprb_submul_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_fmpr(t, x, y, prec);
    fmprb_sub(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_submul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    ufloat_t xr, ar, yr, br, zr;
    int xsign, ysign, tsign, zsign;
    mp_ptr xptr, yptr, tptr, zptr;
    mp_limb_t xtmp, ytmp, ztmp;
    mp_limb_t tmp_stack[MUL_STACK_ALLOC];
    long xn, yn, zn, tn, alloc, r;
    fmpz_t texp;

    if (fmprb_is_exact(y))
    {
        fmprb_submul_fmpr(z, x, fmprb_midref(y), prec);
        return;
    }
    else if (fmprb_is_exact(x))
    {
        fmprb_submul_fmpr(z, y, fmprb_midref(x), prec);
        return;
    }
    else if (fmprb_is_zero(z))
    {
        fmprb_mul(z, x, y, prec);
        fmprb_neg(z, z);
        return;
    }
    else if (fmpr_is_special(fmprb_midref(x)) || fmpr_is_special(fmprb_radref(x))
          || fmpr_is_special(fmprb_midref(y)) || fmpr_is_special(fmprb_radref(y))
          || fmpr_is_special(fmprb_midref(z)) || fmpr_is_special(fmprb_radref(z)))
    {
        fmprb_submul_naive(z, x, y, prec);
        return;
    }

    FMPZ_GET_MPN_READONLY(xsign, xn, xptr, xtmp, *fmpr_manref(fmprb_midref(x)))
    FMPZ_GET_MPN_READONLY(ysign, yn, yptr, ytmp, *fmpr_manref(fmprb_midref(y)))
    FMPZ_GET_MPN_READONLY(zsign, zn, zptr, ztmp, *fmpr_manref(fmprb_midref(z)))

    if (!( ufloat_set_fmpr(ar, fmprb_radref(x))
        && ufloat_set_fmpr(br, fmprb_radref(y))
        && ufloat_set_fmpr(zr, fmprb_radref(z))
        && ufloat_set_mpn_2exp_fmpz(xr, xptr, xn, fmpr_expref(fmprb_midref(x)))
        && ufloat_set_mpn_2exp_fmpz(yr, yptr, yn, fmpr_expref(fmprb_midref(y)))))
    {
        fmprb_submul_naive(z, x, y, prec);
        return;
    }

    ufloat_addmul(zr, xr, br);
    ufloat_addmul(zr, yr, ar);
    ufloat_addmul(zr, ar, br);

    fmpz_init(texp);

    alloc = xn + yn;

    MUL_TMP_ALLOC(tptr, alloc)

    fmpz_add_inline(texp, fmpr_expref(fmprb_midref(x)), fmpr_expref(fmprb_midref(y)));
    tsign = xsign ^ ysign ^ 1;
    FMPR_MPN_MUL(tptr, tn, xptr, xn, yptr, yn)
    FMPR_ADD_MPN(r, fmprb_midref(z), zptr, zn, zsign, fmpr_expref(fmprb_midref(z)), tptr, tn, tsign, texp, prec)

    MUL_TMP_FREE(tptr, alloc)

    fmpz_clear(texp);

    ufloat_add_r_get_error(fmprb_radref(z), fmprb_midref(z), zr, r);
}
