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

#ifndef UFLOAT_H
#define UFLOAT_H

#include "flint.h"
#include "fmpr.h"

#ifdef __cplusplus
extern "C" {
#endif

#define UFLOAT_BITS 30
#define UFLOAT_ONE_HALF (1UL << (UFLOAT_BITS - 1))
#define UFLOAT_MAX_EXP (1L << (FLINT_BITS-6))
#define UFLOAT_MIN_EXP (-UFLOAT_MAX_EXP)

typedef struct
{
    mp_limb_t m;
    long e;
}
ufloat_struct;

typedef ufloat_struct ufloat_t[1];

static __inline__ void
ufloat_zero(ufloat_t x)
{
    x->m = x->e = 0;
}

static __inline__ void
ufloat_get_fmpr(fmpr_t x, const ufloat_t r)
{
    fmpr_set_ui_2exp_si(x, r->m, r->e - UFLOAT_BITS);
}

static __inline__ void
ufloat_print(const ufloat_t x)
{
    fmpr_t t;
    fmpr_init(t);
    ufloat_get_fmpr(t, x);
    fmpr_printd(t, 15);
    fmpr_clear(t);
}

static __inline__ mp_limb_t
__ufloat_fixmul32(mp_limb_t x, mp_limb_t y)
{
    mp_limb_t u, v;
    umul_ppmm(u, v, x, y);
    return (u << (32 - UFLOAT_BITS)) | (v >> UFLOAT_BITS);
}

#if FLINT_BITS == 64
#define UFLOAT_FIXMUL(x, y) (((x) * (y)) >> UFLOAT_BITS)
#else
#define UFLOAT_FIXMUL(x, y) __ufloat_fixmul32((x), (y))
#endif

static __inline__ void
ufloat_mul(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    if (x->m == 0 || y->m == 0)
    {
        ufloat_zero(z);
    }
    else
    {
        z->m = UFLOAT_FIXMUL(x->m, y->m) + 1;
        z->e = x->e + y->e;
    }
}

#define UFLOAT_ADJUST_ONE_TOO_LARGE(m, e) \
    { \
        mp_limb_t __t = (m) >> UFLOAT_BITS; \
        (m) = ((m) >> __t) + __t; \
        (e) += __t; \
    }

#define UFLOAT_CHECK_BITS(rr) \
    if (FLINT_BIT_COUNT(rr->m) > UFLOAT_BITS) \
    { \
        printf("FAIL: overflow in mantissa!\n"); \
        abort(); \
    }

int ufloat_set_mpn_2exp_fmpz(ufloat_t rr, mp_srcptr xp, mp_size_t xn, const fmpz_t exp);

int ufloat_set_fmpr(ufloat_t rr, const fmpr_t x);

int ufloat_set_ui(ufloat_t rr, ulong x);

void ufloat_addmul(ufloat_t z, const ufloat_t x, const ufloat_t y);

void ufloat_add_2exp(ufloat_t z, const ufloat_t x, long e);

#ifdef __cplusplus
}
#endif

#endif

