#include "acb_poly.h"

long _acb_poly_length(const acb_poly_t poly)
{
    return poly->length;
}
long _acb_poly_degree(const acb_poly_t poly)
{
    return poly->length - 1;
}

acb_ptr _acb_poly_arr_get(acb_ptr vec, long i)
{
  acb_ptr a = vec+i;
  return a;
}

int check_accuracy(acb_ptr vec, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (mag_cmp_2exp_si(arb_radref(acb_realref(vec + i)), -prec) >= 0
         || mag_cmp_2exp_si(arb_radref(acb_imagref(vec + i)), -prec) >= 0)
            return 0;
    }

    return 1;
}

acb_ptr
poly_roots(const fmpz_poly_t poly,
    long initial_prec,
    long target_prec);
{
    long i, prec, deg, isolated, maxiter;
    acb_poly_t cpoly;
    acb_ptr roots;

    deg = poly->length - 1;

    acb_poly_init(cpoly);
    roots = _acb_vec_init(deg);

    for (prec = initial_prec; ; prec *= 2)
    {
        acb_poly_set_fmpz_poly(cpoly, poly, prec);
        maxiter = FLINT_MIN(FLINT_MAX(deg, 32), prec);

        isolated = acb_poly_find_roots(roots, cpoly,
            prec == initial_prec ? NULL : roots, maxiter, prec);

        if (isolated == deg && check_accuracy(roots, deg, target_prec))
        {
            break;
        }
    }
    acb_poly_clear(cpoly);
    return roots;
}

