#include "acb.h"

arb_ptr
acb_get_real(acb_t x)
{
  arb_ptr y;
  y = &(x->real);
  return y;
}

arb_ptr
acb_get_imag(acb_t x)
{
  arb_ptr y;
  y = &(x->imag);
  return y;
}
