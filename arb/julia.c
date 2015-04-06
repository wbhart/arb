#include "arb.h"

mag_ptr
arb_get_rad(arb_t x)
{
  mag_ptr y;
  y = &(x->rad);
  return y;
}

arf_ptr
arb_get_mid(arb_t x)
{
  arf_ptr y;
  y = &(x->mid);
  return y;
}
