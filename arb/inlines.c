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
    Copyright (C) 2015 Tommy Hofmann

******************************************************************************/

#define ARB_INLINES_C
#include "arb.h"

mag_ptr arb_get_rad(arb_t x)
{
  mag_ptr y;
  y = &(x->rad);
  return y;
}

arf_ptr arb_get_mid(arb_t x)
{
  arf_ptr y;
  y = &(x->mid);
  return y;
}
