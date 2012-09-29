#include <stdlib.h>
#include "Core.hh"

Vector operator*(ft l, const Vector& r)
{
  return Vector(r.x * l, r.y * l);
}

ft getRand(ft m)
{
  return rand() * m / (RAND_MAX + 1.0);
}
