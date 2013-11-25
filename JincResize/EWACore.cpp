#include "EWACore.h"
#include <math.h>

const int LUT_SIZE = 65536;

void EWACore::InitLutTable()
{
  lut = new float[LUT_SIZE];

  float filter_end = GetSupport()*GetSupport();
  lut_factor = (float)LUT_SIZE / filter_end;

  for (int i = 0; i < LUT_SIZE; i++) {
    lut[i] = factor((float)i / lut_factor);
  }
}

void EWACore::DestroyLutTable()
{
  delete[] lut;
}

float EWACore::GetFactor(float dist)
{
  int index = int(dist*lut_factor);

  if (index >= LUT_SIZE)
      return 0.0;

  return lut[int(dist*lut_factor)];
}
