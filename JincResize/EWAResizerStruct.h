#pragma once

#ifndef __EWARESIZESTRUCT_H
#define __EWARESIZESTRUCT_H

struct EWAPixelCoeffMeta
{
  int start_x, start_y;
  int coeff;
};

struct EWAPixelCoeff
{
  float* factor;
  EWAPixelCoeffMeta* meta;
  int* factor_map;

  int filter_size, quantize_x, quantize_y, coeff_stripe;
};

#endif 

