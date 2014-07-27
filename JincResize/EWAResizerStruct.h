#pragma once

#ifndef __EWARESIZESTRUCT_H
#define __EWARESIZESTRUCT_H

struct EWAPixelCoeffMeta
{
  int start_x, start_y;
  float* coeff;
};

struct EWAPixelCoeff
{
  float* factor;
  EWAPixelCoeffMeta* meta;
  float** factor_map;

  int filter_size, quantize_x, quantize_y, coeff_stripe;
};

#endif 

