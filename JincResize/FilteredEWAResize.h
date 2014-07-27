
#ifndef __FILTERED_EWA_RESIZE_H
#define __FILTERED_EWA_RESIZE_H

#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

#include "EWACore.h"
#include "EWAResizerStruct.h"

typedef void(*EWAResizeCore)(EWAPixelCoeff* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                              int src_width, int src_height, int dst_width, int dst_height,
                              double crop_left, double crop_top, double crop_width, double crop_height);

class FilteredEWAResize : public GenericVideoFilter
{
public:
  FilteredEWAResize(PClip _child, int width, int height,
                    double crop_left, double crop_top, double crop_width, double crop_height,
                    EWACore *func, IScriptEnvironment* env);
  ~FilteredEWAResize();
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);

  static EWAResizeCore GetResizer(int filter_size, IScriptEnvironment* env);

private:
  EWACore *func;
  EWAResizeCore resizer;

  int src_width, src_height;
  double crop_left, crop_top, crop_width, crop_height;

  EWAPixelCoeff *stored_coeff_y, *stored_coeff_u, *stored_coeff_v;
};

#endif
