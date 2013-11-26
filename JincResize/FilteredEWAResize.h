
#ifndef __FILTERED_EWA_RESIZE_H
#define __FILTERED_EWA_RESIZE_H

#include "avisynth.h"
#include "EWACore.h"


typedef void (*EWAResizeCore)(EWACore* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
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

private:
  EWACore *func;
  EWAResizeCore resizer;

  int src_width, src_height;
  double crop_left, crop_top, crop_width, crop_height;
};

#endif
