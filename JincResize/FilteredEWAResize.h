
#ifndef __FILTERED_EWA_RESIZE_H
#define __FILTERED_EWA_RESIZE_H

#include "avisynth.h"
#include "EWACore.h"

class FilteredEWAResize : public GenericVideoFilter
{
public:
  FilteredEWAResize(PClip _child, int width, int height, EWACore *func, IScriptEnvironment* env);
  ~FilteredEWAResize();

  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);

private:
  void ResizePlane(BYTE* dstp, const BYTE* srcp, int dst_pitch, int src_pitch);

  EWACore *func;
  double src_width, src_height;
};

#endif
