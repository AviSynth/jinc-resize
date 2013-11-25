
#include <math.h>
#include "FilteredEWAResize.h"

FilteredEWAResize::FilteredEWAResize(PClip _child, int width, int height, double crop_left, double crop_top, double crop_width, double crop_height, EWACore *func, IScriptEnvironment* env) :
  GenericVideoFilter(_child),
  func( func ),
  crop_left(crop_left), crop_top(crop_top), crop_width(crop_width), crop_height(crop_height)
{
  if (!vi.IsPlanar() || !vi.IsYUV()) {
    env->ThrowError("JincResize: Only planar YUV colorspaces are supported");
  }

  if (width < vi.width || height < vi.height) {
    env->ThrowError("JincResize: java.lang.NotImplementedException");
  }

  if (vi.width < int(ceil(2*func->GetSupport())) || vi.height < int(ceil(2*func->GetSupport()))) {
    env->ThrowError("JincResize: Source image too small.");
  }

  src_width = vi.width;
  src_height = vi.height;

  vi.width = width;
  vi.height = height;
}

FilteredEWAResize::~FilteredEWAResize()
{
  delete func;
}

PVideoFrame __stdcall FilteredEWAResize::GetFrame(int n, IScriptEnvironment* env)
{
  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = env->NewVideoFrame(vi);

  // Luma
  ResizePlane(
    dst->GetWritePtr(), src->GetReadPtr(), dst->GetPitch(), src->GetPitch(),
    src_width, src_height, vi.width, vi.height,
    crop_left, crop_top, crop_width, crop_height
  );

  if (!vi.IsY8()) {
    int subsample_w = vi.GetPlaneWidthSubsampling(PLANAR_U);
    int subsample_h = vi.GetPlaneHeightSubsampling(PLANAR_U);

    double div_w = 1 << subsample_w;
    double div_h = 1 << subsample_h;

    ResizePlane(
      dst->GetWritePtr(PLANAR_U), src->GetReadPtr(PLANAR_U), dst->GetPitch(PLANAR_U), src->GetPitch(PLANAR_U),
      src_width >> subsample_w, src_height >> subsample_h, vi.width >> subsample_w, vi.height >> subsample_h,
      crop_left / div_w, crop_top / div_h, crop_width / div_w, crop_height / div_h
    );

    ResizePlane(
      dst->GetWritePtr(PLANAR_V), src->GetReadPtr(PLANAR_V), dst->GetPitch(PLANAR_V), src->GetPitch(PLANAR_V),
      src_width >> subsample_w, src_height >> subsample_h, vi.width >> subsample_w, vi.height >> subsample_h,
      crop_left / div_w, crop_top / div_h, crop_width / div_w, crop_height / div_h
    );
    
  }

  return dst;
}

void FilteredEWAResize::ResizePlane(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                                    int src_width, int src_height, int dst_width, int dst_height,
                                    double crop_left, double crop_top, double crop_width, double crop_height)
{
  double filter_support = func->GetSupport();
  int filter_size = ceil(filter_support * 2.0);

  double start_x = crop_left + (crop_width - vi.width) / (vi.width*2);
  double start_y = crop_top + (crop_height - vi.height) / (vi.height*2);

  double x_step = crop_width / dst_width;
  double y_step = crop_height / dst_height;

  double ypos = start_y;
  double xpos = start_x;

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      int window_end_x = int(xpos + filter_support);
      int window_end_y = int(ypos + filter_support);

      if (window_end_x >= src_width)
        window_end_x = src_width-1;

      if (window_end_y >= src_height)
        window_end_y = src_height-1;

      int window_begin_x = window_end_x - filter_size + 1;
      int window_begin_y = window_end_y - filter_size + 1;

      if (window_begin_x < 0)
        window_begin_x = 0;

      if (window_begin_y < 0)
        window_begin_y = 0;

      float result = 0.0;
      float divider = 0.0;

      double current_x = xpos < 0 ? 0 : (xpos > (src_width-1) ? (src_width-1) : xpos);
      double current_y = ypos < 0 ? 0 : (ypos > (src_height-1) ? (src_height-1) : ypos);

      int window_y = window_begin_y;
      int window_x = window_begin_x;
      for (int ly = 0; ly < filter_size; ly++) {
        for (int lx = 0; lx < filter_size; lx++) {
          double dx = (current_x-window_x)*(current_x-window_x);
          double dy = (current_y-window_y)*(current_y-window_y);

          double dist_sqr = dx + dy;

          float src_data = (src+window_y*src_pitch)[window_x];
          float factor = func->GetFactor(dist_sqr);
          result += src_data * factor;
          divider += factor;

          window_x++;
        }

        window_x = window_begin_x;
        window_y++;
      }
      
      dst[x] = min(255, max(0, int((result/divider)+0.5)));

      xpos += x_step;
    }

    dst += dst_pitch;
    ypos += y_step;
    xpos = start_x;
  }
}
