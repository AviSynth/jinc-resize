
#include <math.h>
#include <stdint.h>
#include "FilteredEWAResize.h"

// Both ICL and MSVC generates much slower AVX code than SSE
#define USE_AVX2
#define USE_SSE2
#define USE_SSE3
// TODO implement above

// Usually floating point in C code get converted to SSE anyway
#define USE_C

#include "EWAResizer.h"

/************************************************************************/
/* EWACore implementation                                               */
/************************************************************************/

const int LUT_SIZE = 8192;

void EWACore::InitLutTable()
{
  lut = new float[LUT_SIZE];

  float filter_end = GetSupport()*GetSupport();
  lut_factor = ((float) (LUT_SIZE - 16)) / filter_end;

  for (int i = 0; i < LUT_SIZE; i++) {
    lut[i] = factor((float)i / lut_factor);
  }
}

void EWACore::DestroyLutTable()
{
  delete[] lut;
}

inline float EWACore::GetFactor(float dist)
{
  int index = int(dist*lut_factor);

  if (index >= LUT_SIZE)
    return 0;

  return lut[index];
}

FilteredEWAResize::FilteredEWAResize(PClip _child, int width, int height, double crop_left, double crop_top, double crop_width, double crop_height, EWACore *func, IScriptEnvironment* env) :
  GenericVideoFilter(_child),
  func(func),
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

  func->InitLutTable();

  // Select the EWA Resize core
  resizer = GetResizer(int(ceil(func->GetSupport() * 2.0)), env);
}

FilteredEWAResize::~FilteredEWAResize()
{
  func->DestroyLutTable();
  delete func;
}

PVideoFrame __stdcall FilteredEWAResize::GetFrame(int n, IScriptEnvironment* env)
{
  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = env->NewVideoFrame(vi);

  try {
    // Luma
    resizer(func,
      dst->GetWritePtr(), src->GetReadPtr(), dst->GetPitch(), src->GetPitch(),
      src_width, src_height, vi.width, vi.height,
      crop_left, crop_top, crop_width, crop_height
      );

    if (!vi.IsY8()) {
      int subsample_w = vi.GetPlaneWidthSubsampling(PLANAR_U);
      int subsample_h = vi.GetPlaneHeightSubsampling(PLANAR_U);

      double div_w = 1 << subsample_w;
      double div_h = 1 << subsample_h;

      resizer(func,
        dst->GetWritePtr(PLANAR_U), src->GetReadPtr(PLANAR_U), dst->GetPitch(PLANAR_U), src->GetPitch(PLANAR_U),
        src_width >> subsample_w, src_height >> subsample_h, vi.width >> subsample_w, vi.height >> subsample_h,
        crop_left / div_w, crop_top / div_h, crop_width / div_w, crop_height / div_h
        );

      resizer(func,
        dst->GetWritePtr(PLANAR_V), src->GetReadPtr(PLANAR_V), dst->GetPitch(PLANAR_V), src->GetPitch(PLANAR_V),
        src_width >> subsample_w, src_height >> subsample_h, vi.width >> subsample_w, vi.height >> subsample_h,
        crop_left / div_w, crop_top / div_h, crop_width / div_w, crop_height / div_h
        );

    }
  } catch (int err) {
    env->ThrowError("JincResize: Internal error, code=", err);
  }

  return dst;
}


EWAResizeCore FilteredEWAResize::GetResizer(int filter_size, IScriptEnvironment* env)
{
#ifdef USE_AVX2
  if (env->GetCPUFlags() & CPUF_AVX) {
    switch (filter_size) {

#define size(n)  \
    case n: return resize_plane_avx<n>; break;

      size(3); size(5); size(7); size(9);
      size(11); size(13); size(15); size(17);

#undef size

    default:
      env->ThrowError("JincResize: Internal error; filter size '%d' is not supported", filter_size);
    }
  }
#endif

  if (env->GetCPUFlags() & CPUF_SSE3) {
    switch (filter_size) {

#define size(n)  \
    case n: return resize_plane_sse<n, CPUF_SSE3>; break;

      size(3); size(5); size(7); size(9);
      size(11); size(13); size(15); size(17);

#undef size

    default:
      env->ThrowError("JincResize: Internal error; filter size '%d' is not supported", filter_size);
    }
  }
  
  if (env->GetCPUFlags() & CPUF_SSE2) {
    switch (filter_size) {

#define size(n)  \
    case n: return resize_plane_sse<n, CPUF_SSE2>; break;

      size(3); size(5); size(7); size(9);
      size(11); size(13); size(15); size(17);

#undef size

    default:
      env->ThrowError("JincResize: Internal error; filter size '%d' is not supported", filter_size);
    }
  }

#ifdef USE_C
  return resize_plane_c;
#else
  env->ThrowError("JincResize: No support instruction set found. Supported: "
# ifdef USE_AVX2
                  "AVX2, "
# endif
# ifdef USE_SSE3
                  "SSE3, "
# endif
# ifdef USE_SSE2
                  "SSE2, "
# endif
# ifdef USE_C
                  "x86"
# else
                  " and no other"
# endif
                  );
#endif
}

#if !defined(USE_AVX2) && !defined(USE_SSE3) && !defined(USE_SSE2) && !defined(USE_C)
# error Why are you compling this plugin without any supported instruction?
#endif

