
#include <math.h>
#include "FilteredEWAResize.h"

// Intrinsics
#include "smmintrin.h"

static void resize_plane_c(EWACore* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                           int src_width, int src_height, int dst_width, int dst_height,
                           double crop_left, double crop_top, double crop_width, double crop_height)
{
  float filter_support = func->GetSupport();
  int filter_size = (int) ceil(filter_support * 2.0);

  float start_x = (float) (crop_left + (crop_width - dst_width) / (dst_width*2));
  float start_y = (float) (crop_top + (crop_height - dst_height) / (dst_height*2));

  float x_step = (float) (crop_width / dst_width);
  float y_step = (float) (crop_height / dst_height);

  float ypos = start_y;
  float xpos = start_x;

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      // Here, the window_*** variable specified a begin/size/end
      // of EWA window to process.
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

      // This is the location of current target pixel in source pixel
      float current_x = clamp((float) 0., xpos, src_width-(float) 1.);
      float current_y = clamp((float) 0., ypos, src_height-(float) 1.);

      int window_y = window_begin_y;
      int window_x = window_begin_x;

      const BYTE* src_begin = src + window_x + window_y*src_pitch;
      for (int ly = 0; ly < filter_size; ly++) {
        const BYTE* src_current = src_begin;
        src_begin += src_pitch;

        for (int lx = 0; lx < filter_size; lx++) {

          float dx = (current_x-window_x)*(current_x-window_x);
          float dy = (current_y-window_y)*(current_y-window_y);

          float dist_sqr = dx + dy;

          float src_data = src_current[lx];
          float factor = func->GetFactor(dist_sqr);

          result += src_data * factor;
          divider += factor;

          window_x++;
        }

        window_x = window_begin_x;
        window_y++;
      }

      dst[x] = clamp(0, int((result/divider)+0.5), 255);

      xpos += x_step;
    }

    dst += dst_pitch;
    ypos += y_step;
    xpos = start_x;
  }
}


template<int filter_size>
static void resize_plane_sse(EWACore* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                             int src_width, int src_height, int dst_width, int dst_height,
                             double crop_left, double crop_top, double crop_width, double crop_height)
{
  float filter_support = func->GetSupport();
  int filter_size2 = ceil(filter_support * 2.0);

  if (filter_size2 != filter_size) {
    throw 0;
  }

  float start_x = (float) (crop_left + (crop_width - dst_width) / (dst_width*2));
  float start_y = (float) (crop_top + (crop_height - dst_height) / (dst_height*2));

  float x_step = (float) (crop_width / dst_width);
  float y_step = (float) (crop_height / dst_height);

  float ypos = start_y;
  float xpos = start_x;

  __m128 zero = _mm_setzero_ps();
  __m128i zeroi = _mm_setzero_si128();

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

      __m128 result = zero;
      __m128 divider = zero;

      float current_x = clamp((float) 0., xpos, src_width-(float) 1.);
      float current_y = clamp((float) 0., ypos, src_height-(float) 1.);
      __m128 curr_x = _mm_set1_ps(current_x);
      __m128 curr_y = _mm_set1_ps(current_y);

      int window_y = window_begin_y;
      int window_x = window_begin_x;

      const BYTE* src_begin = src + window_x + window_y*src_pitch;

      for (int ly = 0; ly < filter_size; ly++) {
        const BYTE* src_current = src_begin;
        src_begin += src_pitch;

        // Same stuff
        __m128 wind_y = _mm_set1_ps(window_y);

        // Process 4 pixel per internal loop (SSE 128bit + single-precision)
        for (int lx = 0; lx < filter_size; lx+=4) {
          // ---------------------------------
          // Calculate coeff

          __declspec(align(16))
          float factor[4];

          __m128 wind_x = _mm_setr_ps(window_x, window_x+1, window_x+2, window_x+3);

          __m128 dx = _mm_sub_ps(curr_x, wind_x);
          __m128 dy = _mm_sub_ps(curr_y, wind_y);

          dx = _mm_mul_ps(dx, dx);
          dy = _mm_mul_ps(dy, dy);

          __m128 dist_simd = _mm_add_ps(dx, dy);

          _mm_store_ps(factor, dist_simd);

          for (int i = 0; i < 4; i++) {
            factor[i] = func->GetFactor(factor[i]);
          }

          __m128 factor_simd = _mm_load_ps(factor);
          window_x += 4;
          // ---------------------------------

          // ---------------------------------
          // Load data and convert to floating point
          __m128i data_epu8 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src_current+lx));
          __m128i data_ep16 = _mm_unpacklo_epi8(data_epu8, zeroi);
          __m128i data_ep32 = _mm_unpacklo_epi16(data_ep16, zeroi);
          __m128  data      = _mm_cvtepi32_ps(data_ep32);
          // ---------------------------------

          // ---------------------------------
          // Process data
          __m128 res = _mm_mul_ps(data, factor_simd);

          result = _mm_add_ps(result, res);
          divider = _mm_add_ps(divider, factor_simd);
          // ---------------------------------
        }

        window_x = window_begin_x;
        window_y++;
      }

      // Add to single float at the lower bit
      result = _mm_hadd_ps(result, zero);
      result = _mm_hadd_ps(result, zero);

      divider = _mm_hadd_ps(divider, zero);
      divider = _mm_hadd_ps(divider, zero);

      result = _mm_div_ss(result, divider);

      int result_i = _mm_cvtss_si32(result);

      dst[x] = clamp(0, result_i, 255);

      xpos += x_step;
    }

    dst += dst_pitch;
    ypos += y_step;
    xpos = start_x;
  }
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
  if (env->GetCPUFlags() & CPUF_SSE3) {
    switch (int(ceil(func->GetSupport() * 2.0))) {

#define size(n)  \
    case n: resizer = resize_plane_sse<n>; break;

      size(3); size(5); size(7); size(9);
      size(11); size(13); size(15); size(17);

#undef size

    default:
      env->ThrowError("JincResize: Internal error; filter size not supported");
    }
  } else {
    resizer = resize_plane_c;
  }
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
