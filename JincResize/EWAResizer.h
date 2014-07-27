
#ifndef __EWARESIZE_H
#define __EWARESIZE_H

// Intrinsics
#include "smmintrin.h"
#ifdef USE_AVX2
# include "immintrin.h"
#endif

#ifdef USE_C

#pragma intel optimization_parameter target_arch=sse
static void resize_plane_c(EWACore* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                           int src_width, int src_height, int dst_width, int dst_height,
                           double crop_left, double crop_top, double crop_width, double crop_height)
{
  float filter_support = func->GetSupport();
  int filter_size = (int) ceil(filter_support * 2.0);

  float start_x = (float) (crop_left + (crop_width - dst_width) / (dst_width * 2));
  float start_y = (float) (crop_top + (crop_height - dst_height) / (dst_height * 2));

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
        window_end_x = src_width - 1;

      if (window_end_y >= src_height)
        window_end_y = src_height - 1;

      int window_begin_x = window_end_x - filter_size + 1;
      int window_begin_y = window_end_y - filter_size + 1;

      if (window_begin_x < 0)
        window_begin_x = 0;

      if (window_begin_y < 0)
        window_begin_y = 0;

      float result = 0.0;
      float divider = 0.0;

      // This is the location of current target pixel in source pixel
      float current_x = clamp((float) 0., xpos, src_width - (float) 1.);
      float current_y = clamp((float) 0., ypos, src_height - (float) 1.);

      int window_y = window_begin_y;
      int window_x = window_begin_x;

      const BYTE* src_begin = src + window_x + window_y*src_pitch;
      for (int ly = 0; ly < filter_size; ly++) {
        const BYTE* src_current = src_begin;
        src_begin += src_pitch;

        for (int lx = 0; lx < filter_size; lx++) {

          float dx = (current_x - window_x)*(current_x - window_x);
          float dy = (current_y - window_y)*(current_y - window_y);

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

      dst[x] = clamp(0, int((result / divider) + 0.5), 255);

      xpos += x_step;
    }

    dst += dst_pitch;
    ypos += y_step;
    xpos = start_x;
  }
}

#endif

template<int filter_size, int CPU>
static void resize_plane_sse(EWACore* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                             int src_width, int src_height, int dst_width, int dst_height,
                             double crop_left, double crop_top, double crop_width, double crop_height)
{
  float filter_support = func->GetSupport();
  int filter_size2 = (int) ceil(filter_support * 2.0);

  if (filter_size2 != filter_size) {
    throw 0;
  }

  float start_x = (float) (crop_left + (crop_width - dst_width) / (dst_width * 2));
  float start_y = (float) (crop_top + (crop_height - dst_height) / (dst_height * 2));

  float x_step = (float) (crop_width / dst_width);
  float y_step = (float) (crop_height / dst_height);

  float ypos = start_y;
  float xpos = start_x;

  __m128 zero = _mm_setzero_ps();
# define zeroi  _mm_castps_si128(zero)

  __m128 src_width_1 = _mm_set1_ps(src_width - (float) 1);
  __m128 src_height_1 = _mm_set1_ps(src_height - (float) 1);

  __m128 lut_factor = _mm_set1_ps(func->lut_factor);
  __m128 lut_size = _mm_set1_ps(LUT_SIZE);

  __m128 window_factor = _mm_setr_ps(0, 1, 2, 3);

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      int window_end_x = int(xpos + filter_support);
      int window_end_y = int(ypos + filter_support);

      if (window_end_x >= src_width)
        window_end_x = src_width - 1;

      if (window_end_y >= src_height)
        window_end_y = src_height - 1;

      int window_begin_x = window_end_x - filter_size + 1;
      int window_begin_y = window_end_y - filter_size + 1;

      if (window_begin_x < 0)
        window_begin_x = 0;

      if (window_begin_y < 0)
        window_begin_y = 0;

      __m128 result = zero;
      __m128 divider = zero;

      __m128 curr_x = _mm_set1_ps(xpos);
      __m128 curr_y = _mm_set1_ps(ypos);

      curr_x = _mm_max_ps(curr_x, zero);
      curr_y = _mm_max_ps(curr_y, zero);

      curr_x = _mm_min_ps(curr_x, src_width_1);
      curr_y = _mm_min_ps(curr_y, src_height_1);

      int window_y = window_begin_y;
      int window_x = window_begin_x;

      const BYTE* src_begin = src + window_x + window_y*src_pitch;

      for (int ly = 0; ly < filter_size; ly++) {
        const BYTE* src_current = src_begin;
        src_begin += src_pitch;

        // Whole-loop stuff
        __m128 wind_y = _mm_set1_ps(window_y);

        // Process 4 pixel per internal loop (SSE 128bit + single-precision)
        for (int lx = 0; lx < filter_size; lx += 4) {

          // ---------------------------------
          // Calculate coeff
          __declspec(align(16)) float factor[4];
          __declspec(align(16)) int factor_pos[4];

          __m128 wind_x = _mm_set1_ps(window_x);
          wind_x = _mm_add_ps(wind_x, window_factor);

          __m128 dx = _mm_sub_ps(curr_x, wind_x);
          __m128 dy = _mm_sub_ps(curr_y, wind_y);

          dx = _mm_mul_ps(dx, dx);
          dy = _mm_mul_ps(dy, dy);

          __m128 dist_simd = _mm_add_ps(dx, dy);
          __m128 lut_pos = _mm_mul_ps(dist_simd, lut_factor);
          lut_pos = _mm_min_ps(lut_pos, lut_size);

          _mm_store_si128(reinterpret_cast<__m128i*>(factor_pos), _mm_cvtps_epi32(lut_pos));

          for (int i = 0; i < 4; i++) {
            factor[i] = func->lut[factor_pos[i]];
          }

          __m128 factor_simd = _mm_load_ps(factor);

          window_x += 4;
          // ---------------------------------

          // ---------------------------------
          // Load data and convert to floating point
          __m128i data_epu8 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src_current + lx));
          __m128i data_ep16 = _mm_unpacklo_epi8(data_epu8, zeroi);
          __m128i data_ep32 = _mm_unpacklo_epi16(data_ep16, zeroi);
          __m128  data = _mm_cvtepi32_ps(data_ep32);
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
      if (CPU == CPUF_SSE3) {

        result = _mm_hadd_ps(result, zero);
        result = _mm_hadd_ps(result, zero);

        divider = _mm_hadd_ps(divider, zero);
        divider = _mm_hadd_ps(divider, zero);

      } else {
        // use 3xshuffle + 2xadd instead of 2xhadd
        __m128 result1 = result;
        __m128 result2 = _mm_shuffle_ps(result1, result1, _MM_SHUFFLE(1, 0, 3, 2));
        __m128 divider1 = divider;
        __m128 divider2 = _mm_shuffle_ps(divider1, divider1, _MM_SHUFFLE(1, 0, 3, 2));

        result1 = _mm_add_ps(result1, result2);
        result2 = _mm_shuffle_ps(result1, result1, _MM_SHUFFLE(2, 3, 0, 1));
        divider1 = _mm_add_ps(divider1, divider2);
        divider2 = _mm_shuffle_ps(divider1, divider1, _MM_SHUFFLE(2, 3, 0, 1));

        result = _mm_add_ss(result1, result2);
        divider = _mm_add_ss(divider1, divider2);
      }

      result = _mm_div_ss(result, divider);

      __m128i result_i = _mm_cvtps_epi32(result);
      result_i = _mm_packus_epi16(result_i, zeroi);

      dst[x] = _mm_cvtsi128_si32(result_i);

      xpos += x_step;
    }

    dst += dst_pitch;
    ypos += y_step;
    xpos = start_x;
  }

#undef zeroi
}

#ifdef USE_AVX2

#pragma intel optimization_parameter target_arch=avx
template<int filter_size>
static void resize_plane_avx(EWACore* func, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                             int src_width, int src_height, int dst_width, int dst_height,
                             double crop_left, double crop_top, double crop_width, double crop_height)
{
  _mm256_zeroupper();

  float filter_support = func->GetSupport();
  int filter_size2 = (int) ceil(filter_support * 2.0);

  if (filter_size2 != filter_size) {
    throw 0;
  }

  float start_x = (float) (crop_left + (crop_width - dst_width) / (dst_width * 2));
  float start_y = (float) (crop_top + (crop_height - dst_height) / (dst_height * 2));

  float x_step = (float) (crop_width / dst_width);
  float y_step = (float) (crop_height / dst_height);

  float ypos = start_y;
  float xpos = start_x;

  __m256 zero = _mm256_setzero_ps();
# define zeroi  _mm256_castps_si256(zero)

  __m256 src_width_1 = _mm256_set1_ps(src_width - (float) 1);
  __m256 src_height_1 = _mm256_set1_ps(src_height - (float) 1);

  __m256 lut_factor = _mm256_set1_ps(func->lut_factor);
  __m256 lut_size = _mm256_set1_ps(LUT_SIZE);

  __m256 window_factor = _mm256_setr_ps(0, 1, 2, 3, 4, 5, 6, 7);

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      int window_end_x = int(xpos + filter_support);
      int window_end_y = int(ypos + filter_support);

      if (window_end_x >= src_width)
        window_end_x = src_width - 1;

      if (window_end_y >= src_height)
        window_end_y = src_height - 1;

      int window_begin_x = window_end_x - filter_size + 1;
      int window_begin_y = window_end_y - filter_size + 1;

      if (window_begin_x < 0)
        window_begin_x = 0;

      if (window_begin_y < 0)
        window_begin_y = 0;

      __m256 result = zero;
      __m256 divider = zero;

      __m256 curr_x = _mm256_set1_ps(xpos);
      __m256 curr_y = _mm256_set1_ps(ypos);

      curr_x = _mm256_max_ps(curr_x, zero);
      curr_y = _mm256_max_ps(curr_y, zero);

      curr_x = _mm256_min_ps(curr_x, src_width_1);
      curr_y = _mm256_min_ps(curr_y, src_height_1);

      int window_y = window_begin_y;
      int window_x = window_begin_x;

      const BYTE* src_begin = src + window_x + window_y*src_pitch;

      for (int ly = 0; ly < filter_size; ly++) {
        const BYTE* src_current = src_begin;
        src_begin += src_pitch;

        // Whole-loop stuff
        __m256 wind_y = _mm256_set1_ps(window_y);

        // Process 4 pixel per internal loop (AVX 256bit + single-precision)
        for (int lx = 0; lx < filter_size; lx += 8) {

          // ---------------------------------
          // Calculate coeff
          __m256 wind_x = _mm256_set1_ps(window_x);
          wind_x = _mm256_add_ps(wind_x, window_factor);

          __m256 dx = _mm256_sub_ps(curr_x, wind_x);
          __m256 dy = _mm256_sub_ps(curr_y, wind_y);

          dx = _mm256_mul_ps(dx, dx);
          dy = _mm256_mul_ps(dy, dy);

          __m256 dist_simd = _mm256_add_ps(dx, dy);
          __m256 lut_pos = _mm256_mul_ps(dist_simd, lut_factor);
          lut_pos = _mm256_min_ps(lut_pos, lut_size);

          __m256 factor_simd = _mm256_i32gather_ps(func->lut, _mm256_cvtps_epi32(lut_pos), 4);

          window_x += 8;
          // ---------------------------------

          // ---------------------------------
          // Load data and convert to floating point
          __m256i data_epu8 = _mm256_castsi128_si256(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(src_current + lx)));
          __m256i data_ep16 = _mm256_unpacklo_epi8(data_epu8, zeroi);
          __m256i data_ep32 = _mm256_unpacklo_epi16(data_ep16, zeroi);
          __m256  data = _mm256_cvtepi32_ps(data_ep32);
          // ---------------------------------

          // ---------------------------------
          // Process data

          result = _mm256_fmadd_ps(data, factor_simd, result);
          divider = _mm256_add_ps(divider, factor_simd);
          // ---------------------------------
        }

        window_x = window_begin_x;
        window_y++;
      }

      // Add to single float at the lower bit

      result = _mm256_hadd_ps(result, zero);
      result = _mm256_hadd_ps(result, zero);

      divider = _mm256_hadd_ps(divider, zero);
      divider = _mm256_hadd_ps(divider, zero);

      __m128 result_128 = _mm_div_ss(
        _mm256_castps256_ps128(result),
        _mm256_castps256_ps128(divider)
        );

      __m128i result_i = _mm_cvtps_epi32(result_128);
      result_i = _mm_packus_epi16(result_i, _mm_setzero_si128());

      dst[x] = _mm_cvtsi128_si32(result_i);

      xpos += x_step;
    }

    dst += dst_pitch;
    ypos += y_step;
    xpos = start_x;
  }

  _mm256_zeroupper();

#undef zeroi
}
#endif

#endif // Include guard
