
#ifndef __EWARESIZE_H
#define __EWARESIZE_H

#include "cpuid.h"
#include "EWAResizerStruct.h"

// Intrinsics
#if defined(USE_SSE3) || defined(USE_SSE2)
# include "smmintrin.h"
#endif

#if defined(USE_AVX2) || defined(USE_FMA3)
# include "immintrin.h"
#endif

static void init_coeff_table(EWACore* func, EWAPixelCoeff* out, int quantize_x, int quantize_y,
                             int src_width, int src_height, int dst_width, int dst_height,
                             double crop_left, double crop_top, double crop_width, double crop_height)
{
  const float filter_support = func->GetSupport();
  const int filter_size = (int) ceil(filter_support * 2.0);

  out->filter_size = filter_size;
  out->quantize_x = quantize_x;
  out->quantize_y = quantize_y;
  out->coeff_stripe = ((filter_size + 7) / 8) * 8;

  // Allocate metadata
  out->meta = new EWAPixelCoeffMeta[dst_width * dst_height];

  // Alocate factor map
  if (quantize_x*quantize_y > 0)
    out->factor_map = new float*[quantize_x*quantize_y];
  else
    out->factor_map = nullptr;

  // Calculate how many pixel do we need to calculate pixel for

  const float x_step = (float) (crop_width / dst_width);
  const float y_step = (float) (crop_height / dst_height);

  // Calculate number of border pixels
  int total_coeff = 0;

  const int coeff_per_pixel = out->coeff_stripe * filter_size;
  
  if (quantize_x*quantize_y > 0) {
    const int border_size_x = int(ceil(ceil(filter_support) / x_step));
    const int border_size_y = int(ceil(ceil(filter_support) / y_step));

    int border_count = 0;
    border_count += (2 * border_size_x) * dst_height;
    border_count += (2 * border_size_y) * (dst_width - 2 * border_size_x);

    const int pixel_count = border_count + quantize_x*quantize_y;

    total_coeff = pixel_count * coeff_per_pixel;
  } else {
    total_coeff = dst_width * dst_height * coeff_per_pixel;
  }

  out->factor = (float*) _aligned_malloc(total_coeff * sizeof(float), 64); // align to cache line

  // Zeroed memory
  memset(out->factor, 0, total_coeff * sizeof(float));
  if (quantize_x*quantize_y > 0)
    memset(out->factor_map, 0, quantize_x*quantize_y * sizeof(float));
  memset(out->meta, 0, dst_width * dst_height);
}

static void delete_coeff_table(EWAPixelCoeff* out) {
  if (out == nullptr)
    return;

  _aligned_free(out->factor);
  delete[] out->meta;
  delete[] out->factor_map;
}


/************************************************************************/
/* Coefficient table generation                                         */
/************************************************************************/

#pragma intel optimization_parameter target_arch=sse
static void generate_coeff_table_c(EWACore* func, EWAPixelCoeff* out, int quantize_x, int quantize_y,
                                   int src_width, int src_height, int dst_width, int dst_height,
                                   double crop_left, double crop_top, double crop_width, double crop_height)
{
  const float filter_support = func->GetSupport();
  const int filter_size = (int) ceil(filter_support * 2.0);

  const float start_x = (float) (crop_left + (crop_width - dst_width) / (dst_width * 2));
  const float start_y = (float) (crop_top + (crop_height - dst_height) / (dst_height * 2));

  const float x_step = (float) (crop_width / dst_width);
  const float y_step = (float) (crop_height / dst_height);

  float ypos = start_y;
  float xpos = start_x;

  // Initialize EWAPixelCoeff data structure
  init_coeff_table(func, out, quantize_x, quantize_y, src_width, src_height, dst_width, dst_height, crop_left, crop_top, crop_width, crop_height);
  float* coeff_pointer = out->factor;

  // Use to advance the coeff pointer
  const int coeff_per_pixel = out->coeff_stripe * filter_size;

  bool is_quantizing = (quantize_x * quantize_y) > 0;

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      bool is_border = false;

      EWAPixelCoeffMeta* meta = &out->meta[y*dst_width + x];

      // Here, the window_*** variable specified a begin/size/end
      // of EWA window to process.
      int window_end_x = int(xpos + filter_support);
      int window_end_y = int(ypos + filter_support);

      if (window_end_x >= src_width) {
        window_end_x = src_width - 1;
        is_border = true;
      }

      if (window_end_y >= src_height) {
        window_end_y = src_height - 1;
        is_border = true;
      }

      int window_begin_x = window_end_x - filter_size + 1;
      int window_begin_y = window_end_y - filter_size + 1;

      if (window_begin_x < 0) {
        window_begin_x = 0;
        is_border = true;
      }

      if (window_begin_y < 0) {
        window_begin_y = 0;
        is_border = true;
      }

      meta->start_x = window_begin_x;
      meta->start_y = window_begin_y;

      // Quantize xpos and ypos
      const int quantized_x_int = int(xpos * quantize_x + 0.5f);
      const int quantized_y_int = int(ypos * quantize_y + 0.5f);
      const int quantized_x_value = quantized_x_int % quantize_x;
      const int quantized_y_value = quantized_y_int % quantize_y;
      const float quantized_xpos = float(quantized_x_int) / quantize_x;
      const float quantized_ypos = float(quantized_y_int) / quantize_y;

      //if (is_border && x > 500) __asm int 3;

      if (!is_border && out->factor_map[quantized_y_value*quantize_x + quantized_x_value] != nullptr) {
        // Not border pixel and already have coefficient calculated at this quantized position
        meta->coeff = out->factor_map[quantized_y_value*quantize_x + quantized_x_value];
      } else { // then need computation
        float divider = 0.0;

        // This is the location of current target pixel in source pixel
        // Quantized
        float current_x = clamp(0.f, is_border ? xpos : quantized_xpos, src_width - 1.f);
        float current_y = clamp(0.f, is_border ? ypos : quantized_ypos, src_height - 1.f);

        // Windowing position
        int window_y = window_begin_y;
        int window_x = window_begin_x;

        // First loop calculate coeff
        float* curr_factor_ptr = coeff_pointer;

        for (int ly = 0; ly < filter_size; ly++) {
          for (int lx = 0; lx < filter_size; lx++) {
            // Euclidean distance to sampling pixel
            float dx = (current_x - window_x)*(current_x - window_x);
            float dy = (current_y - window_y)*(current_y - window_y);
            float dist_sqr = dx + dy;

            float factor = func->GetFactor(dist_sqr);

            curr_factor_ptr[lx] = factor;
            divider += factor;

            window_x++;
          } // for (lx)

          curr_factor_ptr += out->coeff_stripe;

          window_x = window_begin_x;
          window_y++;
        } // for (ly)

        // Second loop to divide the coeff
        curr_factor_ptr = coeff_pointer;
        for (int ly = 0; ly < filter_size; ly++) {
          for (int lx = 0; lx < filter_size; lx++) {
            curr_factor_ptr[lx] /= divider;
          }

          curr_factor_ptr += out->coeff_stripe;
        }

        if (!is_border) {
          // Save factor to table
          out->factor_map[quantized_y_value*quantize_x + quantized_x_value] = coeff_pointer;
        }

        meta->coeff = coeff_pointer;
        coeff_pointer += coeff_per_pixel;

      } // if (need_compute)

      xpos += x_step;
    } // for (xpos)

    ypos += y_step;
    xpos = start_x;
  } // for (ypos)
}


/************************************************************************/
/* Planar resampling with coeff table                                   */
/************************************************************************/

#pragma intel optimization_parameter target_arch=sse
static void resize_plane_c(EWAPixelCoeff* coeff, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                                 int src_width, int src_height, int dst_width, int dst_height,
                                 double crop_left, double crop_top, double crop_width, double crop_height)
{
  EWAPixelCoeffMeta* meta = coeff->meta;

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      const BYTE* src_ptr = src + (meta->start_y * src_pitch) + meta->start_x;
      const float* coeff_ptr = meta->coeff;
      float result = 0.0;

      for (int ly = 0; ly < coeff->filter_size; ly++) {
        for (int lx = 0; lx < coeff->filter_size; lx++) {
          result += src_ptr[lx] * coeff_ptr[lx];
        }
        coeff_ptr += coeff->coeff_stripe;
        src_ptr += src_pitch;
      }

      // Save data
      dst[x] = clamp(0, int(result), 255);

      meta++;
    } // for (x)
    dst += dst_pitch;
  } // for (y)
}

#if defined(USE_SSE2) || defined(USE_SSE3)

#pragma intel optimization_parameter target_arch=sse2
template<int filter_size, int CPU>
static void resize_plane_sse(EWAPixelCoeff* coeff, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                           int src_width, int src_height, int dst_width, int dst_height,
                           double crop_left, double crop_top, double crop_width, double crop_height)
{
  EWAPixelCoeffMeta* meta = coeff->meta;

//  assert(fitler_size == coeff->filter_size);

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      const BYTE* src_ptr = src + (meta->start_y * src_pitch) + meta->start_x;
      const float* coeff_ptr = meta->coeff;
      
      __m128 result = _mm_setzero_ps();
      __m128i zero = _mm_setzero_si128();

      for (int ly = 0; ly < coeff->filter_size; ly++) {
        for (int lx = 0; lx < coeff->filter_size; lx += 4) {
          __m128i src_epu8 = _mm_cvtsi32_si128(*(int*)(src_ptr+lx));
          __m128i src_ep16 = _mm_unpacklo_epi8(src_epu8, zero);
          __m128i src_ep32 = _mm_unpacklo_epi16(src_ep16, zero);

          __m128 src_ps = _mm_cvtepi32_ps(src_ep32);
          __m128 coeff = _mm_load_ps(coeff_ptr + lx);

          __m128 multiplied = _mm_mul_ps(src_ps, coeff);
          result = _mm_add_ps(result, multiplied);
        }
        coeff_ptr += coeff->coeff_stripe;
        src_ptr += src_pitch;
      }

#ifdef USE_SSE3
      // Add to single float at the lower bit
      if (CPU == EWARESIZE_SSE3) {
        __m128 zero_ps = _mm_setzero_ps();
        result = _mm_hadd_ps(result, zero_ps);
        result = _mm_hadd_ps(result, zero_ps);
      } else
#endif
      {
        // use 3xshuffle + 2xadd instead of 2xhadd
        __m128 result1 = result;
        __m128 result2 = _mm_shuffle_ps(result1, result1, _MM_SHUFFLE(1, 0, 3, 2));

        result1 = _mm_add_ps(result1, result2);
        result2 = _mm_shuffle_ps(result1, result1, _MM_SHUFFLE(2, 3, 0, 1));

        result = _mm_add_ss(result1, result2);
      }

      // Convert back to interger + clamp
      __m128i result_i = _mm_cvtps_epi32(result);
      result_i = _mm_packus_epi16(result_i, zero);

      // Save data
      dst[x] = _mm_cvtsi128_si32(result_i);

      meta++;
    } // for (x)
    dst += dst_pitch;
  } // for (y)
}
#endif

#ifdef USE_AVX2

#pragma intel optimization_parameter target_arch=avx
template<int filter_size, int CPU>
static void resize_plane_avx(EWAPixelCoeff* coeff, BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch,
                             int src_width, int src_height, int dst_width, int dst_height,
                             double crop_left, double crop_top, double crop_width, double crop_height)
{
  EWAPixelCoeffMeta* meta = coeff->meta;

  //  assert(fitler_size == coeff->filter_size);

  for (int y = 0; y < dst_height; y++) {
    for (int x = 0; x < dst_width; x++) {
      const BYTE* src_ptr = src + (meta->start_y * src_pitch) + meta->start_x;
      const float* coeff_ptr = meta->coeff;

      __m256 result = _mm256_setzero_ps();
      __m256i zero = _mm256_setzero_si256();

      for (int ly = 0; ly < coeff->filter_size; ly++) {
        for (int lx = 0; lx < coeff->filter_size; lx += 8) {
          __m256i src_epu8 = _mm256_castsi128_si256(_mm_loadl_epi64((const __m128i*)(src_ptr+lx)));
          __m256i src_ep16 = _mm256_unpacklo_epi8(src_epu8, zero);

          // I don't understand AVX2 at all...
          __m256i src_ep16_shuffled = _mm256_permute4x64_epi64(src_ep16, _MM_SHUFFLE(3, 1, 2, 0));

          __m256i src_ep32 = _mm256_unpacklo_epi16(src_ep16_shuffled, zero);

          __m256 src_ps = _mm256_cvtepi32_ps(src_ep32);
          __m256 coeff = _mm256_load_ps(coeff_ptr + lx);

#ifdef USE_FMA3
          if (CPU == EWARESIZE_FMA3 | EWARESIZE_AVX2) {
            result = _mm256_fmadd_ps(src_ps, coeff, result);
          } else
#endif
          {
            __m256 multiplied = _mm256_mul_ps(src_ps, coeff);
            result = _mm256_add_ps(result, multiplied);
          }
        }
        coeff_ptr += coeff->coeff_stripe;
        src_ptr += src_pitch;
      }

      // Add to single float at the lower bit
      __m256 zero_ps = _mm256_setzero_ps();
      result = _mm256_hadd_ps(result, zero_ps);

      // I don't understand AVX2 at all...
      result = _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(result), _MM_SHUFFLE(3, 1, 2, 0)));

      result = _mm256_hadd_ps(result, zero_ps);
      result = _mm256_hadd_ps(result, zero_ps);

      // Convert back to interger + clamp
      __m256i result_i = _mm256_cvtps_epi32(result);
      result_i = _mm256_packus_epi16(result_i, zero);

      // Save data
      dst[x] = _mm_cvtsi128_si32(_mm256_castsi256_si128(result_i));

      meta++;
    } // for (x)
    dst += dst_pitch;
  } // for (y)
}
#endif

#endif // Include guard
