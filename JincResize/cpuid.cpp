
/************************************************************************/
/* CPUID - because AVS internal doesn't have AVX2 and FMA yet           */
/************************************************************************/

#include "cpuid.h"

#include <intrin.h>
#define IS_BIT_SET(bitfield, bit) ((bitfield) & (1<<(bit)) ? true : false)

static int CPUCheckForExtensions()
{
  int result = 0;
  int cpuinfo[4];

  __cpuid(cpuinfo, 1);
  if (IS_BIT_SET(cpuinfo[3], 25))
    result |= EWARESIZE_SSE;
  if (IS_BIT_SET(cpuinfo[3], 26))
    result |= EWARESIZE_SSE2;
  if (IS_BIT_SET(cpuinfo[2], 0))
    result |= EWARESIZE_SSE3;
  if (IS_BIT_SET(cpuinfo[2], 12))
    result |= EWARESIZE_FMA3;

  // AVX
#if (_MSC_FULL_VER >= 160040219)    // We require VC++2010 SP1 at least
  bool xgetbv_supported = IS_BIT_SET(cpuinfo[2], 27);
  bool avx_supported = IS_BIT_SET(cpuinfo[2], 28);
  if (xgetbv_supported && avx_supported) {
    if ((_xgetbv(_XCR_XFEATURE_ENABLED_MASK) & 0x6ull) == 0x6ull) {
      // Check for avx2
      __cpuid(cpuinfo, 7);
      if (IS_BIT_SET(cpuinfo[1], 5))
        result |= EWARESIZE_AVX2;
    }
  }
#endif

  return result;
}

int get_supported_instruction() {
  static int lCPUExtensionsAvailable = CPUCheckForExtensions();
  return lCPUExtensionsAvailable;
}
