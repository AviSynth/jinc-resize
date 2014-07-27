#pragma once
#ifndef __CPUID_H
#define __CPUID_H

int get_supported_instruction();

#define EWARESIZE_SSE   1
#define EWARESIZE_SSE2  2
#define EWARESIZE_SSE3  4
#define EWARESIZE_AVX2  8
#define EWARESIZE_FMA3  16

#endif
