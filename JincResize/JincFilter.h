
#ifndef __JINC_FILTER_H
#define __JINC_FILTER_H

#include "EWACore.h"

class JincFilter : public EWACore
{
public:
  JincFilter(int taps);
  ~JincFilter() {}
  float GetSupport() { return support; }

protected:
  float factor(float dist);

private:
  int taps;
  float support;
  float window_factor;
};

#endif
