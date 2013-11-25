
#ifndef __EWACORE_H
#define __EWACORE_H

/*
 * Base Class for EWA Resampler Core
 */
class EWACore
{
public:
  virtual float GetSupport() = 0;
  float GetFactor(float dist);

protected:
  virtual float factor(float dist) = 0;
};

#endif
