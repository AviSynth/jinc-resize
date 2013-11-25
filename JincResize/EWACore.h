
#ifndef __EWACORE_H
#define __EWACORE_H

template<typename T>
T clamp(T a, T b, T c)
{
  return a > b ? a : (b > c ? c : b);
}

/*
 * Base Class for EWA Resampler Core
 */
class EWACore
{
public:
  virtual float GetSupport() = 0;
  float GetFactor(float dist);
  void InitLutTable();
  void DestroyLutTable();

protected:
  virtual float factor(float dist) = 0;

private:
  float* lut;
  float lut_factor;
};

#endif
