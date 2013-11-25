
#include <windows.h>
#pragma warning(disable: 4512 4244 4100)
#include "avisynth.h"
#pragma warning(default: 4512 4244 4100)

#include "FilteredEWAResize.h"
#include "JincFilter.h"

AVSValue __cdecl Create_EWAResizer(PClip clip, int target_width, int target_height, const AVSValue* args, EWACore* func, IScriptEnvironment* env)
{
  const VideoInfo& vi = clip->GetVideoInfo();
  const double crop_left = args[0].AsDblDef(0), crop_top = args[1].AsDblDef(0);

  double crop_width = args[2].AsDblDef(vi.width), crop_height = args[3].AsDblDef(vi.height);
  // Crop style syntax
  if (crop_width  <= 0.0) crop_width  = vi.width  - crop_left + crop_width;
  if (crop_height <= 0.0) crop_height = vi.height - crop_top  + crop_height;

  return new FilteredEWAResize(clip, target_width, target_height, crop_left, crop_top, crop_width, crop_height, func, env);
}

template<int tap>
AVSValue __cdecl Create_JincResizer(AVSValue args, void* user_data, IScriptEnvironment* env)
{
  JincFilter* jinc = new JincFilter(tap);
  return Create_EWAResizer(args[0].AsClip(), args[1].AsInt(), args[2].AsInt(), &args[3], jinc, env);
}

const AVS_Linkage *AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
  AVS_linkage = vectors;

  env->AddFunction("Jinc36Resize", "cii[src_left]f[src_top]f[src_width]f[src_height]f", Create_JincResizer<3>, 0);
  env->AddFunction("Jinc64Resize", "cii[src_left]f[src_top]f[src_width]f[src_height]f", Create_JincResizer<4>, 0);
  env->AddFunction("Jinc256Resize", "cii[src_left]f[src_top]f[src_width]f[src_height]f", Create_JincResizer<8>, 0);

  return "Thank you madshi for helping me.";
}
