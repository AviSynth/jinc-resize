#include <windows.h>
#include "avisynth.h"

#include "FilteredEWAResize.h"
#include "JincFilter.h"

AVSValue __cdecl Create_Jinc64Resizer(AVSValue args, void* user_data, IScriptEnvironment* env) {
  JincFilter* jinc = new JincFilter(4);
  return new FilteredEWAResize(args[0].AsClip(), args[1].AsInt(), args[2].AsInt(), jinc, env);
}

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit2(IScriptEnvironment* env) {
  env->AddFunction("Jinc64Resize", "cii", Create_Jinc64Resizer, 0);

  return "EWA Resampler";
}
