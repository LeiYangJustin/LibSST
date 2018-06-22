#pragma once

//// ?
//#pragma warning(disable: 4251)
//#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

#ifdef GENERAL_TOOLS_EXPORTS
#define GENERAL_TOOLS_API extern "C" __declspec(dllexport)
#define GENERAL_TOOLS_CLASS __declspec(dllexport)
#define GENERAL_TOOLS_TEMPLATE __declspec(dllexport)
#else
#define GENERAL_TOOLS_API extern "C" __declspec(dllimport)
#define GENERAL_TOOLS_CLASS __declspec(dllimport)
#define GENERAL_TOOLS_TEMPLATE
#endif