#pragma once

// ?
#pragma warning(disable: 4251)
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

#ifdef DEF_ALGCOLLE_EXPORTS
#define DEF_ALGCOLLE_API extern "C" __declspec(dllexport)
#define DEF_ALGCOLLE_CLASS __declspec(dllexport)
#define DEF_ALGCOLLE_TEMPLATE __declspec(dllexport)
#else
#define DEF_ALGCOLLE_API extern "C" __declspec(dllimport)
#define DEF_ALGCOLLE_CLASS __declspec(dllimport)
#define DEF_ALGCOLLE_TEMPLATE
#endif