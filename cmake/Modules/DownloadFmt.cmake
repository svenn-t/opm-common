include(FetchContent)

macro(DownloadFmt target)
  if(NOT fmt_POPULATED)
    FetchContent_Declare(fmt
                         DOWNLOAD_EXTRACT_TIMESTAMP ON
                         URL https://github.com/fmtlib/fmt/archive/refs/tags/11.0.2.tar.gz
                         URL_HASH SHA512=47ff6d289dcc22681eea6da465b0348172921e7cafff8fd57a1540d3232cc6b53250a4625c954ee0944c87963b17680ecbc3ea123e43c2c822efe0dc6fa6cef3)
    FetchContent_Populate(fmt)
  endif()

  # We do not want to use the target directly as that means it ends up
  # in our depends list for downstream modules and the installation list.
  # Instead, we just download and use header only mode.
  target_compile_definitions(${target} PUBLIC $<BUILD_INTERFACE:FMT_HEADER_ONLY>)
  target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${fmt_SOURCE_DIR}/include>)

  # Target used by binaries
  add_library(fmt::fmt INTERFACE IMPORTED)
  target_compile_definitions(fmt::fmt INTERFACE FMT_HEADER_ONLY)
  target_include_directories(fmt::fmt INTERFACE ${fmt_SOURCE_DIR}/include)

  # Required for super-build
  if(NOT PROJECT_IS_TOP_LEVEL)
    set(fmt_POPULATED ${fmt_POPULATED} PARENT_SCOPE)
    set(fmt_SOURCE_DIR ${fmt_SOURCE_DIR} PARENT_SCOPE)
  endif()
endmacro()
