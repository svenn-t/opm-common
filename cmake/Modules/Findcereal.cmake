include(FetchContent)

if(NOT cereal_POPULATED)
  FetchContent_Declare(
      cereal
      DOWNLOAD_EXTRACT_TIMESTAMP ON
      URL https://github.com/USCiLab/cereal/archive/refs/tags/v1.3.2.tar.gz
      URL_HASH SHA512=98d306d6292789129675f1c5c5aedcb90cfcc1029c4482893a8f9b23f3c9755e5ed4762d7a528f215345cae6392e87cd8d89467115b6f031b41c8673d6b4b109)
  # FetchContent_MakeAvailable(cereal)
  FetchContent_Populate(cereal)
endif()

# Add include folder 
include_directories(${cereal_SOURCE_DIR}/include)
set(cereal_POPULATED ${cereal_POPULATED} PARENT_SCOPE)
set(cereal_SOURCE_DIR ${cereal_SOURCE_DIR} PARENT_SCOPE)
