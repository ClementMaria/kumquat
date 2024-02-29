# A function to add a new module in KUMQUAT

set(KUMQUAT_MODULES_FULL_LIST "")
function(add_kumquat_module file_path)
  set(KUMQUAT_MODULES ${KUMQUAT_MODULES} ${file_path} CACHE INTERNAL "KUMQUAT_MODULES")
  # Required by user_version
  set(KUMQUAT_MODULES_FULL_LIST ${KUMQUAT_MODULES_FULL_LIST} ${file_path} PARENT_SCOPE)
  # Include module headers
  if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/${file_path}/include/)
    include_directories(src/${file_path}/include/)
  endif()

endfunction(add_kumquat_module)

set(KUMQUAT_SUB_DIRECTORIES "${KUMQUAT_SUB_DIRECTORIES};benchmark")
set(KUMQUAT_SUB_DIRECTORIES "${KUMQUAT_SUB_DIRECTORIES};example")
set(KUMQUAT_SUB_DIRECTORIES "${KUMQUAT_SUB_DIRECTORIES};test")

message("++ KUMQUAT_SUB_DIRECTORIES list is:\"${KUMQUAT_SUB_DIRECTORIES}\"")
