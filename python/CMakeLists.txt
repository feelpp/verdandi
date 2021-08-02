find_package(SWIG 4.0 COMPONENTS python OPTIONAL_COMPONENTS CXX)
if(SWIG_FOUND)
  message("SWIG found: ${SWIG_EXECUTABLE}")
  if(NOT SWIG_CXX_FOUND)
    message(WARNING "SWIG CXX bindings cannot be generated")
  endif()
endif()
INCLUDE(${SWIG_USE_FILE})
set (UseSWIG_TARGET_NAME_PREFERENCE STANDARD)

set_property(SOURCE verdandi.i PROPERTY CPLUSPLUS ON )


swig_add_library(_verdandi LANGUAGE python SOURCES verdandi.i)
set_property(TARGET _verdandi PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
set_property(TARGET _verdandi PROPERTY INCLUDE_DIRECTORIES ${libops_SOURCE_DIR} ${seldon_SOURCE_DIR}/ )
set_property(TARGET _verdandi PROPERTY SWIG_COMPILE_DEFINITIONS SELDON_EXCLUDE_FROM_EXPLICIT_INSTANTIATION VERDANDI_DEBUG_LEVEL_4)
target_link_libraries(_verdandi PUBLIC verdandi ops)