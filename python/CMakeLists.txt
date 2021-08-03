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

find_package (Python3 COMPONENTS Interpreter Development)
message(STATUS "python3 include : ${Python3_INCLUDE_DIRS}")

swig_add_library(_verdandi LANGUAGE python SOURCES verdandi.i)
set_property(TARGET _verdandi PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)
set_property(TARGET _verdandi PROPERTY INCLUDE_DIRECTORIES ${libops_SOURCE_DIR} ${seldon_SOURCE_DIR}/ ${Python_INCLUDE_DIRS})
set_property(TARGET _verdandi PROPERTY SWIG_COMPILE_DEFINITIONS SELDON_EXCLUDE_FROM_EXPLICIT_INSTANTIATION VERDANDI_DEBUG_LEVEL_4 VERDANDI_PYTHON_VERSION=${Python_VERSION_MAJOR}.${Python_VERSION_MINOR})
target_link_libraries(_verdandi PUBLIC verdandi ops)
target_include_directories(_verdandi PRIVATE ${Python3_INCLUDE_DIRS} )

if (DEFINED PYTHON_SITE_PACKAGES)
    set (VERDANDI_PYTHON_MODULE_PATH ${PYTHON_SITE_PACKAGES})
else ()
    execute_process (COMMAND ${Python3_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
                      OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
                      RESULT_VARIABLE _PYTHON_pythonlib_result
                      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if (_PYTHON_pythonlib_result)
      message (SEND_ERROR "Could not run ${Python3_EXECUTABLE}")
    endif ()

    get_filename_component (_ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE)
    file (RELATIVE_PATH VERDANDI_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH})
endif ()
set (VERDANDI_PYTHON${PYTHON_VERSION_MAJOR}_MODULE_PATH ${VERDANDI_PYTHON_MODULE_PATH})
message(STATUS "[verdandi] python module path: ${VERDANDI_PYTHON_MODULE_PATH}")
install(TARGETS _verdandi DESTINATION ${VERDANDI_PYTHON_MODULE_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/verdandi.py DESTINATION ${VERDANDI_PYTHON_MODULE_PATH})