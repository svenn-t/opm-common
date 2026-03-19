# - Create config.h based on a list of variables
#
# Synopsis:
#   configure_vars (FILE filename verb varlist)
# where
#   filename      Full path (including name) of config.h
#   verb          WRITE or APPEND if truncating or not
#   varlist       List of variable names that has been defined
#
# In addition, this function will define HAVE_CONFIG_H for the
# following compilations, (only) if the filename is "config.h".
#
# Example:
#   list (APPEND FOO_CONFIG_VARS
#     "HAVE_BAR"
#     "HAVE_BAR_VERSION_2"
#     )
#   configure_vars (
#     FILE  ${PROJECT_BINARY_DIR}/config.h
#     WRITE ${FOO_CONFIG_VARS}
#     )

# Copyright (C) 2012 Uni Research AS
# This file is licensed under the GNU General Public License v3.0

function (configure_vars obj filename verb)
  # this is just to make the syntax look like the build-in commands
  if (NOT ("X Y Z ${obj}" STREQUAL "X Y Z FILE" AND
        (("${verb}" STREQUAL "WRITE") OR ("${verb}" STREQUAL "APPEND"))))
    message (FATAL_ERROR "Syntax error in argument list")
  endif ()

  # truncate the file if the verb was "WRITE"
  if (verb STREQUAL "WRITE")
    file (WRITE "${filename}" "")
  endif (verb STREQUAL "WRITE")

  # whenever we use this, we also signal to the header files that we
  # have "config.h". add this before any other files (known till now)
  # to avoid confusion from other configuration files.
  get_filename_component (_config_path "${filename}" PATH)
  get_filename_component (_config_file "${filename}" NAME)
  if ("${_config_file}" MATCHES "config\\.h(\\..+)?")
    add_definitions (-DHAVE_CONFIG_H=1)
    include_directories (BEFORE "${_config_path}")
  endif ("${_config_file}" MATCHES "config\\.h(\\..+)?")

  # only write the current value of each variable once
  set (_args ${ARGN})
  if (_args)
    list (REMOVE_DUPLICATES _args)
  endif (_args)

  # process each variable
  foreach (_var IN LISTS _args)
    # massage the name to remove source code formatting
    string (REGEX REPLACE "^[\\n\\t\\ ]+" "" _var "${_var}")
    string (REGEX REPLACE "[\\n\\t\\ ]+$" "" _var "${_var}")

    # check for empty variable; variables that are explicitly set to false
    # is not included in this clause
    if ((NOT DEFINED ${_var}) OR ("${${_var}}" STREQUAL "") OR NOT _var)
      file (APPEND "${filename}" "/* #undef ${_var} */\n")
    else ((NOT DEFINED ${_var}) OR ("${${_var}}" STREQUAL ""))
      # write to file using the correct syntax
      if ("${_var}" MATCHES "^HAVE_.*")
        file (APPEND "${filename}" "#define ${_var} 1\n")
      else ()
        file (APPEND "${filename}" "#define ${_var} ${${_var}}\n")
      endif()
    endif ((NOT DEFINED ${_var}) OR ("${${_var}}" STREQUAL "") OR NOT _var)
  endforeach(_var)
endfunction (configure_vars obj filename verb)
