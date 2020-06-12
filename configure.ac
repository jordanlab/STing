# Autoconf requirements
# Force autoconf to be at least this version number:
AC_PREREQ([2.69])
# Initialize, but use more options.  Note parameters brackets and whitespace.
AC_INIT([STing], [1.0.1], [hspitia@gatech.edu])
# Safety check - list a source file that wouldn't be in other directories:
AC_CONFIG_SRCDIR([src/Typer.cpp])
# Put configuration results here, so we can easily #include them:
AC_CONFIG_HEADERS([src/config.h])

# Put autotools auxiliary files in a subdir, so they don't clutter top dir:
AC_CONFIG_AUX_DIR([build-aux])

AC_CONFIG_MACRO_DIR([m4])

# Enable "automake" to simplify creating makefiles:
AM_INIT_AUTOMAKE([1.11 subdir-objects -Wall -Werror])
# AC_CONFIG_FILES([Makefile])

# Look for a C compiler:
#AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CPP
AC_PROG_AWK
AC_PROG_MKDIR_P


# Perhaps check for: (more) programs, libraries, header files, types,
# structures, compiler characteristics, library functions, system services.

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lrt':
AC_CHECK_LIB([rt], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h limits.h malloc.h stddef.h stdlib.h string.h sys/ioctl.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_CHECK_FUNCS([bzero clock_gettime fchdir fdatasync floor ftruncate gettimeofday memmove memset mkdir munmap pow realpath rmdir sqrt strchr strdup strerror strtol])

# Do final output.
AC_OUTPUT(Makefile src/Makefile)