AC_INIT([RGBL], 1.59.2) 
AC_CONFIG_MACRO_DIR([./m4])

# This kludge is to prevent R from inserting a warning 
# about ignoring the value of R_HOME into the BOOSTLOC variable 

# Set R to the default value from the path
R=`which R`


# If R_HOME is set, hide it and use a different name
if test -n "${R_HOME}";
then
  FAKE_R_HOME=${R_HOME}
  unset R_HOME
  R="${FAKE_R_HOME}/bin/R"
fi 

# Check for package
AX_R_PACKAGE([BH],,[${R}])
BOOSTLOC=`${R} --slave -e 'cat(system.file(package="BH", "include")))'`

# If R_HOME was hidden, restore it
if test -n "${FAKE_R_HOME}";
then
  R_HOME=${FAKE_R_HOME}
  unset FAKE_R_HOME
fi

###  PKGLIBS='-L/llvm/clang+llvm-8.0.0-x86_64-linux-gnu-ubuntu-18.04/lib'
### AC_SUBST(PKGLIBS)

BOOSTFLAGS=-I${BOOSTLOC}
AC_SUBST(BOOSTFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

