AC_INIT([RGBL], 1.59.1) 
AC_CONFIG_MACRO_DIR([./m4])
AX_R_PACKAGE([BH])
BOOSTLOC=`Rscript -e 'cat(paste0(system.file(package="BH"),"/include"))'`
BOOSTFLAGS=-I${BOOSTLOC}
AC_SUBST(BOOSTFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT



