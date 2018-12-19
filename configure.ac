AC_INIT([RGBL], 1.59.1) 
BOOSTLOC=`Rscript -e 'cat(paste0(system.file(package="BH"),"/include"))'`
BOOSTFLAGS=-I${BOOSTLOC}
AC_SUBST(BOOSTFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
