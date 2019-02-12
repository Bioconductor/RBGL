#include <Rdefines.h>
// #include "kludge.hpp"

extern "C"  {
void _bgl_abort()
{
    Rf_error
        ("internal: RBGL invoked 'abort'; see warnings() and restart R");
}
}



