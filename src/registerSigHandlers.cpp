#include "RBGL.hpp"
#include <Rdefines.h>

extern "C"  {
void sigabrt_handler(int isig) {
    Rf_error
        ("internal: RBGL invoked 'abort'; see warnings() and restart R");
}
}
