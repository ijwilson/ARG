#ifndef MYENUMS_H__
#define MYENUMS_H__

#include <stdio.h>

#ifdef __cplusplus
    extern "C" {
#endif
/* new stuff*/
typedef enum partype {integer,floating,dist,character} par_type;
typedef enum volum {quiet, loud} volume;

#ifdef __cplusplus
    }
#endif


#endif
