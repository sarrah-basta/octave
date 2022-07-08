#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Defining the content field */
struct _N_VectorContent_Col{
    int length;  /* vector length       */
    int own_data;           /* data ownership flag */
    double *data;            /* data array, will be our column vector */
};

typedef struct _N_VectorContent_Col *N_VectorContent_Col;

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Col (int length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Col(sunindextype length, double *v_data, SUNContext sunctx);

SUNDIALS_STATIC_INLINE
N_Vector_ID N_VGetVectorID_Col(N_Vector v)
{
  return SUNDIALS_NVEC_CUSTOM;
}

#ifdef __cplusplus
}
#endif
