/*********************** dSPACE target specific file *************************

   Header file Pierre_Controller_LTI_trc_ptr.h:

   Declaration of function that initializes the global TRC pointers

   RTI1202 25.2 (02-Nov-2025)
   Wed Nov 18 15:35:18 2026

   Copyright 2026, dSPACE GmbH. All rights reserved.

 *****************************************************************************/
#ifndef RTI_HEADER_Pierre_Controller_LTI_trc_ptr_h_
#define RTI_HEADER_Pierre_Controller_LTI_trc_ptr_h_

/* Include the model header file. */
#include "Pierre_Controller_LTI.h"
#include "Pierre_Controller_LTI_private.h"
#ifdef EXTERN_C
#undef EXTERN_C
#endif

#ifdef __cplusplus
#define EXTERN_C                       extern "C"
#else
#define EXTERN_C                       extern
#endif

/*
 *  Declare the global TRC pointers
 */
EXTERN_C volatile real_T *p_0_Pierre_Controller_LTI_real_T_0;
EXTERN_C volatile real_T *p_1_Pierre_Controller_LTI_real_T_0;
EXTERN_C volatile real_T *p_2_Pierre_Controller_LTI_real_T_0;
EXTERN_C volatile int_T *p_2_Pierre_Controller_LTI_int_T_2;
EXTERN_C volatile real_T *p_3_Pierre_Controller_LTI_real_T_0;

/*
 *  Declare the general function for TRC pointer initialization
 */
EXTERN_C void Pierre_Controller_LTI_rti_init_trc_pointers(void);

#endif                         /* RTI_HEADER_Pierre_Controller_LTI_trc_ptr_h_ */
