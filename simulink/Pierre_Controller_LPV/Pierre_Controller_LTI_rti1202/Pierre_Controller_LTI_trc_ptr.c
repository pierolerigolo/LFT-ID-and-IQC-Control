/***************************************************************************

   Source file Pierre_Controller_LTI_trc_ptr.c:

   Definition of function that initializes the global TRC pointers

   RTI1202 25.2 (02-Nov-2025)
   Wed Nov 18 15:35:18 2026

   Copyright 2026, dSPACE GmbH. All rights reserved.

 *****************************************************************************/

/* Include header file. */
#include "Pierre_Controller_LTI_trc_ptr.h"
#include "Pierre_Controller_LTI.h"
#include "Pierre_Controller_LTI_private.h"

/* Compiler options to turn off optimization. */
#if !defined(DS_OPTIMIZE_INIT_TRC_POINTERS)
#ifdef _MCCPPC

#pragma options -nOt -nOr -nOi -nOx

#endif

#ifdef __GNUC__

#pragma GCC optimize ("O0")

#endif

#ifdef _MSC_VER

#pragma optimize ("", off)

#endif
#endif

/* Definition of Global pointers to data type transitions (for TRC-file access) */
volatile real_T *p_0_Pierre_Controller_LTI_real_T_0 = NULL;
volatile real_T *p_1_Pierre_Controller_LTI_real_T_0 = NULL;
volatile real_T *p_2_Pierre_Controller_LTI_real_T_0 = NULL;
volatile int_T *p_2_Pierre_Controller_LTI_int_T_2 = NULL;
volatile real_T *p_3_Pierre_Controller_LTI_real_T_0 = NULL;

/*
 *  Declare the functions, that initially assign TRC pointers
 */
static void rti_init_trc_pointers_0(void);

/* Global pointers to data type transitions are separated in different functions to avoid overloading */
static void rti_init_trc_pointers_0(void)
{
  p_0_Pierre_Controller_LTI_real_T_0 = &Pierre_Controller_LTI_B.SFunction1;
  p_1_Pierre_Controller_LTI_real_T_0 = &Pierre_Controller_LTI_P.Px;
  p_2_Pierre_Controller_LTI_real_T_0 =
    &Pierre_Controller_LTI_DW.ControllerK_DSTATE[0];
  p_2_Pierre_Controller_LTI_int_T_2 =
    &Pierre_Controller_LTI_DW.SineSweep_IWORK.PrevIndex;
  p_3_Pierre_Controller_LTI_real_T_0 = &Pierre_Controller_LTI_Y.F_Shaker;
}

void Pierre_Controller_LTI_rti_init_trc_pointers(void)
{
  rti_init_trc_pointers_0();
}
