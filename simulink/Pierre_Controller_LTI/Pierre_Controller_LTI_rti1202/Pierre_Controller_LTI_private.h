/*
 * Pierre_Controller_LTI_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Pierre_Controller_LTI".
 *
 * Model version              : 9.4
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C source code generated on : Wed Nov 18 15:35:18 2026
 *
 * Target selection: rti1202.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef Pierre_Controller_LTI_private_h_
#define Pierre_Controller_LTI_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "Pierre_Controller_LTI_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

/* Used by FromWorkspace Block: '<S17>/SineSweep' */
#ifndef rtInterpolate
# define rtInterpolate(v1,v2,f1,f2)    (((v1)==(v2))?((double)(v1)): (((f1)*((double)(v1)))+((f2)*((double)(v2)))))
#endif

#ifndef rtRound
# define rtRound(v)                    ( ((v) >= 0) ? floor((v) + 0.5) : ceil((v) - 0.5) )
#endif

extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_16;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_13;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_14;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_15;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_1;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_2;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_3;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_4;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_5;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_6;
extern DacCl1AnalogOutSDrvObject *pRTIDacC1AnalogOut_Ch_9;
extern AdcCl2AnalogInSDrvObject *pRTIAdcC2AnalogIn_Ch_1;
extern AdcCl2AnalogInSDrvObject *pRTIAdcC2AnalogIn_Ch_2;
extern AdcCl2AnalogInSDrvObject *pRTIAdcC2AnalogIn_Ch_3;
extern AdcCl2AnalogInSDrvObject *pRTIAdcC2AnalogIn_Ch_4;
extern AdcCl2AnalogInSDrvObject *pRTIAdcC2AnalogIn_Ch_5;
extern AdcCl2AnalogInSDrvObject *pRTIAdcC2AnalogIn_Ch_6;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_9;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_10;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_11;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_17;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_18;
extern AdcCl1AnalogInSDrvObject *pRTIAdcC1AnalogIn_Ch_19;

#endif                                 /* Pierre_Controller_LTI_private_h_ */
