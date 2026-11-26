/*
 * Pierre_Controller_LTI.h
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

#ifndef Pierre_Controller_LTI_h_
#define Pierre_Controller_LTI_h_
#ifndef Pierre_Controller_LTI_COMMON_INCLUDES_
#define Pierre_Controller_LTI_COMMON_INCLUDES_
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                              /* Pierre_Controller_LTI_COMMON_INCLUDES_ */

#include "Pierre_Controller_LTI_types.h"
#include <math.h>
#include <string.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T SFunction1;                   /* '<S13>/S-Function1' */
  real_T Clock;                        /* '<S1>/Clock' */
  real_T SFunction1_j;                 /* '<S44>/S-Function1' */
  real_T GainADC_Ch14;                 /* '<S14>/GainADC_Ch14' */
  real_T g2ms2_Wilcoxon;               /* '<S48>/g2ms^2_Wilcoxon' */
  real_T Sum;                          /* '<Root>/Sum' */
  real_T SFunction1_l;                 /* '<S45>/S-Function1' */
  real_T GainADC_Ch15;                 /* '<S14>/GainADC_Ch15' */
  real_T g2ms2_Wilcoxon_a;             /* '<S49>/g2ms^2_Wilcoxon' */
  real_T Sum1;                         /* '<Root>/Sum1' */
  real_T SFunction1_o;                 /* '<S46>/S-Function1' */
  real_T GainADC_Ch16;                 /* '<S14>/GainADC_Ch16' */
  real_T g2ms2_Wilcoxon_o;             /* '<S50>/g2ms^2_Wilcoxon' */
  real_T Sum2;                         /* '<Root>/Sum2' */
  real_T MultiportSwitch2[3];          /* '<S2>/Multiport Switch2' */
  real_T ControllerK[3];               /* '<S2>/Controller K' */
  real_T dSpaceConversion[3];          /* '<S2>/dSpace Conversion' */
  real_T StabilityGain[3];             /* '<S2>/Stability Gain' */
  real_T Gain4;                        /* '<S2>/Gain4' */
  real_T Saturation1;                  /* '<Root>/Saturation1' */
  real_T Clock_c;                      /* '<S18>/Clock' */
  real_T SineSweep;                    /* '<S17>/SineSweep' */
  real_T Product;                      /* '<S17>/Product' */
  real_T PRBS_TMS;                     /* '<Root>/PRBS_TMS' */
  real_T MultiportSwitch2_a;           /* '<Root>/Multiport Switch2' */
  real_T Saturation_TMS;               /* '<Root>/Saturation_TMS' */
  real_T CommandV;                     /* '<S19>/GainDAC_Ch1' */
  real_T Gain3;                        /* '<S2>/Gain3' */
  real_T Add;                          /* '<S2>/Add' */
  real_T Saturation;                   /* '<Root>/Saturation' */
  real_T Gain5;                        /* '<S2>/Gain5' */
  real_T Saturation2;                  /* '<Root>/Saturation2' */
  real_T Gain;                         /* '<S2>/Gain' */
  real_T Gain1;                        /* '<S2>/Gain1' */
  real_T Gain2;                        /* '<S2>/Gain2' */
  real_T SFunction1_m;                 /* '<S20>/S-Function1' */
  real_T GainADC_Ch1_j;                /* '<S3>/GainADC_Ch1' */
  real_T PMA_monitoringA;              /* '<S26>/V2A_PMA_monitoring' */
  real_T SFunction1_f;                 /* '<S21>/S-Function1' */
  real_T GainADC_Ch2;                  /* '<S3>/GainADC_Ch2' */
  real_T PMA_monitoringA_b;            /* '<S27>/V2A_PMA_monitoring' */
  real_T SFunction1_lj;                /* '<S22>/S-Function1' */
  real_T GainADC_Ch3;                  /* '<S3>/GainADC_Ch3' */
  real_T PMA_monitoringA_a;            /* '<S28>/V2A_PMA_monitoring' */
  real_T SFunction1_e;                 /* '<S23>/S-Function1' */
  real_T GainADC_Ch4;                  /* '<S3>/GainADC_Ch4' */
  real_T PMA_monitoringA_bs;           /* '<S29>/V2A_PMA_monitoring' */
  real_T SFunction1_ow;                /* '<S24>/S-Function1' */
  real_T GainADC_Ch5;                  /* '<S3>/GainADC_Ch5' */
  real_T PMA_monitoringA_c;            /* '<S30>/V2A_PMA_monitoring' */
  real_T SFunction1_jx;                /* '<S25>/S-Function1' */
  real_T GainADC_Ch6;                  /* '<S3>/GainADC_Ch6' */
  real_T PMA_monitoringA_cr;           /* '<S31>/V2A_PMA_monitoring' */
  real_T SFunction1_a;                 /* '<S32>/S-Function1' */
  real_T GainADC_ACCE_X;               /* '<S5>/GainADC_ACCE_X' */
  real_T V2g_Wilcoxon;                 /* '<S5>/V2g_Wilcoxon' */
  real_T g2ms2_Wilcoxon_b;             /* '<S39>/g2ms^2_Wilcoxon' */
  real_T SFunction1_c;                 /* '<S33>/S-Function1' */
  real_T GainADC_ACCE_Y;               /* '<S5>/GainADC_ACCE_Y' */
  real_T V2g_Wilcoxon1;                /* '<S5>/V2g_Wilcoxon1' */
  real_T g2ms2_Wilcoxon_d;             /* '<S40>/g2ms^2_Wilcoxon' */
  real_T SFunction1_j0;                /* '<S34>/S-Function1' */
  real_T GainADC_ACCE_Z;               /* '<S5>/GainADC_ACCE_Z' */
  real_T V2g_Wilcoxon2;                /* '<S5>/V2g_Wilcoxon2' */
  real_T g2ms2_Wilcoxon_k;             /* '<S41>/g2ms^2_Wilcoxon' */
  real_T SFunction1_h;                 /* '<S35>/S-Function1' */
  real_T GainADC_GIRO_X;               /* '<S5>/GainADC_GIRO_X' */
  real_T V2g_Wilcoxon3;                /* '<S5>/V2g_Wilcoxon3' */
  real_T SFunction1_en;                /* '<S36>/S-Function1' */
  real_T GainADC_GIRO_Y;               /* '<S5>/GainADC_GIRO_Y' */
  real_T V2g_Wilcoxon4;                /* '<S5>/V2g_Wilcoxon4' */
  real_T SFunction1_cw;                /* '<S37>/S-Function1' */
  real_T GainADC_GIRO_Z;               /* '<S5>/GainADC_GIRO_Z' */
  real_T V2g_Wilcoxon5;                /* '<S5>/V2g_Wilcoxon5' */
  real_T signal;                       /* '<S18>/MATLAB Function' */
  real_T u_exc_cmd;                    /* '<Root>/MATLAB Function' */
  real_T u_exc_dist;                   /* '<Root>/MATLAB Function' */
} B_Pierre_Controller_LTI_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T ControllerK_DSTATE[18];       /* '<S2>/Controller K' */
  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } SineSweep_PWORK;                   /* '<S17>/SineSweep' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } PRBS_TMS_PWORK;                    /* '<Root>/PRBS_TMS' */

  struct {
    int_T PrevIndex;
  } SineSweep_IWORK;                   /* '<S17>/SineSweep' */

  struct {
    int_T PrevIndex;
  } PRBS_TMS_IWORK;                    /* '<Root>/PRBS_TMS' */
} DW_Pierre_Controller_LTI_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T F_Shaker;                     /* '<Root>/F_shaker' */
} ExtY_Pierre_Controller_LTI_T;

/* Parameters (default storage) */
struct P_Pierre_Controller_LTI_T_ {
  real_T Px;                           /* Variable: Px
                                        * Referenced by: '<S2>/Gain3'
                                        */
  real_T Py;                           /* Variable: Py
                                        * Referenced by: '<S2>/Gain4'
                                        */
  real_T Pz;                           /* Variable: Pz
                                        * Referenced by: '<S2>/Gain5'
                                        */
  real_T Sen_AX;                       /* Variable: Sen_AX
                                        * Referenced by: '<S5>/V2g_Wilcoxon'
                                        */
  real_T Sen_AY;                       /* Variable: Sen_AY
                                        * Referenced by: '<S5>/V2g_Wilcoxon1'
                                        */
  real_T Sen_AZ;                       /* Variable: Sen_AZ
                                        * Referenced by: '<S5>/V2g_Wilcoxon2'
                                        */
  real_T Sen_GX;                       /* Variable: Sen_GX
                                        * Referenced by: '<S5>/V2g_Wilcoxon3'
                                        */
  real_T Sen_GY;                       /* Variable: Sen_GY
                                        * Referenced by: '<S5>/V2g_Wilcoxon4'
                                        */
  real_T Sen_GZ;                       /* Variable: Sen_GZ
                                        * Referenced by: '<S5>/V2g_Wilcoxon5'
                                        */
  real_T Sweep_amp;                    /* Variable: Sweep_amp
                                        * Referenced by: '<S17>/Amplitude (Sweep)'
                                        */
  real_T Switch_excit;                 /* Variable: Switch_excit
                                        * Referenced by: '<Root>/Switch Excitation_TMS'
                                        */
  real_T TMS_amp;                      /* Variable: TMS_amp
                                        * Referenced by: '<S18>/Amplitude (Sinus)'
                                        */
  real_T TMS_freq1;                    /* Variable: TMS_freq1
                                        * Referenced by: '<S18>/Frequency 1'
                                        */
  real_T TMS_freq2;                    /* Variable: TMS_freq2
                                        * Referenced by: '<S18>/Frequency 2'
                                        */
  real_T TMS_freq3;                    /* Variable: TMS_freq3
                                        * Referenced by: '<S18>/Frequency 3'
                                        */
  real_T V_sat;                        /* Variable: V_sat
                                        * Referenced by:
                                        *   '<Root>/Saturation'
                                        *   '<Root>/Saturation1'
                                        *   '<Root>/Saturation2'
                                        */
  real_T V_sat_TMS;                    /* Variable: V_sat_TMS
                                        * Referenced by: '<Root>/Saturation_TMS'
                                        */
  real_T Zero_command_Value;           /* Expression: 0
                                        * Referenced by: '<S2>/Zero_command'
                                        */
  real_T GainADC_Ch1_Gain;             /* Expression: 50
                                        * Referenced by: '<S4>/GainADC_Ch1'
                                        */
  real_T Switch_Controller_Value;      /* Expression: 2
                                        * Referenced by: '<S2>/Switch_Controller'
                                        */
  real_T Reference_x_Value;            /* Expression: 0
                                        * Referenced by: '<Root>/Reference_x'
                                        */
  real_T GainADC_Ch14_Gain;            /* Expression: 0.5
                                        * Referenced by: '<S14>/GainADC_Ch14'
                                        */
  real_T g2ms2_Wilcoxon_Gain;          /* Expression: 9.80655
                                        * Referenced by: '<S48>/g2ms^2_Wilcoxon'
                                        */
  real_T Reference_y_Value;            /* Expression: 0
                                        * Referenced by: '<Root>/Reference_y'
                                        */
  real_T GainADC_Ch15_Gain;            /* Expression: 0.5
                                        * Referenced by: '<S14>/GainADC_Ch15'
                                        */
  real_T g2ms2_Wilcoxon_Gain_p;        /* Expression: 9.80655
                                        * Referenced by: '<S49>/g2ms^2_Wilcoxon'
                                        */
  real_T Reference_z_Value;            /* Expression: 0
                                        * Referenced by: '<Root>/Reference_z'
                                        */
  real_T GainADC_Ch16_Gain;            /* Expression: 0.5
                                        * Referenced by: '<S14>/GainADC_Ch16'
                                        */
  real_T g2ms2_Wilcoxon_Gain_i;        /* Expression: 9.80655
                                        * Referenced by: '<S50>/g2ms^2_Wilcoxon'
                                        */
  real_T ControllerK_A[324];           /* Computed Parameter: ControllerK_A
                                        * Referenced by: '<S2>/Controller K'
                                        */
  real_T ControllerK_B[18];            /* Computed Parameter: ControllerK_B
                                        * Referenced by: '<S2>/Controller K'
                                        */
  real_T ControllerK_C[18];            /* Computed Parameter: ControllerK_C
                                        * Referenced by: '<S2>/Controller K'
                                        */
  real_T ControllerK_D;                /* Computed Parameter: ControllerK_D
                                        * Referenced by: '<S2>/Controller K'
                                        */
  real_T ControllerK_InitialCondition; /* Expression: 0
                                        * Referenced by: '<S2>/Controller K'
                                        */
  real_T dSpaceConversion_Gain;        /* Expression: 1/10
                                        * Referenced by: '<S2>/dSpace Conversion'
                                        */
  real_T StabilityGain_Gain;           /* Expression: 1
                                        * Referenced by: '<S2>/Stability Gain'
                                        */
  real_T OFF_TMS_Value;                /* Expression: 0
                                        * Referenced by: '<Root>/OFF_TMS'
                                        */
  real_T GainDAC_Ch1_Gain;             /* Expression: 1/10
                                        * Referenced by: '<S19>/GainDAC_Ch1'
                                        */
  real_T Idonordist_Value;             /* Expression: 0
                                        * Referenced by: '<Root>/Id on or dist'
                                        */
  real_T Gain_Gain;                    /* Expression: 1
                                        * Referenced by: '<S2>/Gain'
                                        */
  real_T Gain1_Gain;                   /* Expression: 1
                                        * Referenced by: '<S2>/Gain1'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1
                                        * Referenced by: '<S2>/Gain2'
                                        */
  real_T GainADC_Ch1_Gain_g;           /* Expression: 10
                                        * Referenced by: '<S3>/GainADC_Ch1'
                                        */
  real_T V2A_PMA_monitoring_Gain;      /* Expression: 1/5
                                        * Referenced by: '<S26>/V2A_PMA_monitoring'
                                        */
  real_T GainADC_Ch2_Gain;             /* Expression: 10
                                        * Referenced by: '<S3>/GainADC_Ch2'
                                        */
  real_T V2A_PMA_monitoring_Gain_d;    /* Expression: 1/5
                                        * Referenced by: '<S27>/V2A_PMA_monitoring'
                                        */
  real_T GainADC_Ch3_Gain;             /* Expression: 10
                                        * Referenced by: '<S3>/GainADC_Ch3'
                                        */
  real_T V2A_PMA_monitoring_Gain_m;    /* Expression: 1/5
                                        * Referenced by: '<S28>/V2A_PMA_monitoring'
                                        */
  real_T GainADC_Ch4_Gain;             /* Expression: 10
                                        * Referenced by: '<S3>/GainADC_Ch4'
                                        */
  real_T V2A_PMA_monitoring_Gain_h;    /* Expression: 1/5
                                        * Referenced by: '<S29>/V2A_PMA_monitoring'
                                        */
  real_T GainADC_Ch5_Gain;             /* Expression: 10
                                        * Referenced by: '<S3>/GainADC_Ch5'
                                        */
  real_T V2A_PMA_monitoring_Gain_n;    /* Expression: 1/5
                                        * Referenced by: '<S30>/V2A_PMA_monitoring'
                                        */
  real_T GainADC_Ch6_Gain;             /* Expression: 10
                                        * Referenced by: '<S3>/GainADC_Ch6'
                                        */
  real_T V2A_PMA_monitoring_Gain_f;    /* Expression: 1/5
                                        * Referenced by: '<S31>/V2A_PMA_monitoring'
                                        */
  real_T GainADC_ACCE_X_Gain;          /* Expression: 10
                                        * Referenced by: '<S5>/GainADC_ACCE_X'
                                        */
  real_T g2ms2_Wilcoxon_Gain_b;        /* Expression: 9.80655
                                        * Referenced by: '<S39>/g2ms^2_Wilcoxon'
                                        */
  real_T GainADC_ACCE_Y_Gain;          /* Expression: 10
                                        * Referenced by: '<S5>/GainADC_ACCE_Y'
                                        */
  real_T g2ms2_Wilcoxon_Gain_h;        /* Expression: 9.80655
                                        * Referenced by: '<S40>/g2ms^2_Wilcoxon'
                                        */
  real_T GainADC_ACCE_Z_Gain;          /* Expression: 10
                                        * Referenced by: '<S5>/GainADC_ACCE_Z'
                                        */
  real_T g2ms2_Wilcoxon_Gain_k;        /* Expression: 9.80655
                                        * Referenced by: '<S41>/g2ms^2_Wilcoxon'
                                        */
  real_T GainADC_GIRO_X_Gain;          /* Expression: 10
                                        * Referenced by: '<S5>/GainADC_GIRO_X'
                                        */
  real_T GainADC_GIRO_Y_Gain;          /* Expression: 10
                                        * Referenced by: '<S5>/GainADC_GIRO_Y'
                                        */
  real_T GainADC_GIRO_Z_Gain;          /* Expression: 10
                                        * Referenced by: '<S5>/GainADC_GIRO_Z'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Pierre_Controller_LTI_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_Pierre_Controller_LTI_T Pierre_Controller_LTI_P;

/* Block signals (default storage) */
extern B_Pierre_Controller_LTI_T Pierre_Controller_LTI_B;

/* Block states (default storage) */
extern DW_Pierre_Controller_LTI_T Pierre_Controller_LTI_DW;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Pierre_Controller_LTI_T Pierre_Controller_LTI_Y;

/* Model entry point functions */
extern void Pierre_Controller_LTI_initialize(void);
extern void Pierre_Controller_LTI_output(void);
extern void Pierre_Controller_LTI_update(void);
extern void Pierre_Controller_LTI_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Pierre_Controller_LT_T *const Pierre_Controller_LTI_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Pierre_Controller_LTI'
 * '<S1>'   : 'Pierre_Controller_LTI/Clock'
 * '<S2>'   : 'Pierre_Controller_LTI/Controller'
 * '<S3>'   : 'Pierre_Controller_LTI/Current monitoring'
 * '<S4>'   : 'Pierre_Controller_LTI/Force_measure'
 * '<S5>'   : 'Pierre_Controller_LTI/IMU'
 * '<S6>'   : 'Pierre_Controller_LTI/MATLAB Function'
 * '<S7>'   : 'Pierre_Controller_LTI/Micromega PMA1'
 * '<S8>'   : 'Pierre_Controller_LTI/Micromega PMA2'
 * '<S9>'   : 'Pierre_Controller_LTI/Micromega PMA3'
 * '<S10>'  : 'Pierre_Controller_LTI/Micromega PMA4'
 * '<S11>'  : 'Pierre_Controller_LTI/Micromega PMA5'
 * '<S12>'  : 'Pierre_Controller_LTI/Micromega PMA6'
 * '<S13>'  : 'Pierre_Controller_LTI/PCB208C01'
 * '<S14>'  : 'Pierre_Controller_LTI/Payload Accelerometers'
 * '<S15>'  : 'Pierre_Controller_LTI/RTI Data'
 * '<S16>'  : 'Pierre_Controller_LTI/Shaker TMS_2007E'
 * '<S17>'  : 'Pierre_Controller_LTI/SineSweep_TMS'
 * '<S18>'  : 'Pierre_Controller_LTI/Sinus_TMS'
 * '<S19>'  : 'Pierre_Controller_LTI/TMS_Command'
 * '<S20>'  : 'Pierre_Controller_LTI/Current monitoring/PMA1_monitoring'
 * '<S21>'  : 'Pierre_Controller_LTI/Current monitoring/PMA2_monitoring'
 * '<S22>'  : 'Pierre_Controller_LTI/Current monitoring/PMA3_monitoring'
 * '<S23>'  : 'Pierre_Controller_LTI/Current monitoring/PMA4_monitoring'
 * '<S24>'  : 'Pierre_Controller_LTI/Current monitoring/PMA5_monitoring'
 * '<S25>'  : 'Pierre_Controller_LTI/Current monitoring/PMA6_monitoring'
 * '<S26>'  : 'Pierre_Controller_LTI/Current monitoring/Subsystem'
 * '<S27>'  : 'Pierre_Controller_LTI/Current monitoring/Subsystem1'
 * '<S28>'  : 'Pierre_Controller_LTI/Current monitoring/Subsystem2'
 * '<S29>'  : 'Pierre_Controller_LTI/Current monitoring/Subsystem3'
 * '<S30>'  : 'Pierre_Controller_LTI/Current monitoring/Subsystem4'
 * '<S31>'  : 'Pierre_Controller_LTI/Current monitoring/Subsystem5'
 * '<S32>'  : 'Pierre_Controller_LTI/IMU/ACCE_X'
 * '<S33>'  : 'Pierre_Controller_LTI/IMU/ACCE_Y'
 * '<S34>'  : 'Pierre_Controller_LTI/IMU/ACCE_Z'
 * '<S35>'  : 'Pierre_Controller_LTI/IMU/GIRO_X'
 * '<S36>'  : 'Pierre_Controller_LTI/IMU/GIRO_Y'
 * '<S37>'  : 'Pierre_Controller_LTI/IMU/GIRO_Z'
 * '<S38>'  : 'Pierre_Controller_LTI/IMU/Subsystem1'
 * '<S39>'  : 'Pierre_Controller_LTI/IMU/Subsystem12'
 * '<S40>'  : 'Pierre_Controller_LTI/IMU/Subsystem13'
 * '<S41>'  : 'Pierre_Controller_LTI/IMU/Subsystem14'
 * '<S42>'  : 'Pierre_Controller_LTI/IMU/Subsystem2'
 * '<S43>'  : 'Pierre_Controller_LTI/IMU/Subsystem3'
 * '<S44>'  : 'Pierre_Controller_LTI/Payload Accelerometers/Accelerometer_Wilcoxon#4546'
 * '<S45>'  : 'Pierre_Controller_LTI/Payload Accelerometers/Accelerometer_Wilcoxon#4786'
 * '<S46>'  : 'Pierre_Controller_LTI/Payload Accelerometers/Accelerometer_Wilcoxon#4787'
 * '<S47>'  : 'Pierre_Controller_LTI/Payload Accelerometers/DocBlock'
 * '<S48>'  : 'Pierre_Controller_LTI/Payload Accelerometers/Subsystem12'
 * '<S49>'  : 'Pierre_Controller_LTI/Payload Accelerometers/Subsystem13'
 * '<S50>'  : 'Pierre_Controller_LTI/Payload Accelerometers/Subsystem14'
 * '<S51>'  : 'Pierre_Controller_LTI/RTI Data/RTI Data Store'
 * '<S52>'  : 'Pierre_Controller_LTI/RTI Data/RTI Data Store/RTI Data Store'
 * '<S53>'  : 'Pierre_Controller_LTI/RTI Data/RTI Data Store/RTI Data Store/RTI Data Store'
 * '<S54>'  : 'Pierre_Controller_LTI/Sinus_TMS/MATLAB Function'
 */
#endif                                 /* Pierre_Controller_LTI_h_ */
