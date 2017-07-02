/* Include files */

#include <stddef.h>
#include "blas.h"
#include "lab5_template_sivom_sfun.h"
#include "c3_lab5_template_sivom.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "lab5_template_sivom_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c3_debug_family_names[64] = { "x_E", "eta", "v_B", "w_B",
  "phi", "theta", "psi", "p", "q", "r", "rho", "R", "r1", "r2", "r3", "r4", "Ct",
  "Cq", "m", "g", "I_B", "eb", "A", "f1", "f2", "f3", "f4", "f_B", "Qr1", "Qr2",
  "Qr3", "Qr4", "m1", "m2", "m3", "m4", "M_B", "R1", "R2", "R3", "R_be", "R_eb",
  "x_E_dot", "T", "eta_dot", "g1", "v_E_dot", "v_B_dot", "i1", "i2", "i3", "i4",
  "Inew", "temp2", "temp3", "w_B_dot", "nargin", "nargout", "w1", "w2", "w3",
  "w4", "states", "output" };

/* Function Declarations */
static void initialize_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void initialize_params_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void enable_c3_lab5_template_sivom(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance);
static void disable_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void set_sim_state_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance, const mxArray *c3_st);
static void finalize_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void sf_c3_lab5_template_sivom(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance);
static void c3_chartstep_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void initSimStructsc3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_output, const char_T *c3_identifier, real_T
  c3_y[12]);
static void c3_b_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[12]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_c_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_d_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[3]);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_e_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[9]);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_emlrt_marshallOut(char * c3_u);
static const mxArray *c3_b_emlrt_marshallOut(uint32_T c3_u);
static void c3_b_info_helper(const mxArray **c3_info);
static real_T c3_mpower(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
  real_T c3_a);
static void c3_eml_scalar_eg(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance);
static real_T c3_sign(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T c3_x);
static void c3_cross(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T c3_a[3], real_T c3_b[3], real_T c3_c[3]);
static void c3_b_eml_scalar_eg(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance);
static void c3_c_eml_scalar_eg(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance);
static real_T c3_sec(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T c3_x);
static real_T c3_eml_div(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
  real_T c3_x, real_T c3_y);
static void c3_inv(SFc3_lab5_template_sivomInstanceStruct *chartInstance, real_T
                   c3_x[9], real_T c3_y[9]);
static void c3_inv3x3(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T c3_x[9], real_T c3_y[9]);
static real_T c3_abs(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T c3_x);
static real_T c3_norm(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T c3_x[9]);
static void c3_eml_warning(SFc3_lab5_template_sivomInstanceStruct *chartInstance);
static void c3_b_eml_warning(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, char_T c3_varargin_2[14]);
static void c3_f_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_sprintf, const char_T *c3_identifier, char_T
  c3_y[14]);
static void c3_g_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  char_T c3_y[14]);
static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_h_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_i_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_b_is_active_c3_lab5_template_sivom, const
  char_T *c3_identifier);
static uint8_T c3_j_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sign(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T *c3_x);
static void c3_b_sec(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T *c3_x);
static void init_dsm_address_info(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c3_is_active_c3_lab5_template_sivom = 0U;
}

static void initialize_params_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
}

static void enable_c3_lab5_template_sivom(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c3_update_debugger_state_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[12];
  const mxArray *c3_b_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_b_u;
  const mxArray *c3_c_y = NULL;
  real_T (*c3_output)[12];
  c3_output = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellarray(2), FALSE);
  for (c3_i0 = 0; c3_i0 < 12; c3_i0++) {
    c3_u[c3_i0] = (*c3_output)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 12), FALSE);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  c3_hoistedGlobal = chartInstance->c3_is_active_c3_lab5_template_sivom;
  c3_b_u = c3_hoistedGlobal;
  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", &c3_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  sf_mex_assign(&c3_st, c3_y, FALSE);
  return c3_st;
}

static void set_sim_state_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance, const mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[12];
  int32_T c3_i1;
  real_T (*c3_output)[12];
  c3_output = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c3_doneDoubleBufferReInit = TRUE;
  c3_u = sf_mex_dup(c3_st);
  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)),
                      "output", c3_dv0);
  for (c3_i1 = 0; c3_i1 < 12; c3_i1++) {
    (*c3_output)[c3_i1] = c3_dv0[c3_i1];
  }

  chartInstance->c3_is_active_c3_lab5_template_sivom = c3_i_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 1)),
     "is_active_c3_lab5_template_sivom");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_lab5_template_sivom(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
}

static void sf_c3_lab5_template_sivom(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance)
{
  int32_T c3_i2;
  int32_T c3_i3;
  real_T *c3_w1;
  real_T *c3_w2;
  real_T *c3_w3;
  real_T *c3_w4;
  real_T (*c3_states)[12];
  real_T (*c3_output)[12];
  c3_states = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 4);
  c3_w4 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c3_w3 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c3_w2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c3_output = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_w1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c3_w1, 0U);
  for (c3_i2 = 0; c3_i2 < 12; c3_i2++) {
    _SFD_DATA_RANGE_CHECK((*c3_output)[c3_i2], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c3_w2, 2U);
  _SFD_DATA_RANGE_CHECK(*c3_w3, 3U);
  _SFD_DATA_RANGE_CHECK(*c3_w4, 4U);
  for (c3_i3 = 0; c3_i3 < 12; c3_i3++) {
    _SFD_DATA_RANGE_CHECK((*c3_states)[c3_i3], 5U);
  }

  chartInstance->c3_sfEvent = CALL_EVENT;
  c3_chartstep_c3_lab5_template_sivom(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_lab5_template_sivomMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c3_chartstep_c3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_c_hoistedGlobal;
  real_T c3_d_hoistedGlobal;
  real_T c3_w1;
  real_T c3_w2;
  real_T c3_w3;
  real_T c3_w4;
  int32_T c3_i4;
  real_T c3_states[12];
  uint32_T c3_debug_family_var_map[64];
  real_T c3_x_E[3];
  real_T c3_eta[3];
  real_T c3_v_B[3];
  real_T c3_w_B[3];
  real_T c3_phi;
  real_T c3_theta;
  real_T c3_psi;
  real_T c3_p;
  real_T c3_q;
  real_T c3_r;
  real_T c3_rho;
  real_T c3_R;
  real_T c3_r1[3];
  real_T c3_r2[3];
  real_T c3_r3[3];
  real_T c3_r4[3];
  real_T c3_Ct;
  real_T c3_Cq;
  real_T c3_m;
  real_T c3_g;
  real_T c3_I_B[9];
  real_T c3_eb[3];
  real_T c3_A;
  real_T c3_f1[3];
  real_T c3_f2[3];
  real_T c3_f3[3];
  real_T c3_f4[3];
  real_T c3_f_B[3];
  real_T c3_Qr1[3];
  real_T c3_Qr2[3];
  real_T c3_Qr3[3];
  real_T c3_Qr4[3];
  real_T c3_m1[3];
  real_T c3_m2[3];
  real_T c3_m3[3];
  real_T c3_m4[3];
  real_T c3_M_B[3];
  real_T c3_R1[9];
  real_T c3_R2[9];
  real_T c3_R3[9];
  real_T c3_R_be[9];
  real_T c3_R_eb[9];
  real_T c3_x_E_dot[3];
  real_T c3_T[9];
  real_T c3_eta_dot[3];
  real_T c3_g1[3];
  real_T c3_v_E_dot[3];
  real_T c3_v_B_dot[3];
  real_T c3_i1[3];
  real_T c3_i2[3];
  real_T c3_i3[3];
  real_T c3_b_i4[3];
  real_T c3_Inew[3];
  real_T c3_temp2[3];
  real_T c3_temp3[3];
  real_T c3_w_B_dot[3];
  real_T c3_nargin = 5.0;
  real_T c3_nargout = 1.0;
  real_T c3_output[12];
  int32_T c3_i5;
  int32_T c3_i6;
  int32_T c3_i7;
  int32_T c3_i8;
  int32_T c3_i9;
  static real_T c3_dv1[3] = { 0.3, 0.0, 0.0 };

  int32_T c3_i10;
  static real_T c3_dv2[3] = { 0.0, 0.3, 0.0 };

  int32_T c3_i11;
  static real_T c3_dv3[3] = { -0.3, 0.0, 0.0 };

  int32_T c3_i12;
  static real_T c3_dv4[3] = { 0.0, -0.3, 0.0 };

  int32_T c3_i13;
  static real_T c3_a[9] = { 0.06, 0.0, 0.0, 0.0, 0.06, 0.0, 0.0, 0.0, 0.12 };

  int32_T c3_i14;
  static real_T c3_b[3] = { 0.0, 0.0, -1.0 };

  real_T c3_b_a;
  real_T c3_y;
  real_T c3_b_b;
  real_T c3_b_y;
  real_T c3_c_a;
  real_T c3_c_y;
  real_T c3_d_a;
  int32_T c3_i15;
  real_T c3_e_a;
  real_T c3_d_y;
  real_T c3_c_b;
  real_T c3_e_y;
  real_T c3_f_a;
  real_T c3_f_y;
  real_T c3_g_a;
  int32_T c3_i16;
  real_T c3_h_a;
  real_T c3_g_y;
  real_T c3_d_b;
  real_T c3_h_y;
  real_T c3_i_a;
  real_T c3_i_y;
  real_T c3_j_a;
  int32_T c3_i17;
  real_T c3_k_a;
  real_T c3_j_y;
  real_T c3_e_b;
  real_T c3_k_y;
  real_T c3_l_a;
  real_T c3_l_y;
  real_T c3_m_a;
  int32_T c3_i18;
  int32_T c3_i19;
  real_T c3_n_a;
  real_T c3_m_y;
  real_T c3_f_b;
  real_T c3_n_y;
  real_T c3_o_a;
  real_T c3_o_y;
  real_T c3_p_a;
  real_T c3_p_y;
  real_T c3_q_a;
  real_T c3_g_b;
  real_T c3_q_y;
  int32_T c3_i20;
  real_T c3_r_a;
  real_T c3_r_y;
  real_T c3_h_b;
  real_T c3_s_y;
  real_T c3_s_a;
  real_T c3_t_y;
  real_T c3_t_a;
  real_T c3_u_y;
  real_T c3_u_a;
  real_T c3_i_b;
  real_T c3_v_y;
  int32_T c3_i21;
  real_T c3_v_a;
  real_T c3_w_y;
  real_T c3_j_b;
  real_T c3_x_y;
  real_T c3_w_a;
  real_T c3_y_y;
  real_T c3_x_a;
  real_T c3_ab_y;
  real_T c3_y_a;
  real_T c3_k_b;
  real_T c3_bb_y;
  int32_T c3_i22;
  real_T c3_ab_a;
  real_T c3_cb_y;
  real_T c3_l_b;
  real_T c3_db_y;
  real_T c3_bb_a;
  real_T c3_eb_y;
  real_T c3_cb_a;
  real_T c3_fb_y;
  real_T c3_db_a;
  real_T c3_m_b;
  real_T c3_gb_y;
  int32_T c3_i23;
  int32_T c3_i24;
  real_T c3_dv5[3];
  int32_T c3_i25;
  real_T c3_b_f1[3];
  real_T c3_C[3];
  int32_T c3_i26;
  int32_T c3_i27;
  real_T c3_dv6[3];
  int32_T c3_i28;
  real_T c3_b_f2[3];
  int32_T c3_i29;
  int32_T c3_i30;
  real_T c3_dv7[3];
  int32_T c3_i31;
  real_T c3_b_f3[3];
  int32_T c3_i32;
  int32_T c3_i33;
  real_T c3_dv8[3];
  int32_T c3_i34;
  real_T c3_b_f4[3];
  int32_T c3_i35;
  int32_T c3_i36;
  real_T c3_x;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_g_x;
  real_T c3_h_x;
  int32_T c3_i37;
  int32_T c3_i38;
  static real_T c3_dv9[3] = { 1.0, 0.0, 0.0 };

  real_T c3_i_x;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_p_x;
  int32_T c3_i39;
  int32_T c3_i40;
  static real_T c3_dv10[3] = { 0.0, 1.0, 0.0 };

  real_T c3_q_x;
  real_T c3_r_x;
  real_T c3_s_x;
  real_T c3_t_x;
  real_T c3_u_x;
  real_T c3_v_x;
  real_T c3_w_x;
  real_T c3_x_x;
  int32_T c3_i41;
  int32_T c3_i42;
  static real_T c3_dv11[3] = { 0.0, 0.0, 1.0 };

  int32_T c3_i43;
  real_T c3_eb_a[9];
  int32_T c3_i44;
  real_T c3_n_b[9];
  int32_T c3_i45;
  int32_T c3_i46;
  int32_T c3_i47;
  real_T c3_hb_y[9];
  int32_T c3_i48;
  int32_T c3_i49;
  int32_T c3_i50;
  int32_T c3_i51;
  int32_T c3_i52;
  int32_T c3_i53;
  int32_T c3_i54;
  int32_T c3_i55;
  int32_T c3_i56;
  int32_T c3_i57;
  int32_T c3_i58;
  int32_T c3_i59;
  int32_T c3_i60;
  int32_T c3_i61;
  int32_T c3_i62;
  int32_T c3_i63;
  int32_T c3_i64;
  int32_T c3_i65;
  int32_T c3_i66;
  int32_T c3_i67;
  real_T c3_o_b[3];
  int32_T c3_i68;
  int32_T c3_i69;
  int32_T c3_i70;
  int32_T c3_i71;
  int32_T c3_i72;
  int32_T c3_i73;
  int32_T c3_i74;
  int32_T c3_i75;
  int32_T c3_i76;
  real_T c3_y_x;
  real_T c3_ab_x;
  real_T c3_bb_x;
  real_T c3_cb_x;
  real_T c3_fb_a;
  real_T c3_p_b;
  real_T c3_ib_y;
  real_T c3_db_x;
  real_T c3_eb_x;
  real_T c3_fb_x;
  real_T c3_gb_x;
  real_T c3_gb_a;
  real_T c3_q_b;
  real_T c3_jb_y;
  real_T c3_hb_x;
  real_T c3_ib_x;
  real_T c3_jb_x;
  real_T c3_kb_x;
  real_T c3_lb_x;
  real_T c3_mb_x;
  real_T c3_hb_a;
  real_T c3_r_b;
  real_T c3_kb_y;
  real_T c3_nb_x;
  real_T c3_ob_x;
  real_T c3_ib_a;
  real_T c3_s_b;
  real_T c3_lb_y;
  int32_T c3_i77;
  int32_T c3_i78;
  int32_T c3_i79;
  int32_T c3_i80;
  int32_T c3_i81;
  int32_T c3_i82;
  int32_T c3_i83;
  int32_T c3_i84;
  int32_T c3_i85;
  int32_T c3_i86;
  int32_T c3_i87;
  int32_T c3_i88;
  static real_T c3_dv12[3] = { 0.0, 0.0, 9.8 };

  int32_T c3_i89;
  int32_T c3_i90;
  int32_T c3_i91;
  int32_T c3_i92;
  real_T c3_mb_y[3];
  int32_T c3_i93;
  int32_T c3_i94;
  int32_T c3_i95;
  int32_T c3_i96;
  int32_T c3_i97;
  int32_T c3_i98;
  int32_T c3_i99;
  int32_T c3_i100;
  int32_T c3_i101;
  int32_T c3_i102;
  int32_T c3_i103;
  int32_T c3_i104;
  int32_T c3_i105;
  int32_T c3_i106;
  real_T c3_t_b;
  real_T c3_nb_y;
  real_T c3_jb_a;
  real_T c3_u_b;
  real_T c3_ob_y;
  real_T c3_v_b;
  real_T c3_pb_y;
  real_T c3_kb_a;
  real_T c3_w_b;
  real_T c3_qb_y;
  real_T c3_x_b;
  real_T c3_rb_y;
  real_T c3_lb_a;
  real_T c3_y_b;
  real_T c3_sb_y;
  real_T c3_ab_b;
  real_T c3_tb_y;
  real_T c3_mb_a;
  real_T c3_bb_b;
  real_T c3_ub_y;
  real_T c3_cb_b;
  real_T c3_vb_y;
  real_T c3_nb_a;
  real_T c3_db_b;
  real_T c3_wb_y;
  real_T c3_eb_b;
  real_T c3_xb_y;
  real_T c3_ob_a;
  real_T c3_fb_b;
  real_T c3_yb_y;
  real_T c3_gb_b;
  real_T c3_ac_y;
  real_T c3_pb_a;
  real_T c3_hb_b;
  real_T c3_bc_y;
  real_T c3_ib_b;
  real_T c3_cc_y;
  real_T c3_qb_a;
  real_T c3_jb_b;
  real_T c3_dc_y;
  int32_T c3_i107;
  int32_T c3_i108;
  int32_T c3_i109;
  int32_T c3_i110;
  int32_T c3_i111;
  int32_T c3_i112;
  int32_T c3_i113;
  int32_T c3_i114;
  int32_T c3_i115;
  int32_T c3_i116;
  int32_T c3_i117;
  int32_T c3_i118;
  real_T c3_b_w_B[3];
  int32_T c3_i119;
  real_T c3_b_temp2[3];
  real_T c3_dv13[3];
  int32_T c3_i120;
  int32_T c3_i121;
  real_T c3_rb_a[9];
  int32_T c3_i122;
  int32_T c3_i123;
  int32_T c3_i124;
  int32_T c3_i125;
  int32_T c3_i126;
  int32_T c3_i127;
  int32_T c3_i128;
  int32_T c3_i129;
  int32_T c3_i130;
  int32_T c3_i131;
  int32_T c3_i132;
  real_T *c3_b_w1;
  real_T *c3_b_w2;
  real_T *c3_b_w3;
  real_T *c3_b_w4;
  real_T (*c3_b_output)[12];
  real_T (*c3_b_states)[12];
  c3_b_states = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 4);
  c3_b_w4 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c3_b_w3 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c3_b_w2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c3_b_output = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_b_w1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
  c3_hoistedGlobal = *c3_b_w1;
  c3_b_hoistedGlobal = *c3_b_w2;
  c3_c_hoistedGlobal = *c3_b_w3;
  c3_d_hoistedGlobal = *c3_b_w4;
  c3_w1 = c3_hoistedGlobal;
  c3_w2 = c3_b_hoistedGlobal;
  c3_w3 = c3_c_hoistedGlobal;
  c3_w4 = c3_d_hoistedGlobal;
  for (c3_i4 = 0; c3_i4 < 12; c3_i4++) {
    c3_states[c3_i4] = (*c3_b_states)[c3_i4];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 64U, 64U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_x_E, 0U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_eta, 1U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_v_B, 2U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_w_B, 3U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_phi, 4U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta, 5U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_psi, 6U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p, 7U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q, 8U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r, 9U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_rho, 10U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_R, 11U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_r1, 12U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_r2, 13U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_r3, 14U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_r4, 15U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Ct, 16U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Cq, 17U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_m, 18U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_g, 19U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_I_B, 20U, c3_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_eb, 21U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_A, 22U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_f1, 23U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_f2, 24U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_f3, 25U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_f4, 26U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_f_B, 27U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Qr1, 28U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Qr2, 29U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Qr3, 30U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Qr4, 31U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_m1, 32U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_m2, 33U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_m3, 34U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_m4, 35U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_M_B, 36U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R1, 37U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R2, 38U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R3, 39U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R_be, 40U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_R_eb, 41U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_x_E_dot, 42U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_T, 43U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_eta_dot, 44U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_g1, 45U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_v_E_dot, 46U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_v_B_dot, 47U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_i1, 48U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_i2, 49U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_i3, 50U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_b_i4, 51U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Inew, 52U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_temp2, 53U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_temp3, 54U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_w_B_dot, 55U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 56U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 57U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_w1, 58U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_w2, 59U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_w3, 60U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_w4, 61U, c3_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_states, 62U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_output, 63U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 9);
  for (c3_i5 = 0; c3_i5 < 3; c3_i5++) {
    c3_x_E[c3_i5] = c3_states[c3_i5];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 10);
  for (c3_i6 = 0; c3_i6 < 3; c3_i6++) {
    c3_eta[c3_i6] = c3_states[c3_i6 + 3];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 11);
  for (c3_i7 = 0; c3_i7 < 3; c3_i7++) {
    c3_v_B[c3_i7] = c3_states[c3_i7 + 6];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 12);
  for (c3_i8 = 0; c3_i8 < 3; c3_i8++) {
    c3_w_B[c3_i8] = c3_states[c3_i8 + 9];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 14);
  c3_phi = c3_eta[0];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 15);
  c3_theta = c3_eta[1];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 16);
  c3_psi = c3_eta[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 19);
  c3_p = c3_w_B[0];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 20);
  c3_q = c3_w_B[1];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 21);
  c3_r = c3_w_B[2];
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 24);
  c3_rho = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 25);
  c3_R = 0.15;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 26);
  for (c3_i9 = 0; c3_i9 < 3; c3_i9++) {
    c3_r1[c3_i9] = c3_dv1[c3_i9];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 27);
  for (c3_i10 = 0; c3_i10 < 3; c3_i10++) {
    c3_r2[c3_i10] = c3_dv2[c3_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 28);
  for (c3_i11 = 0; c3_i11 < 3; c3_i11++) {
    c3_r3[c3_i11] = c3_dv3[c3_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 29);
  for (c3_i12 = 0; c3_i12 < 3; c3_i12++) {
    c3_r4[c3_i12] = c3_dv4[c3_i12];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 30);
  c3_Ct = 0.012;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 31);
  c3_Cq = 0.0026;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 32);
  c3_m = 2.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 33);
  c3_g = 9.8;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 34);
  for (c3_i13 = 0; c3_i13 < 9; c3_i13++) {
    c3_I_B[c3_i13] = c3_a[c3_i13];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 35);
  for (c3_i14 = 0; c3_i14 < 3; c3_i14++) {
    c3_eb[c3_i14] = c3_b[c3_i14];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 36);
  c3_A = 0.070685834705770348;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 37);
  c3_b_a = c3_w1;
  c3_y = c3_b_a * 0.15;
  c3_b_b = c3_mpower(chartInstance, c3_y);
  c3_b_y = 0.070685834705770348 * c3_b_b;
  c3_c_a = c3_b_y;
  c3_c_y = c3_c_a * 0.012;
  c3_d_a = c3_c_y;
  for (c3_i15 = 0; c3_i15 < 3; c3_i15++) {
    c3_f1[c3_i15] = c3_d_a * c3_b[c3_i15];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 38);
  c3_e_a = c3_w2;
  c3_d_y = c3_e_a * 0.15;
  c3_c_b = c3_mpower(chartInstance, c3_d_y);
  c3_e_y = 0.070685834705770348 * c3_c_b;
  c3_f_a = c3_e_y;
  c3_f_y = c3_f_a * 0.012;
  c3_g_a = c3_f_y;
  for (c3_i16 = 0; c3_i16 < 3; c3_i16++) {
    c3_f2[c3_i16] = c3_g_a * c3_b[c3_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 39);
  c3_h_a = c3_w3;
  c3_g_y = c3_h_a * 0.15;
  c3_d_b = c3_mpower(chartInstance, c3_g_y);
  c3_h_y = 0.070685834705770348 * c3_d_b;
  c3_i_a = c3_h_y;
  c3_i_y = c3_i_a * 0.012;
  c3_j_a = c3_i_y;
  for (c3_i17 = 0; c3_i17 < 3; c3_i17++) {
    c3_f3[c3_i17] = c3_j_a * c3_b[c3_i17];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 40);
  c3_k_a = c3_w4;
  c3_j_y = c3_k_a * 0.15;
  c3_e_b = c3_mpower(chartInstance, c3_j_y);
  c3_k_y = 0.070685834705770348 * c3_e_b;
  c3_l_a = c3_k_y;
  c3_l_y = c3_l_a * 0.012;
  c3_m_a = c3_l_y;
  for (c3_i18 = 0; c3_i18 < 3; c3_i18++) {
    c3_f4[c3_i18] = c3_m_a * c3_b[c3_i18];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 41);
  for (c3_i19 = 0; c3_i19 < 3; c3_i19++) {
    c3_f_B[c3_i19] = ((c3_f1[c3_i19] + c3_f2[c3_i19]) + c3_f3[c3_i19]) +
      c3_f4[c3_i19];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 42);
  c3_n_a = c3_w1;
  c3_m_y = c3_n_a * 0.15;
  c3_f_b = c3_mpower(chartInstance, c3_m_y);
  c3_n_y = -0.070685834705770348 * c3_f_b;
  c3_o_a = c3_n_y;
  c3_o_y = c3_o_a * 0.15;
  c3_p_a = c3_o_y;
  c3_p_y = c3_p_a * 0.0026;
  c3_q_a = c3_p_y;
  c3_g_b = c3_w1;
  c3_b_sign(chartInstance, &c3_g_b);
  c3_q_y = c3_q_a * c3_g_b;
  for (c3_i20 = 0; c3_i20 < 3; c3_i20++) {
    c3_Qr1[c3_i20] = c3_q_y * c3_eb[c3_i20];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 43);
  c3_r_a = c3_w2;
  c3_r_y = c3_r_a * 0.15;
  c3_h_b = c3_mpower(chartInstance, c3_r_y);
  c3_s_y = -0.070685834705770348 * c3_h_b;
  c3_s_a = c3_s_y;
  c3_t_y = c3_s_a * 0.15;
  c3_t_a = c3_t_y;
  c3_u_y = c3_t_a * 0.0026;
  c3_u_a = c3_u_y;
  c3_i_b = c3_w2;
  c3_b_sign(chartInstance, &c3_i_b);
  c3_v_y = c3_u_a * c3_i_b;
  for (c3_i21 = 0; c3_i21 < 3; c3_i21++) {
    c3_Qr2[c3_i21] = c3_v_y * c3_eb[c3_i21];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 44);
  c3_v_a = c3_w3;
  c3_w_y = c3_v_a * 0.15;
  c3_j_b = c3_mpower(chartInstance, c3_w_y);
  c3_x_y = -0.070685834705770348 * c3_j_b;
  c3_w_a = c3_x_y;
  c3_y_y = c3_w_a * 0.15;
  c3_x_a = c3_y_y;
  c3_ab_y = c3_x_a * 0.0026;
  c3_y_a = c3_ab_y;
  c3_k_b = c3_w3;
  c3_b_sign(chartInstance, &c3_k_b);
  c3_bb_y = c3_y_a * c3_k_b;
  for (c3_i22 = 0; c3_i22 < 3; c3_i22++) {
    c3_Qr3[c3_i22] = c3_bb_y * c3_eb[c3_i22];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 45);
  c3_ab_a = c3_w4;
  c3_cb_y = c3_ab_a * 0.15;
  c3_l_b = c3_mpower(chartInstance, c3_cb_y);
  c3_db_y = -0.070685834705770348 * c3_l_b;
  c3_bb_a = c3_db_y;
  c3_eb_y = c3_bb_a * 0.15;
  c3_cb_a = c3_eb_y;
  c3_fb_y = c3_cb_a * 0.0026;
  c3_db_a = c3_fb_y;
  c3_m_b = c3_w4;
  c3_b_sign(chartInstance, &c3_m_b);
  c3_gb_y = c3_db_a * c3_m_b;
  for (c3_i23 = 0; c3_i23 < 3; c3_i23++) {
    c3_Qr4[c3_i23] = c3_gb_y * c3_eb[c3_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 47);
  for (c3_i24 = 0; c3_i24 < 3; c3_i24++) {
    c3_dv5[c3_i24] = c3_dv1[c3_i24];
  }

  for (c3_i25 = 0; c3_i25 < 3; c3_i25++) {
    c3_b_f1[c3_i25] = c3_f1[c3_i25];
  }

  c3_cross(chartInstance, c3_dv5, c3_b_f1, c3_C);
  for (c3_i26 = 0; c3_i26 < 3; c3_i26++) {
    c3_m1[c3_i26] = c3_Qr1[c3_i26] + c3_C[c3_i26];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 48);
  for (c3_i27 = 0; c3_i27 < 3; c3_i27++) {
    c3_dv6[c3_i27] = c3_dv2[c3_i27];
  }

  for (c3_i28 = 0; c3_i28 < 3; c3_i28++) {
    c3_b_f2[c3_i28] = c3_f2[c3_i28];
  }

  c3_cross(chartInstance, c3_dv6, c3_b_f2, c3_C);
  for (c3_i29 = 0; c3_i29 < 3; c3_i29++) {
    c3_m2[c3_i29] = c3_Qr2[c3_i29] + c3_C[c3_i29];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 49);
  for (c3_i30 = 0; c3_i30 < 3; c3_i30++) {
    c3_dv7[c3_i30] = c3_dv3[c3_i30];
  }

  for (c3_i31 = 0; c3_i31 < 3; c3_i31++) {
    c3_b_f3[c3_i31] = c3_f3[c3_i31];
  }

  c3_cross(chartInstance, c3_dv7, c3_b_f3, c3_C);
  for (c3_i32 = 0; c3_i32 < 3; c3_i32++) {
    c3_m3[c3_i32] = c3_Qr3[c3_i32] + c3_C[c3_i32];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 50);
  for (c3_i33 = 0; c3_i33 < 3; c3_i33++) {
    c3_dv8[c3_i33] = c3_dv4[c3_i33];
  }

  for (c3_i34 = 0; c3_i34 < 3; c3_i34++) {
    c3_b_f4[c3_i34] = c3_f4[c3_i34];
  }

  c3_cross(chartInstance, c3_dv8, c3_b_f4, c3_C);
  for (c3_i35 = 0; c3_i35 < 3; c3_i35++) {
    c3_m4[c3_i35] = c3_Qr4[c3_i35] + c3_C[c3_i35];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 52);
  for (c3_i36 = 0; c3_i36 < 3; c3_i36++) {
    c3_M_B[c3_i36] = ((c3_m1[c3_i36] + c3_m2[c3_i36]) + c3_m3[c3_i36]) +
      c3_m4[c3_i36];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 56);
  c3_x = c3_phi;
  c3_b_x = c3_x;
  c3_b_x = muDoubleScalarCos(c3_b_x);
  c3_c_x = c3_phi;
  c3_d_x = c3_c_x;
  c3_d_x = muDoubleScalarSin(c3_d_x);
  c3_e_x = c3_phi;
  c3_f_x = c3_e_x;
  c3_f_x = muDoubleScalarSin(c3_f_x);
  c3_g_x = c3_phi;
  c3_h_x = c3_g_x;
  c3_h_x = muDoubleScalarCos(c3_h_x);
  c3_i37 = 0;
  for (c3_i38 = 0; c3_i38 < 3; c3_i38++) {
    c3_R1[c3_i37] = c3_dv9[c3_i38];
    c3_i37 += 3;
  }

  c3_R1[1] = 0.0;
  c3_R1[4] = c3_b_x;
  c3_R1[7] = c3_d_x;
  c3_R1[2] = 0.0;
  c3_R1[5] = -c3_f_x;
  c3_R1[8] = c3_h_x;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 57);
  c3_i_x = c3_theta;
  c3_j_x = c3_i_x;
  c3_j_x = muDoubleScalarCos(c3_j_x);
  c3_k_x = c3_theta;
  c3_l_x = c3_k_x;
  c3_l_x = muDoubleScalarSin(c3_l_x);
  c3_m_x = c3_theta;
  c3_n_x = c3_m_x;
  c3_n_x = muDoubleScalarSin(c3_n_x);
  c3_o_x = c3_theta;
  c3_p_x = c3_o_x;
  c3_p_x = muDoubleScalarCos(c3_p_x);
  c3_R2[0] = c3_j_x;
  c3_R2[3] = 0.0;
  c3_R2[6] = -c3_l_x;
  c3_i39 = 0;
  for (c3_i40 = 0; c3_i40 < 3; c3_i40++) {
    c3_R2[c3_i39 + 1] = c3_dv10[c3_i40];
    c3_i39 += 3;
  }

  c3_R2[2] = c3_n_x;
  c3_R2[5] = 0.0;
  c3_R2[8] = c3_p_x;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 58);
  c3_q_x = c3_psi;
  c3_r_x = c3_q_x;
  c3_r_x = muDoubleScalarCos(c3_r_x);
  c3_s_x = c3_psi;
  c3_t_x = c3_s_x;
  c3_t_x = muDoubleScalarSin(c3_t_x);
  c3_u_x = c3_psi;
  c3_v_x = c3_u_x;
  c3_v_x = muDoubleScalarSin(c3_v_x);
  c3_w_x = c3_psi;
  c3_x_x = c3_w_x;
  c3_x_x = muDoubleScalarCos(c3_x_x);
  c3_R3[0] = c3_r_x;
  c3_R3[3] = c3_t_x;
  c3_R3[6] = 0.0;
  c3_R3[1] = -c3_v_x;
  c3_R3[4] = c3_x_x;
  c3_R3[7] = 0.0;
  c3_i41 = 0;
  for (c3_i42 = 0; c3_i42 < 3; c3_i42++) {
    c3_R3[c3_i41 + 2] = c3_dv11[c3_i42];
    c3_i41 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 59);
  for (c3_i43 = 0; c3_i43 < 9; c3_i43++) {
    c3_eb_a[c3_i43] = c3_R1[c3_i43];
  }

  for (c3_i44 = 0; c3_i44 < 9; c3_i44++) {
    c3_n_b[c3_i44] = c3_R2[c3_i44];
  }

  c3_b_eml_scalar_eg(chartInstance);
  c3_b_eml_scalar_eg(chartInstance);
  for (c3_i45 = 0; c3_i45 < 3; c3_i45++) {
    c3_i46 = 0;
    for (c3_i47 = 0; c3_i47 < 3; c3_i47++) {
      c3_hb_y[c3_i46 + c3_i45] = 0.0;
      c3_i48 = 0;
      for (c3_i49 = 0; c3_i49 < 3; c3_i49++) {
        c3_hb_y[c3_i46 + c3_i45] += c3_eb_a[c3_i48 + c3_i45] * c3_n_b[c3_i49 +
          c3_i46];
        c3_i48 += 3;
      }

      c3_i46 += 3;
    }
  }

  for (c3_i50 = 0; c3_i50 < 9; c3_i50++) {
    c3_n_b[c3_i50] = c3_R3[c3_i50];
  }

  c3_b_eml_scalar_eg(chartInstance);
  c3_b_eml_scalar_eg(chartInstance);
  for (c3_i51 = 0; c3_i51 < 9; c3_i51++) {
    c3_R_be[c3_i51] = 0.0;
  }

  for (c3_i52 = 0; c3_i52 < 9; c3_i52++) {
    c3_R_be[c3_i52] = 0.0;
  }

  for (c3_i53 = 0; c3_i53 < 9; c3_i53++) {
    c3_eb_a[c3_i53] = c3_R_be[c3_i53];
  }

  for (c3_i54 = 0; c3_i54 < 9; c3_i54++) {
    c3_R_be[c3_i54] = c3_eb_a[c3_i54];
  }

  for (c3_i55 = 0; c3_i55 < 9; c3_i55++) {
    c3_eb_a[c3_i55] = c3_R_be[c3_i55];
  }

  for (c3_i56 = 0; c3_i56 < 9; c3_i56++) {
    c3_R_be[c3_i56] = c3_eb_a[c3_i56];
  }

  for (c3_i57 = 0; c3_i57 < 3; c3_i57++) {
    c3_i58 = 0;
    for (c3_i59 = 0; c3_i59 < 3; c3_i59++) {
      c3_R_be[c3_i58 + c3_i57] = 0.0;
      c3_i60 = 0;
      for (c3_i61 = 0; c3_i61 < 3; c3_i61++) {
        c3_R_be[c3_i58 + c3_i57] += c3_hb_y[c3_i60 + c3_i57] * c3_n_b[c3_i61 +
          c3_i58];
        c3_i60 += 3;
      }

      c3_i58 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 62);
  c3_i62 = 0;
  for (c3_i63 = 0; c3_i63 < 3; c3_i63++) {
    c3_i64 = 0;
    for (c3_i65 = 0; c3_i65 < 3; c3_i65++) {
      c3_R_eb[c3_i65 + c3_i62] = c3_R_be[c3_i64 + c3_i63];
      c3_i64 += 3;
    }

    c3_i62 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 64);
  for (c3_i66 = 0; c3_i66 < 9; c3_i66++) {
    c3_eb_a[c3_i66] = c3_R_eb[c3_i66];
  }

  for (c3_i67 = 0; c3_i67 < 3; c3_i67++) {
    c3_o_b[c3_i67] = c3_v_B[c3_i67];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i68 = 0; c3_i68 < 3; c3_i68++) {
    c3_x_E_dot[c3_i68] = 0.0;
  }

  for (c3_i69 = 0; c3_i69 < 3; c3_i69++) {
    c3_x_E_dot[c3_i69] = 0.0;
  }

  for (c3_i70 = 0; c3_i70 < 3; c3_i70++) {
    c3_C[c3_i70] = c3_x_E_dot[c3_i70];
  }

  for (c3_i71 = 0; c3_i71 < 3; c3_i71++) {
    c3_x_E_dot[c3_i71] = c3_C[c3_i71];
  }

  for (c3_i72 = 0; c3_i72 < 3; c3_i72++) {
    c3_C[c3_i72] = c3_x_E_dot[c3_i72];
  }

  for (c3_i73 = 0; c3_i73 < 3; c3_i73++) {
    c3_x_E_dot[c3_i73] = c3_C[c3_i73];
  }

  for (c3_i74 = 0; c3_i74 < 3; c3_i74++) {
    c3_x_E_dot[c3_i74] = 0.0;
    c3_i75 = 0;
    for (c3_i76 = 0; c3_i76 < 3; c3_i76++) {
      c3_x_E_dot[c3_i74] += c3_eb_a[c3_i75 + c3_i74] * c3_o_b[c3_i76];
      c3_i75 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 66);
  c3_y_x = c3_phi;
  c3_ab_x = c3_y_x;
  c3_ab_x = muDoubleScalarSin(c3_ab_x);
  c3_bb_x = c3_theta;
  c3_cb_x = c3_bb_x;
  c3_cb_x = muDoubleScalarTan(c3_cb_x);
  c3_fb_a = c3_ab_x;
  c3_p_b = c3_cb_x;
  c3_ib_y = c3_fb_a * c3_p_b;
  c3_db_x = c3_phi;
  c3_eb_x = c3_db_x;
  c3_eb_x = muDoubleScalarCos(c3_eb_x);
  c3_fb_x = c3_theta;
  c3_gb_x = c3_fb_x;
  c3_gb_x = muDoubleScalarTan(c3_gb_x);
  c3_gb_a = c3_eb_x;
  c3_q_b = c3_gb_x;
  c3_jb_y = c3_gb_a * c3_q_b;
  c3_hb_x = c3_phi;
  c3_ib_x = c3_hb_x;
  c3_ib_x = muDoubleScalarCos(c3_ib_x);
  c3_jb_x = c3_phi;
  c3_kb_x = c3_jb_x;
  c3_kb_x = muDoubleScalarSin(c3_kb_x);
  c3_lb_x = c3_phi;
  c3_mb_x = c3_lb_x;
  c3_mb_x = muDoubleScalarSin(c3_mb_x);
  c3_hb_a = c3_mb_x;
  c3_r_b = c3_theta;
  c3_b_sec(chartInstance, &c3_r_b);
  c3_kb_y = c3_hb_a * c3_r_b;
  c3_nb_x = c3_phi;
  c3_ob_x = c3_nb_x;
  c3_ob_x = muDoubleScalarCos(c3_ob_x);
  c3_ib_a = c3_ob_x;
  c3_s_b = c3_theta;
  c3_b_sec(chartInstance, &c3_s_b);
  c3_lb_y = c3_ib_a * c3_s_b;
  c3_T[0] = 1.0;
  c3_T[3] = c3_ib_y;
  c3_T[6] = c3_jb_y;
  c3_T[1] = 0.0;
  c3_T[4] = c3_ib_x;
  c3_T[7] = -c3_kb_x;
  c3_T[2] = 0.0;
  c3_T[5] = c3_kb_y;
  c3_T[8] = c3_lb_y;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 67);
  for (c3_i77 = 0; c3_i77 < 9; c3_i77++) {
    c3_eb_a[c3_i77] = c3_T[c3_i77];
  }

  for (c3_i78 = 0; c3_i78 < 3; c3_i78++) {
    c3_o_b[c3_i78] = c3_w_B[c3_i78];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i79 = 0; c3_i79 < 3; c3_i79++) {
    c3_eta_dot[c3_i79] = 0.0;
  }

  for (c3_i80 = 0; c3_i80 < 3; c3_i80++) {
    c3_eta_dot[c3_i80] = 0.0;
  }

  for (c3_i81 = 0; c3_i81 < 3; c3_i81++) {
    c3_C[c3_i81] = c3_eta_dot[c3_i81];
  }

  for (c3_i82 = 0; c3_i82 < 3; c3_i82++) {
    c3_eta_dot[c3_i82] = c3_C[c3_i82];
  }

  for (c3_i83 = 0; c3_i83 < 3; c3_i83++) {
    c3_C[c3_i83] = c3_eta_dot[c3_i83];
  }

  for (c3_i84 = 0; c3_i84 < 3; c3_i84++) {
    c3_eta_dot[c3_i84] = c3_C[c3_i84];
  }

  for (c3_i85 = 0; c3_i85 < 3; c3_i85++) {
    c3_eta_dot[c3_i85] = 0.0;
    c3_i86 = 0;
    for (c3_i87 = 0; c3_i87 < 3; c3_i87++) {
      c3_eta_dot[c3_i85] += c3_eb_a[c3_i86 + c3_i85] * c3_o_b[c3_i87];
      c3_i86 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 69);
  for (c3_i88 = 0; c3_i88 < 3; c3_i88++) {
    c3_g1[c3_i88] = c3_dv12[c3_i88];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 73);
  for (c3_i89 = 0; c3_i89 < 3; c3_i89++) {
    c3_C[c3_i89] = c3_f_B[c3_i89];
  }

  for (c3_i90 = 0; c3_i90 < 3; c3_i90++) {
    c3_C[c3_i90] /= 2.0;
  }

  for (c3_i91 = 0; c3_i91 < 9; c3_i91++) {
    c3_eb_a[c3_i91] = c3_R_eb[c3_i91];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i92 = 0; c3_i92 < 3; c3_i92++) {
    c3_mb_y[c3_i92] = 0.0;
    c3_i93 = 0;
    for (c3_i94 = 0; c3_i94 < 3; c3_i94++) {
      c3_mb_y[c3_i92] += c3_eb_a[c3_i93 + c3_i92] * c3_C[c3_i94];
      c3_i93 += 3;
    }
  }

  for (c3_i95 = 0; c3_i95 < 3; c3_i95++) {
    c3_v_E_dot[c3_i95] = c3_g1[c3_i95] + c3_mb_y[c3_i95];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 74);
  for (c3_i96 = 0; c3_i96 < 9; c3_i96++) {
    c3_eb_a[c3_i96] = c3_R_be[c3_i96];
  }

  for (c3_i97 = 0; c3_i97 < 3; c3_i97++) {
    c3_o_b[c3_i97] = c3_v_E_dot[c3_i97];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i98 = 0; c3_i98 < 3; c3_i98++) {
    c3_v_B_dot[c3_i98] = 0.0;
  }

  for (c3_i99 = 0; c3_i99 < 3; c3_i99++) {
    c3_v_B_dot[c3_i99] = 0.0;
  }

  for (c3_i100 = 0; c3_i100 < 3; c3_i100++) {
    c3_C[c3_i100] = c3_v_B_dot[c3_i100];
  }

  for (c3_i101 = 0; c3_i101 < 3; c3_i101++) {
    c3_v_B_dot[c3_i101] = c3_C[c3_i101];
  }

  for (c3_i102 = 0; c3_i102 < 3; c3_i102++) {
    c3_C[c3_i102] = c3_v_B_dot[c3_i102];
  }

  for (c3_i103 = 0; c3_i103 < 3; c3_i103++) {
    c3_v_B_dot[c3_i103] = c3_C[c3_i103];
  }

  for (c3_i104 = 0; c3_i104 < 3; c3_i104++) {
    c3_v_B_dot[c3_i104] = 0.0;
    c3_i105 = 0;
    for (c3_i106 = 0; c3_i106 < 3; c3_i106++) {
      c3_v_B_dot[c3_i104] += c3_eb_a[c3_i105 + c3_i104] * c3_o_b[c3_i106];
      c3_i105 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 77);
  c3_t_b = c3_w1;
  c3_nb_y = 0.005 * c3_t_b;
  c3_jb_a = c3_nb_y;
  c3_u_b = c3_q;
  c3_ob_y = c3_jb_a * c3_u_b;
  c3_v_b = c3_w1;
  c3_pb_y = -0.005 * c3_v_b;
  c3_kb_a = c3_pb_y;
  c3_w_b = c3_p;
  c3_qb_y = c3_kb_a * c3_w_b;
  c3_i1[0] = c3_ob_y;
  c3_i1[1] = c3_qb_y;
  c3_i1[2] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 78);
  c3_x_b = c3_w2;
  c3_rb_y = 0.005 * c3_x_b;
  c3_lb_a = c3_rb_y;
  c3_y_b = c3_q;
  c3_sb_y = c3_lb_a * c3_y_b;
  c3_ab_b = c3_w2;
  c3_tb_y = -0.005 * c3_ab_b;
  c3_mb_a = c3_tb_y;
  c3_bb_b = c3_p;
  c3_ub_y = c3_mb_a * c3_bb_b;
  c3_i2[0] = c3_sb_y;
  c3_i2[1] = c3_ub_y;
  c3_i2[2] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 79);
  c3_cb_b = c3_w3;
  c3_vb_y = 0.005 * c3_cb_b;
  c3_nb_a = c3_vb_y;
  c3_db_b = c3_q;
  c3_wb_y = c3_nb_a * c3_db_b;
  c3_eb_b = c3_w3;
  c3_xb_y = -0.005 * c3_eb_b;
  c3_ob_a = c3_xb_y;
  c3_fb_b = c3_p;
  c3_yb_y = c3_ob_a * c3_fb_b;
  c3_i3[0] = c3_wb_y;
  c3_i3[1] = c3_yb_y;
  c3_i3[2] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 80);
  c3_gb_b = c3_w4;
  c3_ac_y = 0.005 * c3_gb_b;
  c3_pb_a = c3_ac_y;
  c3_hb_b = c3_q;
  c3_bc_y = c3_pb_a * c3_hb_b;
  c3_ib_b = c3_w4;
  c3_cc_y = -0.005 * c3_ib_b;
  c3_qb_a = c3_cc_y;
  c3_jb_b = c3_p;
  c3_dc_y = c3_qb_a * c3_jb_b;
  c3_b_i4[0] = c3_bc_y;
  c3_b_i4[1] = c3_dc_y;
  c3_b_i4[2] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 81);
  for (c3_i107 = 0; c3_i107 < 3; c3_i107++) {
    c3_Inew[c3_i107] = ((c3_i1[c3_i107] + c3_i2[c3_i107]) + c3_i3[c3_i107]) +
      c3_b_i4[c3_i107];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 82);
  for (c3_i108 = 0; c3_i108 < 3; c3_i108++) {
    c3_o_b[c3_i108] = c3_w_B[c3_i108];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i109 = 0; c3_i109 < 3; c3_i109++) {
    c3_temp2[c3_i109] = 0.0;
  }

  for (c3_i110 = 0; c3_i110 < 3; c3_i110++) {
    c3_temp2[c3_i110] = 0.0;
  }

  for (c3_i111 = 0; c3_i111 < 3; c3_i111++) {
    c3_C[c3_i111] = c3_temp2[c3_i111];
  }

  for (c3_i112 = 0; c3_i112 < 3; c3_i112++) {
    c3_temp2[c3_i112] = c3_C[c3_i112];
  }

  for (c3_i113 = 0; c3_i113 < 3; c3_i113++) {
    c3_C[c3_i113] = c3_temp2[c3_i113];
  }

  for (c3_i114 = 0; c3_i114 < 3; c3_i114++) {
    c3_temp2[c3_i114] = c3_C[c3_i114];
  }

  for (c3_i115 = 0; c3_i115 < 3; c3_i115++) {
    c3_temp2[c3_i115] = 0.0;
    c3_i116 = 0;
    for (c3_i117 = 0; c3_i117 < 3; c3_i117++) {
      c3_temp2[c3_i115] += c3_a[c3_i116 + c3_i115] * c3_o_b[c3_i117];
      c3_i116 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 84);
  for (c3_i118 = 0; c3_i118 < 3; c3_i118++) {
    c3_b_w_B[c3_i118] = c3_w_B[c3_i118];
  }

  for (c3_i119 = 0; c3_i119 < 3; c3_i119++) {
    c3_b_temp2[c3_i119] = c3_temp2[c3_i119];
  }

  c3_cross(chartInstance, c3_b_w_B, c3_b_temp2, c3_dv13);
  for (c3_i120 = 0; c3_i120 < 3; c3_i120++) {
    c3_temp3[c3_i120] = c3_dv13[c3_i120];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 85);
  for (c3_i121 = 0; c3_i121 < 9; c3_i121++) {
    c3_rb_a[c3_i121] = c3_a[c3_i121];
  }

  c3_inv(chartInstance, c3_rb_a, c3_eb_a);
  for (c3_i122 = 0; c3_i122 < 3; c3_i122++) {
    c3_o_b[c3_i122] = c3_M_B[c3_i122] - c3_temp3[c3_i122];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i123 = 0; c3_i123 < 3; c3_i123++) {
    c3_mb_y[c3_i123] = 0.0;
    c3_i124 = 0;
    for (c3_i125 = 0; c3_i125 < 3; c3_i125++) {
      c3_mb_y[c3_i123] += c3_eb_a[c3_i124 + c3_i123] * c3_o_b[c3_i125];
      c3_i124 += 3;
    }
  }

  for (c3_i126 = 0; c3_i126 < 3; c3_i126++) {
    c3_w_B_dot[c3_i126] = c3_mb_y[c3_i126] - c3_Inew[c3_i126];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 87);
  for (c3_i127 = 0; c3_i127 < 12; c3_i127++) {
    c3_output[c3_i127] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 88);
  for (c3_i128 = 0; c3_i128 < 3; c3_i128++) {
    c3_output[c3_i128] = c3_x_E_dot[c3_i128];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 89);
  for (c3_i129 = 0; c3_i129 < 3; c3_i129++) {
    c3_output[c3_i129 + 3] = c3_eta_dot[c3_i129];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 90);
  for (c3_i130 = 0; c3_i130 < 3; c3_i130++) {
    c3_output[c3_i130 + 6] = c3_v_B_dot[c3_i130];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 91);
  for (c3_i131 = 0; c3_i131 < 3; c3_i131++) {
    c3_output[c3_i131 + 9] = c3_w_B_dot[c3_i131];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -91);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i132 = 0; c3_i132 < 12; c3_i132++) {
    (*c3_b_output)[c3_i132] = c3_output[c3_i132];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
}

static void initSimStructsc3_lab5_template_sivom
  (SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber)
{
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i133;
  real_T c3_b_inData[12];
  int32_T c3_i134;
  real_T c3_u[12];
  const mxArray *c3_y = NULL;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i133 = 0; c3_i133 < 12; c3_i133++) {
    c3_b_inData[c3_i133] = (*(real_T (*)[12])c3_inData)[c3_i133];
  }

  for (c3_i134 = 0; c3_i134 < 12; c3_i134++) {
    c3_u[c3_i134] = c3_b_inData[c3_i134];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 12), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static void c3_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_output, const char_T *c3_identifier, real_T
  c3_y[12])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_output), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_output);
}

static void c3_b_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[12])
{
  real_T c3_dv14[12];
  int32_T c3_i135;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv14, 1, 0, 0U, 1, 0U, 1, 12);
  for (c3_i135 = 0; c3_i135 < 12; c3_i135++) {
    c3_y[c3_i135] = c3_dv14[c3_i135];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_output;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[12];
  int32_T c3_i136;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_output = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_output), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_output);
  for (c3_i136 = 0; c3_i136 < 12; c3_i136++) {
    (*(real_T (*)[12])c3_outData)[c3_i136] = c3_y[c3_i136];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static real_T c3_c_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d0;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d0, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i137;
  real_T c3_b_inData[3];
  int32_T c3_i138;
  real_T c3_u[3];
  const mxArray *c3_y = NULL;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i137 = 0; c3_i137 < 3; c3_i137++) {
    c3_b_inData[c3_i137] = (*(real_T (*)[3])c3_inData)[c3_i137];
  }

  for (c3_i138 = 0; c3_i138 < 3; c3_i138++) {
    c3_u[c3_i138] = c3_b_inData[c3_i138];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static void c3_d_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[3])
{
  real_T c3_dv15[3];
  int32_T c3_i139;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv15, 1, 0, 0U, 1, 0U, 1, 3);
  for (c3_i139 = 0; c3_i139 < 3; c3_i139++) {
    c3_y[c3_i139] = c3_dv15[c3_i139];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_w_B_dot;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[3];
  int32_T c3_i140;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_w_B_dot = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_w_B_dot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_w_B_dot);
  for (c3_i140 = 0; c3_i140 < 3; c3_i140++) {
    (*(real_T (*)[3])c3_outData)[c3_i140] = c3_y[c3_i140];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i141;
  int32_T c3_i142;
  int32_T c3_i143;
  real_T c3_b_inData[9];
  int32_T c3_i144;
  int32_T c3_i145;
  int32_T c3_i146;
  real_T c3_u[9];
  const mxArray *c3_y = NULL;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i141 = 0;
  for (c3_i142 = 0; c3_i142 < 3; c3_i142++) {
    for (c3_i143 = 0; c3_i143 < 3; c3_i143++) {
      c3_b_inData[c3_i143 + c3_i141] = (*(real_T (*)[9])c3_inData)[c3_i143 +
        c3_i141];
    }

    c3_i141 += 3;
  }

  c3_i144 = 0;
  for (c3_i145 = 0; c3_i145 < 3; c3_i145++) {
    for (c3_i146 = 0; c3_i146 < 3; c3_i146++) {
      c3_u[c3_i146 + c3_i144] = c3_b_inData[c3_i146 + c3_i144];
    }

    c3_i144 += 3;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static void c3_e_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  real_T c3_y[9])
{
  real_T c3_dv16[9];
  int32_T c3_i147;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv16, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c3_i147 = 0; c3_i147 < 9; c3_i147++) {
    c3_y[c3_i147] = c3_dv16[c3_i147];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_T;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[9];
  int32_T c3_i148;
  int32_T c3_i149;
  int32_T c3_i150;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_T = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_T), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_T);
  c3_i148 = 0;
  for (c3_i149 = 0; c3_i149 < 3; c3_i149++) {
    for (c3_i150 = 0; c3_i150 < 3; c3_i150++) {
      (*(real_T (*)[9])c3_outData)[c3_i150 + c3_i148] = c3_y[c3_i150 + c3_i148];
    }

    c3_i148 += 3;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_lab5_template_sivom_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 67, 1),
                FALSE);
  c3_info_helper(&c3_nameCaptureInfo);
  c3_b_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  const mxArray *c3_rhs16 = NULL;
  const mxArray *c3_lhs16 = NULL;
  const mxArray *c3_rhs17 = NULL;
  const mxArray *c3_lhs17 = NULL;
  const mxArray *c3_rhs18 = NULL;
  const mxArray *c3_lhs18 = NULL;
  const mxArray *c3_rhs19 = NULL;
  const mxArray *c3_lhs19 = NULL;
  const mxArray *c3_rhs20 = NULL;
  const mxArray *c3_lhs20 = NULL;
  const mxArray *c3_rhs21 = NULL;
  const mxArray *c3_lhs21 = NULL;
  const mxArray *c3_rhs22 = NULL;
  const mxArray *c3_lhs22 = NULL;
  const mxArray *c3_rhs23 = NULL;
  const mxArray *c3_lhs23 = NULL;
  const mxArray *c3_rhs24 = NULL;
  const mxArray *c3_lhs24 = NULL;
  const mxArray *c3_rhs25 = NULL;
  const mxArray *c3_lhs25 = NULL;
  const mxArray *c3_rhs26 = NULL;
  const mxArray *c3_lhs26 = NULL;
  const mxArray *c3_rhs27 = NULL;
  const mxArray *c3_lhs27 = NULL;
  const mxArray *c3_rhs28 = NULL;
  const mxArray *c3_lhs28 = NULL;
  const mxArray *c3_rhs29 = NULL;
  const mxArray *c3_lhs29 = NULL;
  const mxArray *c3_rhs30 = NULL;
  const mxArray *c3_lhs30 = NULL;
  const mxArray *c3_rhs31 = NULL;
  const mxArray *c3_lhs31 = NULL;
  const mxArray *c3_rhs32 = NULL;
  const mxArray *c3_lhs32 = NULL;
  const mxArray *c3_rhs33 = NULL;
  const mxArray *c3_lhs33 = NULL;
  const mxArray *c3_rhs34 = NULL;
  const mxArray *c3_lhs34 = NULL;
  const mxArray *c3_rhs35 = NULL;
  const mxArray *c3_lhs35 = NULL;
  const mxArray *c3_rhs36 = NULL;
  const mxArray *c3_lhs36 = NULL;
  const mxArray *c3_rhs37 = NULL;
  const mxArray *c3_lhs37 = NULL;
  const mxArray *c3_rhs38 = NULL;
  const mxArray *c3_lhs38 = NULL;
  const mxArray *c3_rhs39 = NULL;
  const mxArray *c3_lhs39 = NULL;
  const mxArray *c3_rhs40 = NULL;
  const mxArray *c3_lhs40 = NULL;
  const mxArray *c3_rhs41 = NULL;
  const mxArray *c3_lhs41 = NULL;
  const mxArray *c3_rhs42 = NULL;
  const mxArray *c3_lhs42 = NULL;
  const mxArray *c3_rhs43 = NULL;
  const mxArray *c3_lhs43 = NULL;
  const mxArray *c3_rhs44 = NULL;
  const mxArray *c3_lhs44 = NULL;
  const mxArray *c3_rhs45 = NULL;
  const mxArray *c3_lhs45 = NULL;
  const mxArray *c3_rhs46 = NULL;
  const mxArray *c3_lhs46 = NULL;
  const mxArray *c3_rhs47 = NULL;
  const mxArray *c3_lhs47 = NULL;
  const mxArray *c3_rhs48 = NULL;
  const mxArray *c3_lhs48 = NULL;
  const mxArray *c3_rhs49 = NULL;
  const mxArray *c3_lhs49 = NULL;
  const mxArray *c3_rhs50 = NULL;
  const mxArray *c3_lhs50 = NULL;
  const mxArray *c3_rhs51 = NULL;
  const mxArray *c3_lhs51 = NULL;
  const mxArray *c3_rhs52 = NULL;
  const mxArray *c3_lhs52 = NULL;
  const mxArray *c3_rhs53 = NULL;
  const mxArray *c3_lhs53 = NULL;
  const mxArray *c3_rhs54 = NULL;
  const mxArray *c3_lhs54 = NULL;
  const mxArray *c3_rhs55 = NULL;
  const mxArray *c3_lhs55 = NULL;
  const mxArray *c3_rhs56 = NULL;
  const mxArray *c3_lhs56 = NULL;
  const mxArray *c3_rhs57 = NULL;
  const mxArray *c3_lhs57 = NULL;
  const mxArray *c3_rhs58 = NULL;
  const mxArray *c3_lhs58 = NULL;
  const mxArray *c3_rhs59 = NULL;
  const mxArray *c3_lhs59 = NULL;
  const mxArray *c3_rhs60 = NULL;
  const mxArray *c3_lhs60 = NULL;
  const mxArray *c3_rhs61 = NULL;
  const mxArray *c3_lhs61 = NULL;
  const mxArray *c3_rhs62 = NULL;
  const mxArray *c3_lhs62 = NULL;
  const mxArray *c3_rhs63 = NULL;
  const mxArray *c3_lhs63 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mpower"), "name", "name", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("ismatrix"), "name", "name", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1331288658U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("power"), "name", "name", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697680U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806196U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1358169940U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("floor"), "name", "name", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697654U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806126U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806196U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mtimes"), "name", "name", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mtimes"), "name", "name", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sign"), "name", "name", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697656U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sign"), "name",
                  "name", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sign.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1356525294U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c3_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cross"), "name", "name", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806242U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c3_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mtimes"), "name", "name", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c3_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("cos"), "name", "name", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343817772U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c3_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806122U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c3_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sin"), "name", "name", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343817786U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c3_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806136U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c3_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323154378U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c3_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806196U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c3_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697670U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c3_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1299060568U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c3_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mtimes"), "name", "name", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c3_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323154378U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c3_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806196U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c3_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1360266150U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c3_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("tan"), "name", "name", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tan.m"), "resolved",
                  "resolved", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343817786U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c3_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tan.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_tan"), "name",
                  "name", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_tan.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806138U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c3_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("sec"), "name", "name", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sec.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343817784U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c3_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sec.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_sec"), "name",
                  "name", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sec.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806136U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c3_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sec.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806122U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c3_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sec.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697666U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c3_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mrdivide"), "name", "name", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1373293908U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1319717366U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c3_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("rdivide"), "name", "name", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697680U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c3_rhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c3_rhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806196U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c3_rhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697666U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c3_rhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("inv"), "name", "name", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1305305400U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c3_rhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323154378U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c3_rhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("abs"), "name", "name", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697652U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c3_rhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c3_rhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806112U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c3_rhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697666U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c3_rhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs47, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mtimes"), "name", "name", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c3_rhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs48, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806178U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c3_rhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs49, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323154378U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c3_rhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs50, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("norm"), "name", "name", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697668U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c3_rhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs51, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c3_rhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs52, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("abs"), "name", "name", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697652U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c3_rhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs53, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isnan"), "name", "name", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697658U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c3_rhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs54, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363698356U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c3_rhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs55, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806176U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c3_rhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs56, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806182U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c3_rhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs57, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mtimes"), "name", "name", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697678U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c3_rhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs58, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_warning"), "name", "name",
                  59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806202U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c3_rhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs59, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isnan"), "name", "name", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363697658U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c3_rhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs60, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eps"), "name", "name", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326711796U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c3_rhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs61, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286806182U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c3_rhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs62, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_eps"), "name", "name", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326711796U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c3_rhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs63, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
  sf_mex_destroy(&c3_rhs16);
  sf_mex_destroy(&c3_lhs16);
  sf_mex_destroy(&c3_rhs17);
  sf_mex_destroy(&c3_lhs17);
  sf_mex_destroy(&c3_rhs18);
  sf_mex_destroy(&c3_lhs18);
  sf_mex_destroy(&c3_rhs19);
  sf_mex_destroy(&c3_lhs19);
  sf_mex_destroy(&c3_rhs20);
  sf_mex_destroy(&c3_lhs20);
  sf_mex_destroy(&c3_rhs21);
  sf_mex_destroy(&c3_lhs21);
  sf_mex_destroy(&c3_rhs22);
  sf_mex_destroy(&c3_lhs22);
  sf_mex_destroy(&c3_rhs23);
  sf_mex_destroy(&c3_lhs23);
  sf_mex_destroy(&c3_rhs24);
  sf_mex_destroy(&c3_lhs24);
  sf_mex_destroy(&c3_rhs25);
  sf_mex_destroy(&c3_lhs25);
  sf_mex_destroy(&c3_rhs26);
  sf_mex_destroy(&c3_lhs26);
  sf_mex_destroy(&c3_rhs27);
  sf_mex_destroy(&c3_lhs27);
  sf_mex_destroy(&c3_rhs28);
  sf_mex_destroy(&c3_lhs28);
  sf_mex_destroy(&c3_rhs29);
  sf_mex_destroy(&c3_lhs29);
  sf_mex_destroy(&c3_rhs30);
  sf_mex_destroy(&c3_lhs30);
  sf_mex_destroy(&c3_rhs31);
  sf_mex_destroy(&c3_lhs31);
  sf_mex_destroy(&c3_rhs32);
  sf_mex_destroy(&c3_lhs32);
  sf_mex_destroy(&c3_rhs33);
  sf_mex_destroy(&c3_lhs33);
  sf_mex_destroy(&c3_rhs34);
  sf_mex_destroy(&c3_lhs34);
  sf_mex_destroy(&c3_rhs35);
  sf_mex_destroy(&c3_lhs35);
  sf_mex_destroy(&c3_rhs36);
  sf_mex_destroy(&c3_lhs36);
  sf_mex_destroy(&c3_rhs37);
  sf_mex_destroy(&c3_lhs37);
  sf_mex_destroy(&c3_rhs38);
  sf_mex_destroy(&c3_lhs38);
  sf_mex_destroy(&c3_rhs39);
  sf_mex_destroy(&c3_lhs39);
  sf_mex_destroy(&c3_rhs40);
  sf_mex_destroy(&c3_lhs40);
  sf_mex_destroy(&c3_rhs41);
  sf_mex_destroy(&c3_lhs41);
  sf_mex_destroy(&c3_rhs42);
  sf_mex_destroy(&c3_lhs42);
  sf_mex_destroy(&c3_rhs43);
  sf_mex_destroy(&c3_lhs43);
  sf_mex_destroy(&c3_rhs44);
  sf_mex_destroy(&c3_lhs44);
  sf_mex_destroy(&c3_rhs45);
  sf_mex_destroy(&c3_lhs45);
  sf_mex_destroy(&c3_rhs46);
  sf_mex_destroy(&c3_lhs46);
  sf_mex_destroy(&c3_rhs47);
  sf_mex_destroy(&c3_lhs47);
  sf_mex_destroy(&c3_rhs48);
  sf_mex_destroy(&c3_lhs48);
  sf_mex_destroy(&c3_rhs49);
  sf_mex_destroy(&c3_lhs49);
  sf_mex_destroy(&c3_rhs50);
  sf_mex_destroy(&c3_lhs50);
  sf_mex_destroy(&c3_rhs51);
  sf_mex_destroy(&c3_lhs51);
  sf_mex_destroy(&c3_rhs52);
  sf_mex_destroy(&c3_lhs52);
  sf_mex_destroy(&c3_rhs53);
  sf_mex_destroy(&c3_lhs53);
  sf_mex_destroy(&c3_rhs54);
  sf_mex_destroy(&c3_lhs54);
  sf_mex_destroy(&c3_rhs55);
  sf_mex_destroy(&c3_lhs55);
  sf_mex_destroy(&c3_rhs56);
  sf_mex_destroy(&c3_lhs56);
  sf_mex_destroy(&c3_rhs57);
  sf_mex_destroy(&c3_lhs57);
  sf_mex_destroy(&c3_rhs58);
  sf_mex_destroy(&c3_lhs58);
  sf_mex_destroy(&c3_rhs59);
  sf_mex_destroy(&c3_lhs59);
  sf_mex_destroy(&c3_rhs60);
  sf_mex_destroy(&c3_lhs60);
  sf_mex_destroy(&c3_rhs61);
  sf_mex_destroy(&c3_lhs61);
  sf_mex_destroy(&c3_rhs62);
  sf_mex_destroy(&c3_lhs62);
  sf_mex_destroy(&c3_rhs63);
  sf_mex_destroy(&c3_lhs63);
}

static const mxArray *c3_emlrt_marshallOut(char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), FALSE);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c3_y;
}

static void c3_b_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs64 = NULL;
  const mxArray *c3_lhs64 = NULL;
  const mxArray *c3_rhs65 = NULL;
  const mxArray *c3_lhs65 = NULL;
  const mxArray *c3_rhs66 = NULL;
  const mxArray *c3_lhs66 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326711796U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c3_rhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs64, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1360266150U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c3_rhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs65, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "name", "name", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1319717368U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c3_rhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c3_lhs66, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs66), "lhs", "lhs",
                  66);
  sf_mex_destroy(&c3_rhs64);
  sf_mex_destroy(&c3_lhs64);
  sf_mex_destroy(&c3_rhs65);
  sf_mex_destroy(&c3_lhs65);
  sf_mex_destroy(&c3_rhs66);
  sf_mex_destroy(&c3_lhs66);
}

static real_T c3_mpower(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
  real_T c3_a)
{
  real_T c3_b_a;
  real_T c3_c_a;
  real_T c3_ak;
  real_T c3_d_a;
  real_T c3_e_a;
  real_T c3_b;
  c3_b_a = c3_a;
  c3_c_a = c3_b_a;
  c3_eml_scalar_eg(chartInstance);
  c3_ak = c3_c_a;
  c3_d_a = c3_ak;
  c3_eml_scalar_eg(chartInstance);
  c3_e_a = c3_d_a;
  c3_b = c3_d_a;
  return c3_e_a * c3_b;
}

static void c3_eml_scalar_eg(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance)
{
}

static real_T c3_sign(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  c3_b_sign(chartInstance, &c3_b_x);
  return c3_b_x;
}

static void c3_cross(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T c3_a[3], real_T c3_b[3], real_T c3_c[3])
{
  real_T c3_b_a;
  real_T c3_b_b;
  real_T c3_y;
  real_T c3_c_a;
  real_T c3_c_b;
  real_T c3_b_y;
  real_T c3_c1;
  real_T c3_d_a;
  real_T c3_d_b;
  real_T c3_c_y;
  real_T c3_e_a;
  real_T c3_e_b;
  real_T c3_d_y;
  real_T c3_c2;
  real_T c3_f_a;
  real_T c3_f_b;
  real_T c3_e_y;
  real_T c3_g_a;
  real_T c3_g_b;
  real_T c3_f_y;
  real_T c3_c3;
  c3_b_a = c3_a[1];
  c3_b_b = c3_b[2];
  c3_y = c3_b_a * c3_b_b;
  c3_c_a = c3_a[2];
  c3_c_b = c3_b[1];
  c3_b_y = c3_c_a * c3_c_b;
  c3_c1 = c3_y - c3_b_y;
  c3_d_a = c3_a[2];
  c3_d_b = c3_b[0];
  c3_c_y = c3_d_a * c3_d_b;
  c3_e_a = c3_a[0];
  c3_e_b = c3_b[2];
  c3_d_y = c3_e_a * c3_e_b;
  c3_c2 = c3_c_y - c3_d_y;
  c3_f_a = c3_a[0];
  c3_f_b = c3_b[1];
  c3_e_y = c3_f_a * c3_f_b;
  c3_g_a = c3_a[1];
  c3_g_b = c3_b[0];
  c3_f_y = c3_g_a * c3_g_b;
  c3_c3 = c3_e_y - c3_f_y;
  c3_c[0] = c3_c1;
  c3_c[1] = c3_c2;
  c3_c[2] = c3_c3;
}

static void c3_b_eml_scalar_eg(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance)
{
}

static void c3_c_eml_scalar_eg(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance)
{
}

static real_T c3_sec(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  c3_b_sec(chartInstance, &c3_b_x);
  return c3_b_x;
}

static real_T c3_eml_div(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
  real_T c3_x, real_T c3_y)
{
  return c3_x / c3_y;
}

static void c3_inv(SFc3_lab5_template_sivomInstanceStruct *chartInstance, real_T
                   c3_x[9], real_T c3_y[9])
{
  int32_T c3_i151;
  real_T c3_b_x[9];
  int32_T c3_i152;
  real_T c3_c_x[9];
  real_T c3_n1x;
  int32_T c3_i153;
  real_T c3_b_y[9];
  real_T c3_n1xinv;
  real_T c3_a;
  real_T c3_b;
  real_T c3_c_y;
  real_T c3_rc;
  real_T c3_d_x;
  boolean_T c3_b_b;
  real_T c3_e_x;
  int32_T c3_i154;
  static char_T c3_cv0[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c3_u[8];
  const mxArray *c3_d_y = NULL;
  real_T c3_b_u;
  const mxArray *c3_e_y = NULL;
  real_T c3_c_u;
  const mxArray *c3_f_y = NULL;
  real_T c3_d_u;
  const mxArray *c3_g_y = NULL;
  char_T c3_str[14];
  int32_T c3_i155;
  char_T c3_b_str[14];
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  boolean_T guard3 = FALSE;
  for (c3_i151 = 0; c3_i151 < 9; c3_i151++) {
    c3_b_x[c3_i151] = c3_x[c3_i151];
  }

  c3_inv3x3(chartInstance, c3_b_x, c3_y);
  for (c3_i152 = 0; c3_i152 < 9; c3_i152++) {
    c3_c_x[c3_i152] = c3_x[c3_i152];
  }

  c3_n1x = c3_norm(chartInstance, c3_c_x);
  for (c3_i153 = 0; c3_i153 < 9; c3_i153++) {
    c3_b_y[c3_i153] = c3_y[c3_i153];
  }

  c3_n1xinv = c3_norm(chartInstance, c3_b_y);
  c3_a = c3_n1x;
  c3_b = c3_n1xinv;
  c3_c_y = c3_a * c3_b;
  c3_rc = 1.0 / c3_c_y;
  guard1 = FALSE;
  guard2 = FALSE;
  if (c3_n1x == 0.0) {
    guard2 = TRUE;
  } else if (c3_n1xinv == 0.0) {
    guard2 = TRUE;
  } else if (c3_rc == 0.0) {
    guard1 = TRUE;
  } else {
    c3_d_x = c3_rc;
    c3_b_b = muDoubleScalarIsNaN(c3_d_x);
    guard3 = FALSE;
    if (c3_b_b) {
      guard3 = TRUE;
    } else {
      if (c3_rc < 2.2204460492503131E-16) {
        guard3 = TRUE;
      }
    }

    if (guard3 == TRUE) {
      c3_e_x = c3_rc;
      for (c3_i154 = 0; c3_i154 < 8; c3_i154++) {
        c3_u[c3_i154] = c3_cv0[c3_i154];
      }

      c3_d_y = NULL;
      sf_mex_assign(&c3_d_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    FALSE);
      c3_b_u = 14.0;
      c3_e_y = NULL;
      sf_mex_assign(&c3_e_y, sf_mex_create("y", &c3_b_u, 0, 0U, 0U, 0U, 0),
                    FALSE);
      c3_c_u = 6.0;
      c3_f_y = NULL;
      sf_mex_assign(&c3_f_y, sf_mex_create("y", &c3_c_u, 0, 0U, 0U, 0U, 0),
                    FALSE);
      c3_d_u = c3_e_x;
      c3_g_y = NULL;
      sf_mex_assign(&c3_g_y, sf_mex_create("y", &c3_d_u, 0, 0U, 0U, 0U, 0),
                    FALSE);
      c3_f_emlrt_marshallIn(chartInstance, sf_mex_call_debug("sprintf", 1U, 2U,
        14, sf_mex_call_debug("sprintf", 1U, 3U, 14, c3_d_y, 14, c3_e_y, 14,
        c3_f_y), 14, c3_g_y), "sprintf", c3_str);
      for (c3_i155 = 0; c3_i155 < 14; c3_i155++) {
        c3_b_str[c3_i155] = c3_str[c3_i155];
      }

      c3_b_eml_warning(chartInstance, c3_b_str);
    }
  }

  if (guard2 == TRUE) {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    c3_eml_warning(chartInstance);
  }
}

static void c3_inv3x3(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T c3_x[9], real_T c3_y[9])
{
  int32_T c3_p1;
  int32_T c3_p2;
  int32_T c3_p3;
  real_T c3_absx11;
  real_T c3_absx21;
  real_T c3_absx31;
  real_T c3_t1;
  real_T c3_a;
  real_T c3_b;
  real_T c3_b_y;
  real_T c3_b_a;
  real_T c3_b_b;
  real_T c3_c_y;
  real_T c3_c_a;
  real_T c3_c_b;
  real_T c3_d_y;
  real_T c3_d_a;
  real_T c3_d_b;
  real_T c3_e_y;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_f_y;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_g_y;
  int32_T c3_itmp;
  real_T c3_f_x;
  real_T c3_h_y;
  real_T c3_z;
  real_T c3_e_a;
  real_T c3_e_b;
  real_T c3_i_y;
  real_T c3_f_a;
  real_T c3_f_b;
  real_T c3_j_y;
  real_T c3_g_x;
  real_T c3_k_y;
  real_T c3_t3;
  real_T c3_g_a;
  real_T c3_g_b;
  real_T c3_l_y;
  real_T c3_h_x;
  real_T c3_m_y;
  real_T c3_t2;
  int32_T c3_h_a;
  int32_T c3_c;
  real_T c3_i_a;
  real_T c3_h_b;
  real_T c3_n_y;
  real_T c3_j_a;
  real_T c3_i_b;
  real_T c3_o_y;
  real_T c3_i_x;
  real_T c3_p_y;
  real_T c3_b_z;
  int32_T c3_k_a;
  int32_T c3_b_c;
  int32_T c3_l_a;
  int32_T c3_c_c;
  real_T c3_j_x;
  real_T c3_q_y;
  real_T c3_m_a;
  real_T c3_j_b;
  real_T c3_r_y;
  real_T c3_k_x;
  real_T c3_s_y;
  int32_T c3_n_a;
  int32_T c3_d_c;
  real_T c3_o_a;
  real_T c3_k_b;
  real_T c3_t_y;
  real_T c3_p_a;
  real_T c3_l_b;
  real_T c3_u_y;
  real_T c3_l_x;
  real_T c3_v_y;
  real_T c3_c_z;
  int32_T c3_q_a;
  int32_T c3_e_c;
  int32_T c3_r_a;
  int32_T c3_f_c;
  real_T c3_w_y;
  real_T c3_s_a;
  real_T c3_m_b;
  real_T c3_x_y;
  real_T c3_m_x;
  real_T c3_y_y;
  int32_T c3_t_a;
  int32_T c3_g_c;
  real_T c3_u_a;
  real_T c3_n_b;
  real_T c3_ab_y;
  real_T c3_v_a;
  real_T c3_o_b;
  real_T c3_bb_y;
  real_T c3_n_x;
  real_T c3_cb_y;
  real_T c3_d_z;
  int32_T c3_w_a;
  int32_T c3_h_c;
  int32_T c3_x_a;
  int32_T c3_i_c;
  boolean_T guard1 = FALSE;
  c3_p1 = 0;
  c3_p2 = 3;
  c3_p3 = 6;
  c3_absx11 = c3_abs(chartInstance, c3_x[0]);
  c3_absx21 = c3_abs(chartInstance, c3_x[1]);
  c3_absx31 = c3_abs(chartInstance, c3_x[2]);
  guard1 = FALSE;
  if (c3_absx21 > c3_absx11) {
    if (c3_absx21 > c3_absx31) {
      c3_p1 = 3;
      c3_p2 = 0;
      c3_t1 = c3_x[0];
      c3_x[0] = c3_x[1];
      c3_x[1] = c3_t1;
      c3_t1 = c3_x[3];
      c3_x[3] = c3_x[4];
      c3_x[4] = c3_t1;
      c3_t1 = c3_x[6];
      c3_x[6] = c3_x[7];
      c3_x[7] = c3_t1;
    } else {
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    if (c3_absx31 > c3_absx11) {
      c3_p1 = 6;
      c3_p3 = 0;
      c3_t1 = c3_x[0];
      c3_x[0] = c3_x[2];
      c3_x[2] = c3_t1;
      c3_t1 = c3_x[3];
      c3_x[3] = c3_x[5];
      c3_x[5] = c3_t1;
      c3_t1 = c3_x[6];
      c3_x[6] = c3_x[8];
      c3_x[8] = c3_t1;
    }
  }

  c3_x[1] = c3_eml_div(chartInstance, c3_x[1], c3_x[0]);
  c3_x[2] = c3_eml_div(chartInstance, c3_x[2], c3_x[0]);
  c3_a = c3_x[1];
  c3_b = c3_x[3];
  c3_b_y = c3_a * c3_b;
  c3_x[4] -= c3_b_y;
  c3_b_a = c3_x[2];
  c3_b_b = c3_x[3];
  c3_c_y = c3_b_a * c3_b_b;
  c3_x[5] -= c3_c_y;
  c3_c_a = c3_x[1];
  c3_c_b = c3_x[6];
  c3_d_y = c3_c_a * c3_c_b;
  c3_x[7] -= c3_d_y;
  c3_d_a = c3_x[2];
  c3_d_b = c3_x[6];
  c3_e_y = c3_d_a * c3_d_b;
  c3_x[8] -= c3_e_y;
  c3_b_x = c3_x[5];
  c3_c_x = c3_b_x;
  c3_f_y = muDoubleScalarAbs(c3_c_x);
  c3_d_x = c3_x[4];
  c3_e_x = c3_d_x;
  c3_g_y = muDoubleScalarAbs(c3_e_x);
  if (c3_f_y > c3_g_y) {
    c3_itmp = c3_p2;
    c3_p2 = c3_p3;
    c3_p3 = c3_itmp;
    c3_t1 = c3_x[1];
    c3_x[1] = c3_x[2];
    c3_x[2] = c3_t1;
    c3_t1 = c3_x[4];
    c3_x[4] = c3_x[5];
    c3_x[5] = c3_t1;
    c3_t1 = c3_x[7];
    c3_x[7] = c3_x[8];
    c3_x[8] = c3_t1;
  }

  c3_f_x = c3_x[5];
  c3_h_y = c3_x[4];
  c3_z = c3_f_x / c3_h_y;
  c3_x[5] = c3_z;
  c3_e_a = c3_x[5];
  c3_e_b = c3_x[7];
  c3_i_y = c3_e_a * c3_e_b;
  c3_x[8] -= c3_i_y;
  c3_f_a = c3_x[5];
  c3_f_b = c3_x[1];
  c3_j_y = c3_f_a * c3_f_b;
  c3_g_x = c3_j_y - c3_x[2];
  c3_k_y = c3_x[8];
  c3_t3 = c3_g_x / c3_k_y;
  c3_g_a = c3_x[7];
  c3_g_b = c3_t3;
  c3_l_y = c3_g_a * c3_g_b;
  c3_h_x = -(c3_x[1] + c3_l_y);
  c3_m_y = c3_x[4];
  c3_t2 = c3_h_x / c3_m_y;
  c3_h_a = c3_p1 + 1;
  c3_c = c3_h_a;
  c3_i_a = c3_x[3];
  c3_h_b = c3_t2;
  c3_n_y = c3_i_a * c3_h_b;
  c3_j_a = c3_x[6];
  c3_i_b = c3_t3;
  c3_o_y = c3_j_a * c3_i_b;
  c3_i_x = (1.0 - c3_n_y) - c3_o_y;
  c3_p_y = c3_x[0];
  c3_b_z = c3_i_x / c3_p_y;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_c), 1, 9, 1, 0) - 1] = c3_b_z;
  c3_k_a = c3_p1 + 2;
  c3_b_c = c3_k_a;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_b_c), 1, 9, 1, 0) - 1] = c3_t2;
  c3_l_a = c3_p1 + 3;
  c3_c_c = c3_l_a;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_c_c), 1, 9, 1, 0) - 1] = c3_t3;
  c3_j_x = -c3_x[5];
  c3_q_y = c3_x[8];
  c3_t3 = c3_j_x / c3_q_y;
  c3_m_a = c3_x[7];
  c3_j_b = c3_t3;
  c3_r_y = c3_m_a * c3_j_b;
  c3_k_x = 1.0 - c3_r_y;
  c3_s_y = c3_x[4];
  c3_t2 = c3_k_x / c3_s_y;
  c3_n_a = c3_p2 + 1;
  c3_d_c = c3_n_a;
  c3_o_a = c3_x[3];
  c3_k_b = c3_t2;
  c3_t_y = c3_o_a * c3_k_b;
  c3_p_a = c3_x[6];
  c3_l_b = c3_t3;
  c3_u_y = c3_p_a * c3_l_b;
  c3_l_x = -(c3_t_y + c3_u_y);
  c3_v_y = c3_x[0];
  c3_c_z = c3_l_x / c3_v_y;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_d_c), 1, 9, 1, 0) - 1] = c3_c_z;
  c3_q_a = c3_p2 + 2;
  c3_e_c = c3_q_a;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_e_c), 1, 9, 1, 0) - 1] = c3_t2;
  c3_r_a = c3_p2 + 3;
  c3_f_c = c3_r_a;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_f_c), 1, 9, 1, 0) - 1] = c3_t3;
  c3_w_y = c3_x[8];
  c3_t3 = 1.0 / c3_w_y;
  c3_s_a = -c3_x[7];
  c3_m_b = c3_t3;
  c3_x_y = c3_s_a * c3_m_b;
  c3_m_x = c3_x_y;
  c3_y_y = c3_x[4];
  c3_t2 = c3_m_x / c3_y_y;
  c3_t_a = c3_p3 + 1;
  c3_g_c = c3_t_a;
  c3_u_a = c3_x[3];
  c3_n_b = c3_t2;
  c3_ab_y = c3_u_a * c3_n_b;
  c3_v_a = c3_x[6];
  c3_o_b = c3_t3;
  c3_bb_y = c3_v_a * c3_o_b;
  c3_n_x = -(c3_ab_y + c3_bb_y);
  c3_cb_y = c3_x[0];
  c3_d_z = c3_n_x / c3_cb_y;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_g_c), 1, 9, 1, 0) - 1] = c3_d_z;
  c3_w_a = c3_p3 + 2;
  c3_h_c = c3_w_a;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_h_c), 1, 9, 1, 0) - 1] = c3_t2;
  c3_x_a = c3_p3 + 3;
  c3_i_c = c3_x_a;
  c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c3_i_c), 1, 9, 1, 0) - 1] = c3_t3;
}

static real_T c3_abs(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  return muDoubleScalarAbs(c3_b_x);
}

static real_T c3_norm(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T c3_x[9])
{
  real_T c3_y;
  int32_T c3_j;
  real_T c3_b_j;
  real_T c3_s;
  int32_T c3_i;
  real_T c3_b_i;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_b_y;
  real_T c3_d_x;
  boolean_T c3_b;
  boolean_T exitg1;
  c3_y = 0.0;
  c3_j = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (c3_j < 3)) {
    c3_b_j = 1.0 + (real_T)c3_j;
    c3_s = 0.0;
    for (c3_i = 0; c3_i < 3; c3_i++) {
      c3_b_i = 1.0 + (real_T)c3_i;
      c3_b_x = c3_x[((int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c3_b_i), 1, 3, 1, 0) + 3 * ((int32_T)(real_T)
        _SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", c3_b_j),
        1, 3, 2, 0) - 1)) - 1];
      c3_c_x = c3_b_x;
      c3_b_y = muDoubleScalarAbs(c3_c_x);
      c3_s += c3_b_y;
    }

    c3_d_x = c3_s;
    c3_b = muDoubleScalarIsNaN(c3_d_x);
    if (c3_b) {
      c3_y = rtNaN;
      exitg1 = TRUE;
    } else {
      if (c3_s > c3_y) {
        c3_y = c3_s;
      }

      c3_j++;
    }
  }

  return c3_y;
}

static void c3_eml_warning(SFc3_lab5_template_sivomInstanceStruct *chartInstance)
{
  int32_T c3_i156;
  static char_T c3_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c3_u[27];
  const mxArray *c3_y = NULL;
  for (c3_i156 = 0; c3_i156 < 27; c3_i156++) {
    c3_u[c3_i156] = c3_varargin_1[c3_i156];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c3_y));
}

static void c3_b_eml_warning(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, char_T c3_varargin_2[14])
{
  int32_T c3_i157;
  static char_T c3_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c3_u[33];
  const mxArray *c3_y = NULL;
  int32_T c3_i158;
  char_T c3_b_u[14];
  const mxArray *c3_b_y = NULL;
  for (c3_i157 = 0; c3_i157 < 33; c3_i157++) {
    c3_u[c3_i157] = c3_varargin_1[c3_i157];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 33), FALSE);
  for (c3_i158 = 0; c3_i158 < 14; c3_i158++) {
    c3_b_u[c3_i158] = c3_varargin_2[c3_i158];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U,
    14, c3_y, 14, c3_b_y));
}

static void c3_f_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_sprintf, const char_T *c3_identifier, char_T
  c3_y[14])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_sprintf), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_sprintf);
}

static void c3_g_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId,
  char_T c3_y[14])
{
  char_T c3_cv1[14];
  int32_T c3_i159;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_cv1, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c3_i159 = 0; c3_i159 < 14; c3_i159++) {
    c3_y[c3_i159] = c3_cv1[c3_i159];
  }

  sf_mex_destroy(&c3_u);
}

static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static int32_T c3_h_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i160;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i160, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i160;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_i_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_b_is_active_c3_lab5_template_sivom, const
  char_T *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_j_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_b_is_active_c3_lab5_template_sivom), &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_lab5_template_sivom);
  return c3_y;
}

static uint8_T c3_j_emlrt_marshallIn(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sign(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                      real_T *c3_x)
{
  *c3_x = muDoubleScalarSign(*c3_x);
}

static void c3_b_sec(SFc3_lab5_template_sivomInstanceStruct *chartInstance,
                     real_T *c3_x)
{
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_y;
  c3_b_x = *c3_x;
  c3_c_x = c3_b_x;
  c3_c_x = muDoubleScalarCos(c3_c_x);
  c3_y = c3_c_x;
  *c3_x = 1.0 / c3_y;
}

static void init_dsm_address_info(SFc3_lab5_template_sivomInstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c3_lab5_template_sivom_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(475828561U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(4148236346U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1846840142U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2311060526U);
}

mxArray *sf_c3_lab5_template_sivom_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("DMfTTpQ8VzlKosSIHh7rKG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c3_lab5_template_sivom_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_lab5_template_sivom_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c3_lab5_template_sivom(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"output\",},{M[8],M[0],T\"is_active_c3_lab5_template_sivom\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_lab5_template_sivom_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_lab5_template_sivomInstanceStruct *chartInstance;
    chartInstance = (SFc3_lab5_template_sivomInstanceStruct *) ((ChartInfoStruct
      *)(ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _lab5_template_sivomMachineNumber_,
           3,
           1,
           1,
           6,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_lab5_template_sivomMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_lab5_template_sivomMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _lab5_template_sivomMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"w1");
          _SFD_SET_DATA_PROPS(1,2,0,1,"output");
          _SFD_SET_DATA_PROPS(2,1,1,0,"w2");
          _SFD_SET_DATA_PROPS(3,1,1,0,"w3");
          _SFD_SET_DATA_PROPS(4,1,1,0,"w4");
          _SFD_SET_DATA_PROPS(5,1,1,0,"states");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,1810);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c3_w1;
          real_T *c3_w2;
          real_T *c3_w3;
          real_T *c3_w4;
          real_T (*c3_output)[12];
          real_T (*c3_states)[12];
          c3_states = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 4);
          c3_w4 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c3_w3 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c3_w2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c3_output = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
          c3_w1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c3_w1);
          _SFD_SET_DATA_VALUE_PTR(1U, *c3_output);
          _SFD_SET_DATA_VALUE_PTR(2U, c3_w2);
          _SFD_SET_DATA_VALUE_PTR(3U, c3_w3);
          _SFD_SET_DATA_VALUE_PTR(4U, c3_w4);
          _SFD_SET_DATA_VALUE_PTR(5U, *c3_states);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _lab5_template_sivomMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "XzHnY7mYP8yigbH9Lzf8tF";
}

static void sf_opaque_initialize_c3_lab5_template_sivom(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_lab5_template_sivomInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c3_lab5_template_sivom
    ((SFc3_lab5_template_sivomInstanceStruct*) chartInstanceVar);
  initialize_c3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c3_lab5_template_sivom(void *chartInstanceVar)
{
  enable_c3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c3_lab5_template_sivom(void *chartInstanceVar)
{
  disable_c3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c3_lab5_template_sivom(void *chartInstanceVar)
{
  sf_c3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c3_lab5_template_sivom(SimStruct*
  S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c3_lab5_template_sivom
    ((SFc3_lab5_template_sivomInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_lab5_template_sivom();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c3_lab5_template_sivom(SimStruct* S, const
  mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_lab5_template_sivom();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c3_lab5_template_sivom(SimStruct*
  S)
{
  return sf_internal_get_sim_state_c3_lab5_template_sivom(S);
}

static void sf_opaque_set_sim_state_c3_lab5_template_sivom(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c3_lab5_template_sivom(S, st);
}

static void sf_opaque_terminate_c3_lab5_template_sivom(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_lab5_template_sivomInstanceStruct*) chartInstanceVar
      )->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_lab5_template_sivom_optimization_info();
    }

    finalize_c3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
      chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_lab5_template_sivom((SFc3_lab5_template_sivomInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_lab5_template_sivom(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c3_lab5_template_sivom
      ((SFc3_lab5_template_sivomInstanceStruct*)(((ChartInfoStruct *)
         ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_lab5_template_sivom(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_lab5_template_sivom_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,3,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1894577519U));
  ssSetChecksum1(S,(3596603238U));
  ssSetChecksum2(S,(40856549U));
  ssSetChecksum3(S,(3392902966U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_lab5_template_sivom(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_lab5_template_sivom(SimStruct *S)
{
  SFc3_lab5_template_sivomInstanceStruct *chartInstance;
  chartInstance = (SFc3_lab5_template_sivomInstanceStruct *)utMalloc(sizeof
    (SFc3_lab5_template_sivomInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_lab5_template_sivomInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c3_lab5_template_sivom;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c3_lab5_template_sivom;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c3_lab5_template_sivom;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c3_lab5_template_sivom;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c3_lab5_template_sivom;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c3_lab5_template_sivom;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c3_lab5_template_sivom;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c3_lab5_template_sivom;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_lab5_template_sivom;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_lab5_template_sivom;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c3_lab5_template_sivom;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c3_lab5_template_sivom_method_dispatcher(SimStruct *S, int_T method, void
  *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_lab5_template_sivom(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_lab5_template_sivom(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_lab5_template_sivom(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_lab5_template_sivom_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
