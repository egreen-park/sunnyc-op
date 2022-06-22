#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8535675879233886937) {
   out_8535675879233886937[0] = delta_x[0] + nom_x[0];
   out_8535675879233886937[1] = delta_x[1] + nom_x[1];
   out_8535675879233886937[2] = delta_x[2] + nom_x[2];
   out_8535675879233886937[3] = delta_x[3] + nom_x[3];
   out_8535675879233886937[4] = delta_x[4] + nom_x[4];
   out_8535675879233886937[5] = delta_x[5] + nom_x[5];
   out_8535675879233886937[6] = delta_x[6] + nom_x[6];
   out_8535675879233886937[7] = delta_x[7] + nom_x[7];
   out_8535675879233886937[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6876297271512678565) {
   out_6876297271512678565[0] = -nom_x[0] + true_x[0];
   out_6876297271512678565[1] = -nom_x[1] + true_x[1];
   out_6876297271512678565[2] = -nom_x[2] + true_x[2];
   out_6876297271512678565[3] = -nom_x[3] + true_x[3];
   out_6876297271512678565[4] = -nom_x[4] + true_x[4];
   out_6876297271512678565[5] = -nom_x[5] + true_x[5];
   out_6876297271512678565[6] = -nom_x[6] + true_x[6];
   out_6876297271512678565[7] = -nom_x[7] + true_x[7];
   out_6876297271512678565[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3429589102358787544) {
   out_3429589102358787544[0] = 1.0;
   out_3429589102358787544[1] = 0;
   out_3429589102358787544[2] = 0;
   out_3429589102358787544[3] = 0;
   out_3429589102358787544[4] = 0;
   out_3429589102358787544[5] = 0;
   out_3429589102358787544[6] = 0;
   out_3429589102358787544[7] = 0;
   out_3429589102358787544[8] = 0;
   out_3429589102358787544[9] = 0;
   out_3429589102358787544[10] = 1.0;
   out_3429589102358787544[11] = 0;
   out_3429589102358787544[12] = 0;
   out_3429589102358787544[13] = 0;
   out_3429589102358787544[14] = 0;
   out_3429589102358787544[15] = 0;
   out_3429589102358787544[16] = 0;
   out_3429589102358787544[17] = 0;
   out_3429589102358787544[18] = 0;
   out_3429589102358787544[19] = 0;
   out_3429589102358787544[20] = 1.0;
   out_3429589102358787544[21] = 0;
   out_3429589102358787544[22] = 0;
   out_3429589102358787544[23] = 0;
   out_3429589102358787544[24] = 0;
   out_3429589102358787544[25] = 0;
   out_3429589102358787544[26] = 0;
   out_3429589102358787544[27] = 0;
   out_3429589102358787544[28] = 0;
   out_3429589102358787544[29] = 0;
   out_3429589102358787544[30] = 1.0;
   out_3429589102358787544[31] = 0;
   out_3429589102358787544[32] = 0;
   out_3429589102358787544[33] = 0;
   out_3429589102358787544[34] = 0;
   out_3429589102358787544[35] = 0;
   out_3429589102358787544[36] = 0;
   out_3429589102358787544[37] = 0;
   out_3429589102358787544[38] = 0;
   out_3429589102358787544[39] = 0;
   out_3429589102358787544[40] = 1.0;
   out_3429589102358787544[41] = 0;
   out_3429589102358787544[42] = 0;
   out_3429589102358787544[43] = 0;
   out_3429589102358787544[44] = 0;
   out_3429589102358787544[45] = 0;
   out_3429589102358787544[46] = 0;
   out_3429589102358787544[47] = 0;
   out_3429589102358787544[48] = 0;
   out_3429589102358787544[49] = 0;
   out_3429589102358787544[50] = 1.0;
   out_3429589102358787544[51] = 0;
   out_3429589102358787544[52] = 0;
   out_3429589102358787544[53] = 0;
   out_3429589102358787544[54] = 0;
   out_3429589102358787544[55] = 0;
   out_3429589102358787544[56] = 0;
   out_3429589102358787544[57] = 0;
   out_3429589102358787544[58] = 0;
   out_3429589102358787544[59] = 0;
   out_3429589102358787544[60] = 1.0;
   out_3429589102358787544[61] = 0;
   out_3429589102358787544[62] = 0;
   out_3429589102358787544[63] = 0;
   out_3429589102358787544[64] = 0;
   out_3429589102358787544[65] = 0;
   out_3429589102358787544[66] = 0;
   out_3429589102358787544[67] = 0;
   out_3429589102358787544[68] = 0;
   out_3429589102358787544[69] = 0;
   out_3429589102358787544[70] = 1.0;
   out_3429589102358787544[71] = 0;
   out_3429589102358787544[72] = 0;
   out_3429589102358787544[73] = 0;
   out_3429589102358787544[74] = 0;
   out_3429589102358787544[75] = 0;
   out_3429589102358787544[76] = 0;
   out_3429589102358787544[77] = 0;
   out_3429589102358787544[78] = 0;
   out_3429589102358787544[79] = 0;
   out_3429589102358787544[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_557623254591872103) {
   out_557623254591872103[0] = state[0];
   out_557623254591872103[1] = state[1];
   out_557623254591872103[2] = state[2];
   out_557623254591872103[3] = state[3];
   out_557623254591872103[4] = state[4];
   out_557623254591872103[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_557623254591872103[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_557623254591872103[7] = state[7];
   out_557623254591872103[8] = state[8];
}
void F_fun(double *state, double dt, double *out_5593547350807029480) {
   out_5593547350807029480[0] = 1;
   out_5593547350807029480[1] = 0;
   out_5593547350807029480[2] = 0;
   out_5593547350807029480[3] = 0;
   out_5593547350807029480[4] = 0;
   out_5593547350807029480[5] = 0;
   out_5593547350807029480[6] = 0;
   out_5593547350807029480[7] = 0;
   out_5593547350807029480[8] = 0;
   out_5593547350807029480[9] = 0;
   out_5593547350807029480[10] = 1;
   out_5593547350807029480[11] = 0;
   out_5593547350807029480[12] = 0;
   out_5593547350807029480[13] = 0;
   out_5593547350807029480[14] = 0;
   out_5593547350807029480[15] = 0;
   out_5593547350807029480[16] = 0;
   out_5593547350807029480[17] = 0;
   out_5593547350807029480[18] = 0;
   out_5593547350807029480[19] = 0;
   out_5593547350807029480[20] = 1;
   out_5593547350807029480[21] = 0;
   out_5593547350807029480[22] = 0;
   out_5593547350807029480[23] = 0;
   out_5593547350807029480[24] = 0;
   out_5593547350807029480[25] = 0;
   out_5593547350807029480[26] = 0;
   out_5593547350807029480[27] = 0;
   out_5593547350807029480[28] = 0;
   out_5593547350807029480[29] = 0;
   out_5593547350807029480[30] = 1;
   out_5593547350807029480[31] = 0;
   out_5593547350807029480[32] = 0;
   out_5593547350807029480[33] = 0;
   out_5593547350807029480[34] = 0;
   out_5593547350807029480[35] = 0;
   out_5593547350807029480[36] = 0;
   out_5593547350807029480[37] = 0;
   out_5593547350807029480[38] = 0;
   out_5593547350807029480[39] = 0;
   out_5593547350807029480[40] = 1;
   out_5593547350807029480[41] = 0;
   out_5593547350807029480[42] = 0;
   out_5593547350807029480[43] = 0;
   out_5593547350807029480[44] = 0;
   out_5593547350807029480[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5593547350807029480[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5593547350807029480[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5593547350807029480[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5593547350807029480[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5593547350807029480[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5593547350807029480[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5593547350807029480[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5593547350807029480[53] = -9.8000000000000007*dt;
   out_5593547350807029480[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5593547350807029480[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5593547350807029480[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5593547350807029480[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5593547350807029480[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5593547350807029480[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5593547350807029480[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5593547350807029480[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5593547350807029480[62] = 0;
   out_5593547350807029480[63] = 0;
   out_5593547350807029480[64] = 0;
   out_5593547350807029480[65] = 0;
   out_5593547350807029480[66] = 0;
   out_5593547350807029480[67] = 0;
   out_5593547350807029480[68] = 0;
   out_5593547350807029480[69] = 0;
   out_5593547350807029480[70] = 1;
   out_5593547350807029480[71] = 0;
   out_5593547350807029480[72] = 0;
   out_5593547350807029480[73] = 0;
   out_5593547350807029480[74] = 0;
   out_5593547350807029480[75] = 0;
   out_5593547350807029480[76] = 0;
   out_5593547350807029480[77] = 0;
   out_5593547350807029480[78] = 0;
   out_5593547350807029480[79] = 0;
   out_5593547350807029480[80] = 1;
}
void h_25(double *state, double *unused, double *out_1879777114914353005) {
   out_1879777114914353005[0] = state[6];
}
void H_25(double *state, double *unused, double *out_456538637935514216) {
   out_456538637935514216[0] = 0;
   out_456538637935514216[1] = 0;
   out_456538637935514216[2] = 0;
   out_456538637935514216[3] = 0;
   out_456538637935514216[4] = 0;
   out_456538637935514216[5] = 0;
   out_456538637935514216[6] = 1;
   out_456538637935514216[7] = 0;
   out_456538637935514216[8] = 0;
}
void h_24(double *state, double *unused, double *out_1062153024039501540) {
   out_1062153024039501540[0] = state[4];
   out_1062153024039501540[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4014464236737837352) {
   out_4014464236737837352[0] = 0;
   out_4014464236737837352[1] = 0;
   out_4014464236737837352[2] = 0;
   out_4014464236737837352[3] = 0;
   out_4014464236737837352[4] = 1;
   out_4014464236737837352[5] = 0;
   out_4014464236737837352[6] = 0;
   out_4014464236737837352[7] = 0;
   out_4014464236737837352[8] = 0;
   out_4014464236737837352[9] = 0;
   out_4014464236737837352[10] = 0;
   out_4014464236737837352[11] = 0;
   out_4014464236737837352[12] = 0;
   out_4014464236737837352[13] = 0;
   out_4014464236737837352[14] = 1;
   out_4014464236737837352[15] = 0;
   out_4014464236737837352[16] = 0;
   out_4014464236737837352[17] = 0;
}
void h_30(double *state, double *unused, double *out_3281741573064809973) {
   out_3281741573064809973[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7373228979427130971) {
   out_7373228979427130971[0] = 0;
   out_7373228979427130971[1] = 0;
   out_7373228979427130971[2] = 0;
   out_7373228979427130971[3] = 0;
   out_7373228979427130971[4] = 1;
   out_7373228979427130971[5] = 0;
   out_7373228979427130971[6] = 0;
   out_7373228979427130971[7] = 0;
   out_7373228979427130971[8] = 0;
}
void h_26(double *state, double *unused, double *out_999822268178850535) {
   out_999822268178850535[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3284964680938542008) {
   out_3284964680938542008[0] = 0;
   out_3284964680938542008[1] = 0;
   out_3284964680938542008[2] = 0;
   out_3284964680938542008[3] = 0;
   out_3284964680938542008[4] = 0;
   out_3284964680938542008[5] = 0;
   out_3284964680938542008[6] = 0;
   out_3284964680938542008[7] = 1;
   out_3284964680938542008[8] = 0;
}
void h_27(double *state, double *unused, double *out_6842157112595646237) {
   out_6842157112595646237[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5198465667626706060) {
   out_5198465667626706060[0] = 0;
   out_5198465667626706060[1] = 0;
   out_5198465667626706060[2] = 0;
   out_5198465667626706060[3] = 1;
   out_5198465667626706060[4] = 0;
   out_5198465667626706060[5] = 0;
   out_5198465667626706060[6] = 0;
   out_5198465667626706060[7] = 0;
   out_5198465667626706060[8] = 0;
}
void h_29(double *state, double *unused, double *out_6759345681000514990) {
   out_6759345681000514990[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7883460323741523155) {
   out_7883460323741523155[0] = 0;
   out_7883460323741523155[1] = 1;
   out_7883460323741523155[2] = 0;
   out_7883460323741523155[3] = 0;
   out_7883460323741523155[4] = 0;
   out_7883460323741523155[5] = 0;
   out_7883460323741523155[6] = 0;
   out_7883460323741523155[7] = 0;
   out_7883460323741523155[8] = 0;
}
void h_28(double *state, double *unused, double *out_7547257068044724173) {
   out_7547257068044724173[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1597296076312375547) {
   out_1597296076312375547[0] = 1;
   out_1597296076312375547[1] = 0;
   out_1597296076312375547[2] = 0;
   out_1597296076312375547[3] = 0;
   out_1597296076312375547[4] = 0;
   out_1597296076312375547[5] = 0;
   out_1597296076312375547[6] = 0;
   out_1597296076312375547[7] = 0;
   out_1597296076312375547[8] = 0;
}
void h_31(double *state, double *unused, double *out_4080565123435819441) {
   out_4080565123435819441[0] = state[8];
}
void H_31(double *state, double *unused, double *out_487184599812474644) {
   out_487184599812474644[0] = 0;
   out_487184599812474644[1] = 0;
   out_487184599812474644[2] = 0;
   out_487184599812474644[3] = 0;
   out_487184599812474644[4] = 0;
   out_487184599812474644[5] = 0;
   out_487184599812474644[6] = 0;
   out_487184599812474644[7] = 0;
   out_487184599812474644[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8535675879233886937) {
  err_fun(nom_x, delta_x, out_8535675879233886937);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6876297271512678565) {
  inv_err_fun(nom_x, true_x, out_6876297271512678565);
}
void car_H_mod_fun(double *state, double *out_3429589102358787544) {
  H_mod_fun(state, out_3429589102358787544);
}
void car_f_fun(double *state, double dt, double *out_557623254591872103) {
  f_fun(state,  dt, out_557623254591872103);
}
void car_F_fun(double *state, double dt, double *out_5593547350807029480) {
  F_fun(state,  dt, out_5593547350807029480);
}
void car_h_25(double *state, double *unused, double *out_1879777114914353005) {
  h_25(state, unused, out_1879777114914353005);
}
void car_H_25(double *state, double *unused, double *out_456538637935514216) {
  H_25(state, unused, out_456538637935514216);
}
void car_h_24(double *state, double *unused, double *out_1062153024039501540) {
  h_24(state, unused, out_1062153024039501540);
}
void car_H_24(double *state, double *unused, double *out_4014464236737837352) {
  H_24(state, unused, out_4014464236737837352);
}
void car_h_30(double *state, double *unused, double *out_3281741573064809973) {
  h_30(state, unused, out_3281741573064809973);
}
void car_H_30(double *state, double *unused, double *out_7373228979427130971) {
  H_30(state, unused, out_7373228979427130971);
}
void car_h_26(double *state, double *unused, double *out_999822268178850535) {
  h_26(state, unused, out_999822268178850535);
}
void car_H_26(double *state, double *unused, double *out_3284964680938542008) {
  H_26(state, unused, out_3284964680938542008);
}
void car_h_27(double *state, double *unused, double *out_6842157112595646237) {
  h_27(state, unused, out_6842157112595646237);
}
void car_H_27(double *state, double *unused, double *out_5198465667626706060) {
  H_27(state, unused, out_5198465667626706060);
}
void car_h_29(double *state, double *unused, double *out_6759345681000514990) {
  h_29(state, unused, out_6759345681000514990);
}
void car_H_29(double *state, double *unused, double *out_7883460323741523155) {
  H_29(state, unused, out_7883460323741523155);
}
void car_h_28(double *state, double *unused, double *out_7547257068044724173) {
  h_28(state, unused, out_7547257068044724173);
}
void car_H_28(double *state, double *unused, double *out_1597296076312375547) {
  H_28(state, unused, out_1597296076312375547);
}
void car_h_31(double *state, double *unused, double *out_4080565123435819441) {
  h_31(state, unused, out_4080565123435819441);
}
void car_H_31(double *state, double *unused, double *out_487184599812474644) {
  H_31(state, unused, out_487184599812474644);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
