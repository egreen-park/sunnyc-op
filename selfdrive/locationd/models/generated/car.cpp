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
void err_fun(double *nom_x, double *delta_x, double *out_4588065020861109859) {
   out_4588065020861109859[0] = delta_x[0] + nom_x[0];
   out_4588065020861109859[1] = delta_x[1] + nom_x[1];
   out_4588065020861109859[2] = delta_x[2] + nom_x[2];
   out_4588065020861109859[3] = delta_x[3] + nom_x[3];
   out_4588065020861109859[4] = delta_x[4] + nom_x[4];
   out_4588065020861109859[5] = delta_x[5] + nom_x[5];
   out_4588065020861109859[6] = delta_x[6] + nom_x[6];
   out_4588065020861109859[7] = delta_x[7] + nom_x[7];
   out_4588065020861109859[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7162724072691955863) {
   out_7162724072691955863[0] = -nom_x[0] + true_x[0];
   out_7162724072691955863[1] = -nom_x[1] + true_x[1];
   out_7162724072691955863[2] = -nom_x[2] + true_x[2];
   out_7162724072691955863[3] = -nom_x[3] + true_x[3];
   out_7162724072691955863[4] = -nom_x[4] + true_x[4];
   out_7162724072691955863[5] = -nom_x[5] + true_x[5];
   out_7162724072691955863[6] = -nom_x[6] + true_x[6];
   out_7162724072691955863[7] = -nom_x[7] + true_x[7];
   out_7162724072691955863[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8163166015551776475) {
   out_8163166015551776475[0] = 1.0;
   out_8163166015551776475[1] = 0;
   out_8163166015551776475[2] = 0;
   out_8163166015551776475[3] = 0;
   out_8163166015551776475[4] = 0;
   out_8163166015551776475[5] = 0;
   out_8163166015551776475[6] = 0;
   out_8163166015551776475[7] = 0;
   out_8163166015551776475[8] = 0;
   out_8163166015551776475[9] = 0;
   out_8163166015551776475[10] = 1.0;
   out_8163166015551776475[11] = 0;
   out_8163166015551776475[12] = 0;
   out_8163166015551776475[13] = 0;
   out_8163166015551776475[14] = 0;
   out_8163166015551776475[15] = 0;
   out_8163166015551776475[16] = 0;
   out_8163166015551776475[17] = 0;
   out_8163166015551776475[18] = 0;
   out_8163166015551776475[19] = 0;
   out_8163166015551776475[20] = 1.0;
   out_8163166015551776475[21] = 0;
   out_8163166015551776475[22] = 0;
   out_8163166015551776475[23] = 0;
   out_8163166015551776475[24] = 0;
   out_8163166015551776475[25] = 0;
   out_8163166015551776475[26] = 0;
   out_8163166015551776475[27] = 0;
   out_8163166015551776475[28] = 0;
   out_8163166015551776475[29] = 0;
   out_8163166015551776475[30] = 1.0;
   out_8163166015551776475[31] = 0;
   out_8163166015551776475[32] = 0;
   out_8163166015551776475[33] = 0;
   out_8163166015551776475[34] = 0;
   out_8163166015551776475[35] = 0;
   out_8163166015551776475[36] = 0;
   out_8163166015551776475[37] = 0;
   out_8163166015551776475[38] = 0;
   out_8163166015551776475[39] = 0;
   out_8163166015551776475[40] = 1.0;
   out_8163166015551776475[41] = 0;
   out_8163166015551776475[42] = 0;
   out_8163166015551776475[43] = 0;
   out_8163166015551776475[44] = 0;
   out_8163166015551776475[45] = 0;
   out_8163166015551776475[46] = 0;
   out_8163166015551776475[47] = 0;
   out_8163166015551776475[48] = 0;
   out_8163166015551776475[49] = 0;
   out_8163166015551776475[50] = 1.0;
   out_8163166015551776475[51] = 0;
   out_8163166015551776475[52] = 0;
   out_8163166015551776475[53] = 0;
   out_8163166015551776475[54] = 0;
   out_8163166015551776475[55] = 0;
   out_8163166015551776475[56] = 0;
   out_8163166015551776475[57] = 0;
   out_8163166015551776475[58] = 0;
   out_8163166015551776475[59] = 0;
   out_8163166015551776475[60] = 1.0;
   out_8163166015551776475[61] = 0;
   out_8163166015551776475[62] = 0;
   out_8163166015551776475[63] = 0;
   out_8163166015551776475[64] = 0;
   out_8163166015551776475[65] = 0;
   out_8163166015551776475[66] = 0;
   out_8163166015551776475[67] = 0;
   out_8163166015551776475[68] = 0;
   out_8163166015551776475[69] = 0;
   out_8163166015551776475[70] = 1.0;
   out_8163166015551776475[71] = 0;
   out_8163166015551776475[72] = 0;
   out_8163166015551776475[73] = 0;
   out_8163166015551776475[74] = 0;
   out_8163166015551776475[75] = 0;
   out_8163166015551776475[76] = 0;
   out_8163166015551776475[77] = 0;
   out_8163166015551776475[78] = 0;
   out_8163166015551776475[79] = 0;
   out_8163166015551776475[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4524055425416236406) {
   out_4524055425416236406[0] = state[0];
   out_4524055425416236406[1] = state[1];
   out_4524055425416236406[2] = state[2];
   out_4524055425416236406[3] = state[3];
   out_4524055425416236406[4] = state[4];
   out_4524055425416236406[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4524055425416236406[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4524055425416236406[7] = state[7];
   out_4524055425416236406[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4795460408306949169) {
   out_4795460408306949169[0] = 1;
   out_4795460408306949169[1] = 0;
   out_4795460408306949169[2] = 0;
   out_4795460408306949169[3] = 0;
   out_4795460408306949169[4] = 0;
   out_4795460408306949169[5] = 0;
   out_4795460408306949169[6] = 0;
   out_4795460408306949169[7] = 0;
   out_4795460408306949169[8] = 0;
   out_4795460408306949169[9] = 0;
   out_4795460408306949169[10] = 1;
   out_4795460408306949169[11] = 0;
   out_4795460408306949169[12] = 0;
   out_4795460408306949169[13] = 0;
   out_4795460408306949169[14] = 0;
   out_4795460408306949169[15] = 0;
   out_4795460408306949169[16] = 0;
   out_4795460408306949169[17] = 0;
   out_4795460408306949169[18] = 0;
   out_4795460408306949169[19] = 0;
   out_4795460408306949169[20] = 1;
   out_4795460408306949169[21] = 0;
   out_4795460408306949169[22] = 0;
   out_4795460408306949169[23] = 0;
   out_4795460408306949169[24] = 0;
   out_4795460408306949169[25] = 0;
   out_4795460408306949169[26] = 0;
   out_4795460408306949169[27] = 0;
   out_4795460408306949169[28] = 0;
   out_4795460408306949169[29] = 0;
   out_4795460408306949169[30] = 1;
   out_4795460408306949169[31] = 0;
   out_4795460408306949169[32] = 0;
   out_4795460408306949169[33] = 0;
   out_4795460408306949169[34] = 0;
   out_4795460408306949169[35] = 0;
   out_4795460408306949169[36] = 0;
   out_4795460408306949169[37] = 0;
   out_4795460408306949169[38] = 0;
   out_4795460408306949169[39] = 0;
   out_4795460408306949169[40] = 1;
   out_4795460408306949169[41] = 0;
   out_4795460408306949169[42] = 0;
   out_4795460408306949169[43] = 0;
   out_4795460408306949169[44] = 0;
   out_4795460408306949169[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4795460408306949169[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4795460408306949169[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4795460408306949169[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4795460408306949169[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4795460408306949169[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4795460408306949169[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4795460408306949169[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4795460408306949169[53] = -9.8000000000000007*dt;
   out_4795460408306949169[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4795460408306949169[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4795460408306949169[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4795460408306949169[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4795460408306949169[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4795460408306949169[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4795460408306949169[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4795460408306949169[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4795460408306949169[62] = 0;
   out_4795460408306949169[63] = 0;
   out_4795460408306949169[64] = 0;
   out_4795460408306949169[65] = 0;
   out_4795460408306949169[66] = 0;
   out_4795460408306949169[67] = 0;
   out_4795460408306949169[68] = 0;
   out_4795460408306949169[69] = 0;
   out_4795460408306949169[70] = 1;
   out_4795460408306949169[71] = 0;
   out_4795460408306949169[72] = 0;
   out_4795460408306949169[73] = 0;
   out_4795460408306949169[74] = 0;
   out_4795460408306949169[75] = 0;
   out_4795460408306949169[76] = 0;
   out_4795460408306949169[77] = 0;
   out_4795460408306949169[78] = 0;
   out_4795460408306949169[79] = 0;
   out_4795460408306949169[80] = 1;
}
void h_25(double *state, double *unused, double *out_4290333048114306844) {
   out_4290333048114306844[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3915141440058088316) {
   out_3915141440058088316[0] = 0;
   out_3915141440058088316[1] = 0;
   out_3915141440058088316[2] = 0;
   out_3915141440058088316[3] = 0;
   out_3915141440058088316[4] = 0;
   out_3915141440058088316[5] = 0;
   out_3915141440058088316[6] = 1;
   out_3915141440058088316[7] = 0;
   out_3915141440058088316[8] = 0;
}
void h_24(double *state, double *unused, double *out_6449532279368977895) {
   out_6449532279368977895[0] = state[4];
   out_6449532279368977895[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8783956305085795168) {
   out_8783956305085795168[0] = 0;
   out_8783956305085795168[1] = 0;
   out_8783956305085795168[2] = 0;
   out_8783956305085795168[3] = 0;
   out_8783956305085795168[4] = 1;
   out_8783956305085795168[5] = 0;
   out_8783956305085795168[6] = 0;
   out_8783956305085795168[7] = 0;
   out_8783956305085795168[8] = 0;
   out_8783956305085795168[9] = 0;
   out_8783956305085795168[10] = 0;
   out_8783956305085795168[11] = 0;
   out_8783956305085795168[12] = 0;
   out_8783956305085795168[13] = 0;
   out_8783956305085795168[14] = 1;
   out_8783956305085795168[15] = 0;
   out_8783956305085795168[16] = 0;
   out_8783956305085795168[17] = 0;
}
void h_30(double *state, double *unused, double *out_6393145768774572331) {
   out_6393145768774572331[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8442837770185696514) {
   out_8442837770185696514[0] = 0;
   out_8442837770185696514[1] = 0;
   out_8442837770185696514[2] = 0;
   out_8442837770185696514[3] = 0;
   out_8442837770185696514[4] = 1;
   out_8442837770185696514[5] = 0;
   out_8442837770185696514[6] = 0;
   out_8442837770185696514[7] = 0;
   out_8442837770185696514[8] = 0;
}
void h_26(double *state, double *unused, double *out_9016812253490646600) {
   out_9016812253490646600[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7656644758932144540) {
   out_7656644758932144540[0] = 0;
   out_7656644758932144540[1] = 0;
   out_7656644758932144540[2] = 0;
   out_7656644758932144540[3] = 0;
   out_7656644758932144540[4] = 0;
   out_7656644758932144540[5] = 0;
   out_7656644758932144540[6] = 0;
   out_7656644758932144540[7] = 1;
   out_7656644758932144540[8] = 0;
}
void h_27(double *state, double *unused, double *out_5233842901015969279) {
   out_5233842901015969279[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6219243699001753297) {
   out_6219243699001753297[0] = 0;
   out_6219243699001753297[1] = 0;
   out_6219243699001753297[2] = 0;
   out_6219243699001753297[3] = 1;
   out_6219243699001753297[4] = 0;
   out_6219243699001753297[5] = 0;
   out_6219243699001753297[6] = 0;
   out_6219243699001753297[7] = 0;
   out_6219243699001753297[8] = 0;
}
void h_29(double *state, double *unused, double *out_5449635915872909896) {
   out_5449635915872909896[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7932606425871304330) {
   out_7932606425871304330[0] = 0;
   out_7932606425871304330[1] = 1;
   out_7932606425871304330[2] = 0;
   out_7932606425871304330[3] = 0;
   out_7932606425871304330[4] = 0;
   out_7932606425871304330[5] = 0;
   out_7932606425871304330[6] = 0;
   out_7932606425871304330[7] = 0;
   out_7932606425871304330[8] = 0;
}
void h_28(double *state, double *unused, double *out_8270628952460301964) {
   out_8270628952460301964[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5431738630768716712) {
   out_5431738630768716712[0] = 1;
   out_5431738630768716712[1] = 0;
   out_5431738630768716712[2] = 0;
   out_5431738630768716712[3] = 0;
   out_5431738630768716712[4] = 0;
   out_5431738630768716712[5] = 0;
   out_5431738630768716712[6] = 0;
   out_5431738630768716712[7] = 0;
   out_5431738630768716712[8] = 0;
}
void h_31(double *state, double *unused, double *out_1623399978556441529) {
   out_1623399978556441529[0] = state[8];
}
void H_31(double *state, double *unused, double *out_8282852861165496016) {
   out_8282852861165496016[0] = 0;
   out_8282852861165496016[1] = 0;
   out_8282852861165496016[2] = 0;
   out_8282852861165496016[3] = 0;
   out_8282852861165496016[4] = 0;
   out_8282852861165496016[5] = 0;
   out_8282852861165496016[6] = 0;
   out_8282852861165496016[7] = 0;
   out_8282852861165496016[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_4588065020861109859) {
  err_fun(nom_x, delta_x, out_4588065020861109859);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7162724072691955863) {
  inv_err_fun(nom_x, true_x, out_7162724072691955863);
}
void car_H_mod_fun(double *state, double *out_8163166015551776475) {
  H_mod_fun(state, out_8163166015551776475);
}
void car_f_fun(double *state, double dt, double *out_4524055425416236406) {
  f_fun(state,  dt, out_4524055425416236406);
}
void car_F_fun(double *state, double dt, double *out_4795460408306949169) {
  F_fun(state,  dt, out_4795460408306949169);
}
void car_h_25(double *state, double *unused, double *out_4290333048114306844) {
  h_25(state, unused, out_4290333048114306844);
}
void car_H_25(double *state, double *unused, double *out_3915141440058088316) {
  H_25(state, unused, out_3915141440058088316);
}
void car_h_24(double *state, double *unused, double *out_6449532279368977895) {
  h_24(state, unused, out_6449532279368977895);
}
void car_H_24(double *state, double *unused, double *out_8783956305085795168) {
  H_24(state, unused, out_8783956305085795168);
}
void car_h_30(double *state, double *unused, double *out_6393145768774572331) {
  h_30(state, unused, out_6393145768774572331);
}
void car_H_30(double *state, double *unused, double *out_8442837770185696514) {
  H_30(state, unused, out_8442837770185696514);
}
void car_h_26(double *state, double *unused, double *out_9016812253490646600) {
  h_26(state, unused, out_9016812253490646600);
}
void car_H_26(double *state, double *unused, double *out_7656644758932144540) {
  H_26(state, unused, out_7656644758932144540);
}
void car_h_27(double *state, double *unused, double *out_5233842901015969279) {
  h_27(state, unused, out_5233842901015969279);
}
void car_H_27(double *state, double *unused, double *out_6219243699001753297) {
  H_27(state, unused, out_6219243699001753297);
}
void car_h_29(double *state, double *unused, double *out_5449635915872909896) {
  h_29(state, unused, out_5449635915872909896);
}
void car_H_29(double *state, double *unused, double *out_7932606425871304330) {
  H_29(state, unused, out_7932606425871304330);
}
void car_h_28(double *state, double *unused, double *out_8270628952460301964) {
  h_28(state, unused, out_8270628952460301964);
}
void car_H_28(double *state, double *unused, double *out_5431738630768716712) {
  H_28(state, unused, out_5431738630768716712);
}
void car_h_31(double *state, double *unused, double *out_1623399978556441529) {
  h_31(state, unused, out_1623399978556441529);
}
void car_H_31(double *state, double *unused, double *out_8282852861165496016) {
  H_31(state, unused, out_8282852861165496016);
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
