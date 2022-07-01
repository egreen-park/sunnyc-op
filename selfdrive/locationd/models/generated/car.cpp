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
void err_fun(double *nom_x, double *delta_x, double *out_1714716606689705258) {
   out_1714716606689705258[0] = delta_x[0] + nom_x[0];
   out_1714716606689705258[1] = delta_x[1] + nom_x[1];
   out_1714716606689705258[2] = delta_x[2] + nom_x[2];
   out_1714716606689705258[3] = delta_x[3] + nom_x[3];
   out_1714716606689705258[4] = delta_x[4] + nom_x[4];
   out_1714716606689705258[5] = delta_x[5] + nom_x[5];
   out_1714716606689705258[6] = delta_x[6] + nom_x[6];
   out_1714716606689705258[7] = delta_x[7] + nom_x[7];
   out_1714716606689705258[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8901916672351263566) {
   out_8901916672351263566[0] = -nom_x[0] + true_x[0];
   out_8901916672351263566[1] = -nom_x[1] + true_x[1];
   out_8901916672351263566[2] = -nom_x[2] + true_x[2];
   out_8901916672351263566[3] = -nom_x[3] + true_x[3];
   out_8901916672351263566[4] = -nom_x[4] + true_x[4];
   out_8901916672351263566[5] = -nom_x[5] + true_x[5];
   out_8901916672351263566[6] = -nom_x[6] + true_x[6];
   out_8901916672351263566[7] = -nom_x[7] + true_x[7];
   out_8901916672351263566[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_2185529571574228716) {
   out_2185529571574228716[0] = 1.0;
   out_2185529571574228716[1] = 0;
   out_2185529571574228716[2] = 0;
   out_2185529571574228716[3] = 0;
   out_2185529571574228716[4] = 0;
   out_2185529571574228716[5] = 0;
   out_2185529571574228716[6] = 0;
   out_2185529571574228716[7] = 0;
   out_2185529571574228716[8] = 0;
   out_2185529571574228716[9] = 0;
   out_2185529571574228716[10] = 1.0;
   out_2185529571574228716[11] = 0;
   out_2185529571574228716[12] = 0;
   out_2185529571574228716[13] = 0;
   out_2185529571574228716[14] = 0;
   out_2185529571574228716[15] = 0;
   out_2185529571574228716[16] = 0;
   out_2185529571574228716[17] = 0;
   out_2185529571574228716[18] = 0;
   out_2185529571574228716[19] = 0;
   out_2185529571574228716[20] = 1.0;
   out_2185529571574228716[21] = 0;
   out_2185529571574228716[22] = 0;
   out_2185529571574228716[23] = 0;
   out_2185529571574228716[24] = 0;
   out_2185529571574228716[25] = 0;
   out_2185529571574228716[26] = 0;
   out_2185529571574228716[27] = 0;
   out_2185529571574228716[28] = 0;
   out_2185529571574228716[29] = 0;
   out_2185529571574228716[30] = 1.0;
   out_2185529571574228716[31] = 0;
   out_2185529571574228716[32] = 0;
   out_2185529571574228716[33] = 0;
   out_2185529571574228716[34] = 0;
   out_2185529571574228716[35] = 0;
   out_2185529571574228716[36] = 0;
   out_2185529571574228716[37] = 0;
   out_2185529571574228716[38] = 0;
   out_2185529571574228716[39] = 0;
   out_2185529571574228716[40] = 1.0;
   out_2185529571574228716[41] = 0;
   out_2185529571574228716[42] = 0;
   out_2185529571574228716[43] = 0;
   out_2185529571574228716[44] = 0;
   out_2185529571574228716[45] = 0;
   out_2185529571574228716[46] = 0;
   out_2185529571574228716[47] = 0;
   out_2185529571574228716[48] = 0;
   out_2185529571574228716[49] = 0;
   out_2185529571574228716[50] = 1.0;
   out_2185529571574228716[51] = 0;
   out_2185529571574228716[52] = 0;
   out_2185529571574228716[53] = 0;
   out_2185529571574228716[54] = 0;
   out_2185529571574228716[55] = 0;
   out_2185529571574228716[56] = 0;
   out_2185529571574228716[57] = 0;
   out_2185529571574228716[58] = 0;
   out_2185529571574228716[59] = 0;
   out_2185529571574228716[60] = 1.0;
   out_2185529571574228716[61] = 0;
   out_2185529571574228716[62] = 0;
   out_2185529571574228716[63] = 0;
   out_2185529571574228716[64] = 0;
   out_2185529571574228716[65] = 0;
   out_2185529571574228716[66] = 0;
   out_2185529571574228716[67] = 0;
   out_2185529571574228716[68] = 0;
   out_2185529571574228716[69] = 0;
   out_2185529571574228716[70] = 1.0;
   out_2185529571574228716[71] = 0;
   out_2185529571574228716[72] = 0;
   out_2185529571574228716[73] = 0;
   out_2185529571574228716[74] = 0;
   out_2185529571574228716[75] = 0;
   out_2185529571574228716[76] = 0;
   out_2185529571574228716[77] = 0;
   out_2185529571574228716[78] = 0;
   out_2185529571574228716[79] = 0;
   out_2185529571574228716[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5480925988431591271) {
   out_5480925988431591271[0] = state[0];
   out_5480925988431591271[1] = state[1];
   out_5480925988431591271[2] = state[2];
   out_5480925988431591271[3] = state[3];
   out_5480925988431591271[4] = state[4];
   out_5480925988431591271[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5480925988431591271[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5480925988431591271[7] = state[7];
   out_5480925988431591271[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3690847359028175658) {
   out_3690847359028175658[0] = 1;
   out_3690847359028175658[1] = 0;
   out_3690847359028175658[2] = 0;
   out_3690847359028175658[3] = 0;
   out_3690847359028175658[4] = 0;
   out_3690847359028175658[5] = 0;
   out_3690847359028175658[6] = 0;
   out_3690847359028175658[7] = 0;
   out_3690847359028175658[8] = 0;
   out_3690847359028175658[9] = 0;
   out_3690847359028175658[10] = 1;
   out_3690847359028175658[11] = 0;
   out_3690847359028175658[12] = 0;
   out_3690847359028175658[13] = 0;
   out_3690847359028175658[14] = 0;
   out_3690847359028175658[15] = 0;
   out_3690847359028175658[16] = 0;
   out_3690847359028175658[17] = 0;
   out_3690847359028175658[18] = 0;
   out_3690847359028175658[19] = 0;
   out_3690847359028175658[20] = 1;
   out_3690847359028175658[21] = 0;
   out_3690847359028175658[22] = 0;
   out_3690847359028175658[23] = 0;
   out_3690847359028175658[24] = 0;
   out_3690847359028175658[25] = 0;
   out_3690847359028175658[26] = 0;
   out_3690847359028175658[27] = 0;
   out_3690847359028175658[28] = 0;
   out_3690847359028175658[29] = 0;
   out_3690847359028175658[30] = 1;
   out_3690847359028175658[31] = 0;
   out_3690847359028175658[32] = 0;
   out_3690847359028175658[33] = 0;
   out_3690847359028175658[34] = 0;
   out_3690847359028175658[35] = 0;
   out_3690847359028175658[36] = 0;
   out_3690847359028175658[37] = 0;
   out_3690847359028175658[38] = 0;
   out_3690847359028175658[39] = 0;
   out_3690847359028175658[40] = 1;
   out_3690847359028175658[41] = 0;
   out_3690847359028175658[42] = 0;
   out_3690847359028175658[43] = 0;
   out_3690847359028175658[44] = 0;
   out_3690847359028175658[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3690847359028175658[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3690847359028175658[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3690847359028175658[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3690847359028175658[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3690847359028175658[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3690847359028175658[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3690847359028175658[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3690847359028175658[53] = -9.8000000000000007*dt;
   out_3690847359028175658[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3690847359028175658[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3690847359028175658[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3690847359028175658[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3690847359028175658[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3690847359028175658[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3690847359028175658[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3690847359028175658[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3690847359028175658[62] = 0;
   out_3690847359028175658[63] = 0;
   out_3690847359028175658[64] = 0;
   out_3690847359028175658[65] = 0;
   out_3690847359028175658[66] = 0;
   out_3690847359028175658[67] = 0;
   out_3690847359028175658[68] = 0;
   out_3690847359028175658[69] = 0;
   out_3690847359028175658[70] = 1;
   out_3690847359028175658[71] = 0;
   out_3690847359028175658[72] = 0;
   out_3690847359028175658[73] = 0;
   out_3690847359028175658[74] = 0;
   out_3690847359028175658[75] = 0;
   out_3690847359028175658[76] = 0;
   out_3690847359028175658[77] = 0;
   out_3690847359028175658[78] = 0;
   out_3690847359028175658[79] = 0;
   out_3690847359028175658[80] = 1;
}
void h_25(double *state, double *unused, double *out_4444842497576461069) {
   out_4444842497576461069[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4527421359871902241) {
   out_4527421359871902241[0] = 0;
   out_4527421359871902241[1] = 0;
   out_4527421359871902241[2] = 0;
   out_4527421359871902241[3] = 0;
   out_4527421359871902241[4] = 0;
   out_4527421359871902241[5] = 0;
   out_4527421359871902241[6] = 1;
   out_4527421359871902241[7] = 0;
   out_4527421359871902241[8] = 0;
}
void h_24(double *state, double *unused, double *out_5924539727079180738) {
   out_5924539727079180738[0] = state[4];
   out_5924539727079180738[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5144430306025478481) {
   out_5144430306025478481[0] = 0;
   out_5144430306025478481[1] = 0;
   out_5144430306025478481[2] = 0;
   out_5144430306025478481[3] = 0;
   out_5144430306025478481[4] = 1;
   out_5144430306025478481[5] = 0;
   out_5144430306025478481[6] = 0;
   out_5144430306025478481[7] = 0;
   out_5144430306025478481[8] = 0;
   out_5144430306025478481[9] = 0;
   out_5144430306025478481[10] = 0;
   out_5144430306025478481[11] = 0;
   out_5144430306025478481[12] = 0;
   out_5144430306025478481[13] = 0;
   out_5144430306025478481[14] = 1;
   out_5144430306025478481[15] = 0;
   out_5144430306025478481[16] = 0;
   out_5144430306025478481[17] = 0;
}
void h_30(double *state, double *unused, double *out_4692874787107390113) {
   out_4692874787107390113[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4656760307015142311) {
   out_4656760307015142311[0] = 0;
   out_4656760307015142311[1] = 0;
   out_4656760307015142311[2] = 0;
   out_4656760307015142311[3] = 0;
   out_4656760307015142311[4] = 1;
   out_4656760307015142311[5] = 0;
   out_4656760307015142311[6] = 0;
   out_4656760307015142311[7] = 0;
   out_4656760307015142311[8] = 0;
}
void h_26(double *state, double *unused, double *out_9196708246802926460) {
   out_9196708246802926460[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8268924678745958465) {
   out_8268924678745958465[0] = 0;
   out_8268924678745958465[1] = 0;
   out_8268924678745958465[2] = 0;
   out_8268924678745958465[3] = 0;
   out_8268924678745958465[4] = 0;
   out_8268924678745958465[5] = 0;
   out_8268924678745958465[6] = 0;
   out_8268924678745958465[7] = 1;
   out_8268924678745958465[8] = 0;
}
void h_27(double *state, double *unused, double *out_387512604507279863) {
   out_387512604507279863[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6831523618815567222) {
   out_6831523618815567222[0] = 0;
   out_6831523618815567222[1] = 0;
   out_6831523618815567222[2] = 0;
   out_6831523618815567222[3] = 1;
   out_6831523618815567222[4] = 0;
   out_6831523618815567222[5] = 0;
   out_6831523618815567222[6] = 0;
   out_6831523618815567222[7] = 0;
   out_6831523618815567222[8] = 0;
}
void h_29(double *state, double *unused, double *out_5138544832043743250) {
   out_5138544832043743250[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4146528962700750127) {
   out_4146528962700750127[0] = 0;
   out_4146528962700750127[1] = 1;
   out_4146528962700750127[2] = 0;
   out_4146528962700750127[3] = 0;
   out_4146528962700750127[4] = 0;
   out_4146528962700750127[5] = 0;
   out_4146528962700750127[6] = 0;
   out_4146528962700750127[7] = 0;
   out_4146528962700750127[8] = 0;
}
void h_28(double *state, double *unused, double *out_6030242475588236698) {
   out_6030242475588236698[0] = state[0];
}
void H_28(double *state, double *unused, double *out_9217816093939270915) {
   out_9217816093939270915[0] = 1;
   out_9217816093939270915[1] = 0;
   out_9217816093939270915[2] = 0;
   out_9217816093939270915[3] = 0;
   out_9217816093939270915[4] = 0;
   out_9217816093939270915[5] = 0;
   out_9217816093939270915[6] = 0;
   out_9217816093939270915[7] = 0;
   out_9217816093939270915[8] = 0;
}
void h_31(double *state, double *unused, double *out_2924870645492627357) {
   out_2924870645492627357[0] = state[8];
}
void H_31(double *state, double *unused, double *out_4496775397994941813) {
   out_4496775397994941813[0] = 0;
   out_4496775397994941813[1] = 0;
   out_4496775397994941813[2] = 0;
   out_4496775397994941813[3] = 0;
   out_4496775397994941813[4] = 0;
   out_4496775397994941813[5] = 0;
   out_4496775397994941813[6] = 0;
   out_4496775397994941813[7] = 0;
   out_4496775397994941813[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_1714716606689705258) {
  err_fun(nom_x, delta_x, out_1714716606689705258);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8901916672351263566) {
  inv_err_fun(nom_x, true_x, out_8901916672351263566);
}
void car_H_mod_fun(double *state, double *out_2185529571574228716) {
  H_mod_fun(state, out_2185529571574228716);
}
void car_f_fun(double *state, double dt, double *out_5480925988431591271) {
  f_fun(state,  dt, out_5480925988431591271);
}
void car_F_fun(double *state, double dt, double *out_3690847359028175658) {
  F_fun(state,  dt, out_3690847359028175658);
}
void car_h_25(double *state, double *unused, double *out_4444842497576461069) {
  h_25(state, unused, out_4444842497576461069);
}
void car_H_25(double *state, double *unused, double *out_4527421359871902241) {
  H_25(state, unused, out_4527421359871902241);
}
void car_h_24(double *state, double *unused, double *out_5924539727079180738) {
  h_24(state, unused, out_5924539727079180738);
}
void car_H_24(double *state, double *unused, double *out_5144430306025478481) {
  H_24(state, unused, out_5144430306025478481);
}
void car_h_30(double *state, double *unused, double *out_4692874787107390113) {
  h_30(state, unused, out_4692874787107390113);
}
void car_H_30(double *state, double *unused, double *out_4656760307015142311) {
  H_30(state, unused, out_4656760307015142311);
}
void car_h_26(double *state, double *unused, double *out_9196708246802926460) {
  h_26(state, unused, out_9196708246802926460);
}
void car_H_26(double *state, double *unused, double *out_8268924678745958465) {
  H_26(state, unused, out_8268924678745958465);
}
void car_h_27(double *state, double *unused, double *out_387512604507279863) {
  h_27(state, unused, out_387512604507279863);
}
void car_H_27(double *state, double *unused, double *out_6831523618815567222) {
  H_27(state, unused, out_6831523618815567222);
}
void car_h_29(double *state, double *unused, double *out_5138544832043743250) {
  h_29(state, unused, out_5138544832043743250);
}
void car_H_29(double *state, double *unused, double *out_4146528962700750127) {
  H_29(state, unused, out_4146528962700750127);
}
void car_h_28(double *state, double *unused, double *out_6030242475588236698) {
  h_28(state, unused, out_6030242475588236698);
}
void car_H_28(double *state, double *unused, double *out_9217816093939270915) {
  H_28(state, unused, out_9217816093939270915);
}
void car_h_31(double *state, double *unused, double *out_2924870645492627357) {
  h_31(state, unused, out_2924870645492627357);
}
void car_H_31(double *state, double *unused, double *out_4496775397994941813) {
  H_31(state, unused, out_4496775397994941813);
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
