#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_1543763515856100254);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5322435719520533216);
void car_H_mod_fun(double *state, double *out_3799576479346241340);
void car_f_fun(double *state, double dt, double *out_505757690443295894);
void car_F_fun(double *state, double dt, double *out_8408676881083716623);
void car_h_25(double *state, double *unused, double *out_2421938134280317232);
void car_H_25(double *state, double *unused, double *out_2138941943651888845);
void car_h_24(double *state, double *unused, double *out_2538555682541301931);
void car_H_24(double *state, double *unused, double *out_4765841301227063104);
void car_h_30(double *state, double *unused, double *out_5471001360604206127);
void car_H_30(double *state, double *unused, double *out_2009602996508648775);
void car_h_26(double *state, double *unused, double *out_2434215309159873442);
void car_H_26(double *state, double *unused, double *out_1602561375222167379);
void car_h_27(double *state, double *unused, double *out_1478428281378654797);
void car_H_27(double *state, double *unused, double *out_165160315291776136);
void car_h_29(double *state, double *unused, double *out_7633875477219314926);
void car_H_29(double *state, double *unused, double *out_2519834340823040959);
void car_h_28(double *state, double *unused, double *out_1558357770065677888);
void car_H_28(double *state, double *unused, double *out_4483464612388367210);
void car_h_31(double *state, double *unused, double *out_3798439048566696737);
void car_H_31(double *state, double *unused, double *out_2169587905528849273);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}