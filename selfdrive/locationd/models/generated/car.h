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
void car_err_fun(double *nom_x, double *delta_x, double *out_4588065020861109859);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7162724072691955863);
void car_H_mod_fun(double *state, double *out_8163166015551776475);
void car_f_fun(double *state, double dt, double *out_4524055425416236406);
void car_F_fun(double *state, double dt, double *out_4795460408306949169);
void car_h_25(double *state, double *unused, double *out_4290333048114306844);
void car_H_25(double *state, double *unused, double *out_3915141440058088316);
void car_h_24(double *state, double *unused, double *out_6449532279368977895);
void car_H_24(double *state, double *unused, double *out_8783956305085795168);
void car_h_30(double *state, double *unused, double *out_6393145768774572331);
void car_H_30(double *state, double *unused, double *out_8442837770185696514);
void car_h_26(double *state, double *unused, double *out_9016812253490646600);
void car_H_26(double *state, double *unused, double *out_7656644758932144540);
void car_h_27(double *state, double *unused, double *out_5233842901015969279);
void car_H_27(double *state, double *unused, double *out_6219243699001753297);
void car_h_29(double *state, double *unused, double *out_5449635915872909896);
void car_H_29(double *state, double *unused, double *out_7932606425871304330);
void car_h_28(double *state, double *unused, double *out_8270628952460301964);
void car_H_28(double *state, double *unused, double *out_5431738630768716712);
void car_h_31(double *state, double *unused, double *out_1623399978556441529);
void car_H_31(double *state, double *unused, double *out_8282852861165496016);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}