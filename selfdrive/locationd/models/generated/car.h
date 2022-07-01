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
void car_err_fun(double *nom_x, double *delta_x, double *out_1714716606689705258);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8901916672351263566);
void car_H_mod_fun(double *state, double *out_2185529571574228716);
void car_f_fun(double *state, double dt, double *out_5480925988431591271);
void car_F_fun(double *state, double dt, double *out_3690847359028175658);
void car_h_25(double *state, double *unused, double *out_4444842497576461069);
void car_H_25(double *state, double *unused, double *out_4527421359871902241);
void car_h_24(double *state, double *unused, double *out_5924539727079180738);
void car_H_24(double *state, double *unused, double *out_5144430306025478481);
void car_h_30(double *state, double *unused, double *out_4692874787107390113);
void car_H_30(double *state, double *unused, double *out_4656760307015142311);
void car_h_26(double *state, double *unused, double *out_9196708246802926460);
void car_H_26(double *state, double *unused, double *out_8268924678745958465);
void car_h_27(double *state, double *unused, double *out_387512604507279863);
void car_H_27(double *state, double *unused, double *out_6831523618815567222);
void car_h_29(double *state, double *unused, double *out_5138544832043743250);
void car_H_29(double *state, double *unused, double *out_4146528962700750127);
void car_h_28(double *state, double *unused, double *out_6030242475588236698);
void car_H_28(double *state, double *unused, double *out_9217816093939270915);
void car_h_31(double *state, double *unused, double *out_2924870645492627357);
void car_H_31(double *state, double *unused, double *out_4496775397994941813);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}