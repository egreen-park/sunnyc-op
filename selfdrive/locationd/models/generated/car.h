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
void car_err_fun(double *nom_x, double *delta_x, double *out_8535675879233886937);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6876297271512678565);
void car_H_mod_fun(double *state, double *out_3429589102358787544);
void car_f_fun(double *state, double dt, double *out_557623254591872103);
void car_F_fun(double *state, double dt, double *out_5593547350807029480);
void car_h_25(double *state, double *unused, double *out_1879777114914353005);
void car_H_25(double *state, double *unused, double *out_456538637935514216);
void car_h_24(double *state, double *unused, double *out_1062153024039501540);
void car_H_24(double *state, double *unused, double *out_4014464236737837352);
void car_h_30(double *state, double *unused, double *out_3281741573064809973);
void car_H_30(double *state, double *unused, double *out_7373228979427130971);
void car_h_26(double *state, double *unused, double *out_999822268178850535);
void car_H_26(double *state, double *unused, double *out_3284964680938542008);
void car_h_27(double *state, double *unused, double *out_6842157112595646237);
void car_H_27(double *state, double *unused, double *out_5198465667626706060);
void car_h_29(double *state, double *unused, double *out_6759345681000514990);
void car_H_29(double *state, double *unused, double *out_7883460323741523155);
void car_h_28(double *state, double *unused, double *out_7547257068044724173);
void car_H_28(double *state, double *unused, double *out_1597296076312375547);
void car_h_31(double *state, double *unused, double *out_4080565123435819441);
void car_H_31(double *state, double *unused, double *out_487184599812474644);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}