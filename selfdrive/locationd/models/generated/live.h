#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8673519326481956153);
void live_err_fun(double *nom_x, double *delta_x, double *out_6509520593120827794);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1854403484088584618);
void live_H_mod_fun(double *state, double *out_2341655123536354807);
void live_f_fun(double *state, double dt, double *out_9130093515614849635);
void live_F_fun(double *state, double dt, double *out_1465107715934980664);
void live_h_4(double *state, double *unused, double *out_5550033602791381692);
void live_H_4(double *state, double *unused, double *out_6513682271191925253);
void live_h_9(double *state, double *unused, double *out_1917081015155043783);
void live_H_9(double *state, double *unused, double *out_6754871917821515898);
void live_h_10(double *state, double *unused, double *out_6306219088906799002);
void live_H_10(double *state, double *unused, double *out_6227195161888097549);
void live_h_12(double *state, double *unused, double *out_2284330100859273191);
void live_H_12(double *state, double *unused, double *out_6913605394485664568);
void live_h_31(double *state, double *unused, double *out_1719538927514605940);
void live_H_31(double *state, double *unused, double *out_8566399745145018987);
void live_h_32(double *state, double *unused, double *out_4339813958746689082);
void live_H_32(double *state, double *unused, double *out_2965479824870416273);
void live_h_13(double *state, double *unused, double *out_7800067890510650184);
void live_H_13(double *state, double *unused, double *out_4915980534915030171);
void live_h_14(double *state, double *unused, double *out_1917081015155043783);
void live_H_14(double *state, double *unused, double *out_6754871917821515898);
void live_h_33(double *state, double *unused, double *out_5524091102416707340);
void live_H_33(double *state, double *unused, double *out_5415842740506161383);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}