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
void live_H(double *in_vec, double *out_7089917002068737761);
void live_err_fun(double *nom_x, double *delta_x, double *out_5097557166370017434);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_2500389703739010254);
void live_H_mod_fun(double *state, double *out_8605363687855840588);
void live_f_fun(double *state, double dt, double *out_7196605423571964472);
void live_F_fun(double *state, double dt, double *out_5329286583971009738);
void live_h_4(double *state, double *unused, double *out_4155706609999190354);
void live_H_4(double *state, double *unused, double *out_5339599855239346172);
void live_h_9(double *state, double *unused, double *out_1429887674170277705);
void live_H_9(double *state, double *unused, double *out_1947619080025101298);
void live_h_10(double *state, double *unused, double *out_84592872543830086);
void live_H_10(double *state, double *unused, double *out_8116331784126176044);
void live_h_12(double *state, double *unused, double *out_5853750055083503413);
void live_H_12(double *state, double *unused, double *out_6725885841427472448);
void live_h_31(double *state, double *unused, double *out_6772459046034304281);
void live_H_31(double *state, double *unused, double *out_8975295199957065459);
void live_h_32(double *state, double *unused, double *out_209189202375476947);
void live_H_32(double *state, double *unused, double *out_3173990827749452901);
void live_h_13(double *state, double *unused, double *out_6889798033787065432);
void live_H_13(double *state, double *unused, double *out_535304116827451560);
void live_h_14(double *state, double *unused, double *out_1429887674170277705);
void live_H_14(double *state, double *unused, double *out_1947619080025101298);
void live_h_33(double *state, double *unused, double *out_2022710884799879907);
void live_H_33(double *state, double *unused, double *out_5824738195318207855);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}