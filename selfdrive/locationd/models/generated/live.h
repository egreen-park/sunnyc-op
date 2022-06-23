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
void live_H(double *in_vec, double *out_9094854505137700114);
void live_err_fun(double *nom_x, double *delta_x, double *out_7189592963740075634);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7508967909468565488);
void live_H_mod_fun(double *state, double *out_6406871339729092844);
void live_f_fun(double *state, double dt, double *out_8326757184263707931);
void live_F_fun(double *state, double dt, double *out_6845140873145997354);
void live_h_4(double *state, double *unused, double *out_1404230682439242376);
void live_H_4(double *state, double *unused, double *out_1958599729245028118);
void live_h_9(double *state, double *unused, double *out_4901740455674426920);
void live_H_9(double *state, double *unused, double *out_2198568007109749365);
void live_h_10(double *state, double *unused, double *out_2989697102392970489);
void live_H_10(double *state, double *unused, double *out_1394276942919263448);
void live_h_12(double *state, double *unused, double *out_7092051119784227734);
void live_H_12(double *state, double *unused, double *out_2579698754292621785);
void live_h_31(double *state, double *unused, double *out_5309002211289508523);
void live_H_31(double *state, double *unused, double *out_5325261786617635494);
void live_h_32(double *state, double *unused, double *out_2904342480900863083);
void live_H_32(double *state, double *unused, double *out_984257010884528606);
void live_h_13(double *state, double *unused, double *out_4995199596212226186);
void live_H_13(double *state, double *unused, double *out_206491193835808968);
void live_h_14(double *state, double *unused, double *out_4901740455674426920);
void live_H_14(double *state, double *unused, double *out_2198568007109749365);
void live_h_33(double *state, double *unused, double *out_3600356268595564734);
void live_H_33(double *state, double *unused, double *out_8475818791256493098);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}