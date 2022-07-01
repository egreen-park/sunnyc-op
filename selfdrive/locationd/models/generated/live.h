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
void live_H(double *in_vec, double *out_7152185926377700488);
void live_err_fun(double *nom_x, double *delta_x, double *out_3950646692998680769);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3626463882087283498);
void live_H_mod_fun(double *state, double *out_7313751890375090654);
void live_f_fun(double *state, double dt, double *out_3447414059194549864);
void live_F_fun(double *state, double dt, double *out_5197492253302578706);
void live_h_4(double *state, double *unused, double *out_7014553232075907042);
void live_H_4(double *state, double *unused, double *out_48671540338297278);
void live_h_9(double *state, double *unused, double *out_108663546386961525);
void live_H_9(double *state, double *unused, double *out_289861186967887923);
void live_h_10(double *state, double *unused, double *out_4906287204831309256);
void live_H_10(double *state, double *unused, double *out_2357670500591332152);
void live_h_12(double *state, double *unused, double *out_5667340936158651549);
void live_H_12(double *state, double *unused, double *out_5068127948370259073);
void live_h_31(double *state, double *unused, double *out_1338021097045339986);
void live_H_31(double *state, double *unused, double *out_3415333597710904654);
void live_h_32(double *state, double *unused, double *out_6899534750071656641);
void live_H_32(double *state, double *unused, double *out_4927209557846818286);
void live_h_13(double *state, double *unused, double *out_5121156266088859451);
void live_H_13(double *state, double *unused, double *out_8738892226839846620);
void live_h_14(double *state, double *unused, double *out_108663546386961525);
void live_H_14(double *state, double *unused, double *out_289861186967887923);
void live_h_33(double *state, double *unused, double *out_4129557182372481520);
void live_H_33(double *state, double *unused, double *out_6565890602349762258);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}