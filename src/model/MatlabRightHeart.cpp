// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "MatlabRightHeart.h"

#include <algorithm>
#include <cmath>

#include "Model.h"

namespace {
constexpr double kPi = 3.14159265358979323846;
}

void MatlabRightHeart::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 5,
                     {"q_tricuspid", "rv_volume", "zeta_pulmonary_valve",
                      "zeta_tricuspid"});
}

double MatlabRightHeart::clamp(double value, double low, double high) {
  return std::max(low, std::min(value, high));
}

double MatlabRightHeart::activation(std::vector<double>& parameters) const {
  const double cardiac_period = model->cardiac_cycle_period;
  const double ts = parameters[global_param_ids[ParamId::TS]];
  const double tr = parameters[global_param_ids[ParamId::TR]];
  const double a_coeff = 2.0 / (kPi + 2.0);
  const double t_mod = std::fmod(model->time, cardiac_period);

  if (t_mod < ts / 2.0) {
    return (1.0 - a_coeff) * 2.0 * t_mod / ts;
  }
  if (t_mod < ts) {
    return 1.0 - a_coeff + a_coeff * std::sin(kPi * (t_mod / ts - 0.5));
  }
  if (t_mod < ts + tr) {
    return 1.0 - std::sin((kPi / 2.0) * (t_mod - ts) / tr);
  }
  return 0.0;
}

void MatlabRightHeart::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  const double p_out = y[global_var_ids[0]];
  const double q_out = y[global_var_ids[1]];
  const double qtc = y[global_var_ids[2]];
  const double vrv = y[global_var_ids[3]];
  const double zeta_pv = y[global_var_ids[4]];
  const double zeta_tc = y[global_var_ids[5]];

  const double ura = parameters[global_param_ids[ParamId::URA]];
  const double fd = parameters[global_param_ids[ParamId::FD]];
  const double fs = parameters[global_param_ids[ParamId::FS]];
  const double rv_slope_d = parameters[global_param_ids[ParamId::RV_SLOPE_D]];
  const double rv_offset_d = parameters[global_param_ids[ParamId::RV_OFFSET_D]];
  const double rv_slope_s = parameters[global_param_ids[ParamId::RV_SLOPE_S]];
  const double rv_offset_s = parameters[global_param_ids[ParamId::RV_OFFSET_S]];
  const double rho = parameters[global_param_ids[ParamId::RHO]];
  const double leff_pv = parameters[global_param_ids[ParamId::LEFF_PV]];
  const double aann_pv = parameters[global_param_ids[ParamId::AANN_PV]];
  const double kvo_pv = parameters[global_param_ids[ParamId::KVO_PV]];
  const double kvc_pv = parameters[global_param_ids[ParamId::KVC_PV]];
  const double dpopen_pv = parameters[global_param_ids[ParamId::DPOPEN_PV]];
  const double dpclose_pv = parameters[global_param_ids[ParamId::DPCLOSE_PV]];
  const double leff_tc = parameters[global_param_ids[ParamId::LEFF_TC]];
  const double aann_tc = parameters[global_param_ids[ParamId::AANN_TC]];
  const double kvo_tc = parameters[global_param_ids[ParamId::KVO_TC]];
  const double kvc_tc = parameters[global_param_ids[ParamId::KVC_TC]];
  const double dpopen_tc = parameters[global_param_ids[ParamId::DPOPEN_TC]];
  const double dpclose_tc = parameters[global_param_ids[ParamId::DPCLOSE_TC]];
  const double p_conv =
      parameters[global_param_ids[ParamId::PRESSURE_CONVERSION]];
  const double zeta_min = parameters[global_param_ids[ParamId::ZETA_MIN]];

  const double p_out_mmhg = p_out / p_conv;
  const double act = activation(parameters);
  const double prv =
      (1.0 - act) * fd * (rv_slope_d * vrv + rv_offset_d) +
      act * fs * (rv_slope_s * vrv + rv_offset_s);

  const double zeta_pv_eff = std::max(zeta_pv, zeta_min);
  const double zeta_tc_eff = std::max(zeta_tc, zeta_min);
  const double aeff_pv = aann_pv * zeta_pv_eff;
  const double aeff_tc = aann_tc * zeta_tc_eff;
  const double b_pv = rho / (2.0 * aeff_pv * aeff_pv);
  const double b_tc = rho / (2.0 * aeff_tc * aeff_tc);
  const double l_pv = rho * leff_pv / aeff_pv;
  const double l_tc = rho * leff_tc / aeff_tc;

  const double delta_p_pv = prv - p_out_mmhg;
  const double delta_p_tc = ura - prv;

  const double dqtc_rhs =
      (aeff_tc < 1e-5) ? 0.0 : (delta_p_tc - b_tc * qtc * std::abs(qtc)) / l_tc;
  const double dqout_rhs =
      (aeff_pv < 1e-5)
          ? 0.0
          : (delta_p_pv - b_pv * q_out * std::abs(q_out)) / l_pv;
  const double dzeta_pv_rhs =
      (delta_p_pv > dpopen_pv)
          ? (1.0 - zeta_pv) * kvo_pv * (delta_p_pv - dpopen_pv)
          : zeta_pv * kvc_pv * (delta_p_pv - dpclose_pv);
  const double dzeta_tc_rhs =
      (delta_p_tc > dpopen_tc)
          ? (1.0 - zeta_tc) * kvo_tc * (delta_p_tc - dpopen_tc)
          : zeta_tc * kvc_tc * (delta_p_tc - dpclose_tc);

  system.E.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[3]) = 1.0;
  system.E.coeffRef(global_eqn_ids[2], global_var_ids[1]) = 1.0;
  system.E.coeffRef(global_eqn_ids[3], global_var_ids[4]) = 1.0;
  system.E.coeffRef(global_eqn_ids[4], global_var_ids[5]) = 1.0;

  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;

  system.C.coeffRef(global_eqn_ids[0]) = -dqtc_rhs;
  system.C.coeffRef(global_eqn_ids[2]) = -dqout_rhs;
  system.C.coeffRef(global_eqn_ids[3]) = -dzeta_pv_rhs;
  system.C.coeffRef(global_eqn_ids[4]) = -dzeta_tc_rhs;

  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[3]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[5]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[0]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[1]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[3]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[0]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[3]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[4]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[3]) = 0.0;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[5]) = 0.0;

  const double dprv_dvrv =
      (1.0 - act) * fd * rv_slope_d + act * fs * rv_slope_s;

  if (aeff_tc >= 1e-5) {
    const double inv_l_tc = 1.0 / l_tc;
    const double abs_qtc = std::abs(qtc);
    const double dl_dzeta_tc =
        -rho * leff_tc / (aann_tc * zeta_tc_eff * zeta_tc_eff);
    const double db_dzeta_tc =
        -rho / (aann_tc * aann_tc * zeta_tc_eff * zeta_tc_eff * zeta_tc_eff);
    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
        inv_l_tc * b_tc * 2.0 * abs_qtc;
    system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[3]) =
        inv_l_tc * dprv_dvrv;
    if (zeta_tc > zeta_min) {
      system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[5]) =
          -(qtc * abs_qtc * db_dzeta_tc * l_tc -
            (delta_p_tc - b_tc * qtc * abs_qtc) * dl_dzeta_tc) /
          (l_tc * l_tc);
    }
  }

  if (aeff_pv >= 1e-5) {
    const double inv_l_pv = 1.0 / l_pv;
    const double abs_qout = std::abs(q_out);
    const double dl_dzeta_pv =
        -rho * leff_pv / (aann_pv * zeta_pv_eff * zeta_pv_eff);
    const double db_dzeta_pv =
        -rho / (aann_pv * aann_pv * zeta_pv_eff * zeta_pv_eff * zeta_pv_eff);
    system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[0]) =
        inv_l_pv / p_conv;
    system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[1]) =
        inv_l_pv * b_pv * 2.0 * abs_qout;
    system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[3]) =
        -inv_l_pv * dprv_dvrv;
    if (zeta_pv > zeta_min) {
      system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[4]) =
          -(q_out * abs_qout * db_dzeta_pv * l_pv -
            (delta_p_pv - b_pv * q_out * abs_qout) * dl_dzeta_pv) /
          (l_pv * l_pv);
    }
  }

  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[0]) =
      (delta_p_pv > dpopen_pv) ? kvo_pv * (1.0 - zeta_pv) / p_conv
                               : kvc_pv * zeta_pv / p_conv;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[3]) =
      (delta_p_pv > dpopen_pv) ? -kvo_pv * (1.0 - zeta_pv) * dprv_dvrv
                               : -kvc_pv * zeta_pv * dprv_dvrv;
  system.dC_dy.coeffRef(global_eqn_ids[3], global_var_ids[4]) =
      (delta_p_pv > dpopen_pv) ? kvo_pv * (delta_p_pv - dpopen_pv)
                               : -kvc_pv * (delta_p_pv - dpclose_pv);

  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[3]) =
      (delta_p_tc > dpopen_tc) ? kvo_tc * (1.0 - zeta_tc) * dprv_dvrv
                               : kvc_tc * zeta_tc * dprv_dvrv;
  system.dC_dy.coeffRef(global_eqn_ids[4], global_var_ids[5]) =
      (delta_p_tc > dpopen_tc) ? kvo_tc * (delta_p_tc - dpopen_tc)
                               : -kvc_tc * (delta_p_tc - dpclose_tc);
}

void MatlabRightHeart::post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) {
  const double zeta_cutoff = clamp(
      model->get_parameter_value(global_param_ids[ParamId::ZETA_CLOSED_THRESHOLD]),
      0.0, 1.0);
  y[global_var_ids[4]] = clamp(y[global_var_ids[4]], 0.0, 1.0);
  y[global_var_ids[5]] = clamp(y[global_var_ids[5]], 0.0, 1.0);
  if (y[global_var_ids[4]] < zeta_cutoff) {
    y[global_var_ids[4]] = 0.0;
  }
  if (y[global_var_ids[5]] < zeta_cutoff) {
    y[global_var_ids[5]] = 0.0;
  }
}
