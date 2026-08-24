// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file MatlabRightHeartPulmonary.h
 * @brief Monolithic right-heart + pulmonary block matching the reference MATLAB model.
 */
#ifndef SVZERODSOLVER_MODEL_MATLAB_RIGHT_HEART_PULMONARY_HPP_
#define SVZERODSOLVER_MODEL_MATLAB_RIGHT_HEART_PULMONARY_HPP_

#include "Block.h"
#include "SparseSystem.h"

class MatlabRightHeartPulmonary : public Block {
 public:
  enum ParamId {
    URA = 0,
    ULA = 1,
    TS = 2,
    TR = 3,
    FD = 4,
    FS = 5,
    RV_SLOPE_D = 6,
    RV_OFFSET_D = 7,
    RV_SLOPE_S = 8,
    RV_OFFSET_S = 9,
    RPA = 10,
    RP1 = 11,
    CPA1 = 12,
    RD1 = 13,
    RP2 = 14,
    CPA2 = 15,
    RD2 = 16,
    RHO = 17,
    LEFF_PV = 18,
    AANN_PV = 19,
    KVO_PV = 20,
    KVC_PV = 21,
    DPOPEN_PV = 22,
    DPCLOSE_PV = 23,
    LEFF_TC = 24,
    AANN_TC = 25,
    KVO_TC = 26,
    KVC_TC = 27,
    DPOPEN_TC = 28,
    DPCLOSE_TC = 29,
    ZETA_MIN = 30
  };

  MatlabRightHeartPulmonary(int id, Model* model)
      : Block(id, model, BlockType::matlab_right_heart_pulmonary,
              BlockClass::closed_loop,
              {{"Ura", InputParameter()},
               {"Ula", InputParameter()},
               {"Ts", InputParameter()},
               {"Tr", InputParameter()},
               {"fD", InputParameter()},
               {"fS", InputParameter()},
               {"rv_slope_d", InputParameter()},
               {"rv_offset_d", InputParameter()},
               {"rv_slope_s", InputParameter()},
               {"rv_offset_s", InputParameter()},
               {"Rpa", InputParameter()},
               {"Rp1", InputParameter()},
               {"Cpa1", InputParameter()},
               {"Rd1", InputParameter()},
               {"Rp2", InputParameter()},
               {"Cpa2", InputParameter()},
               {"Rd2", InputParameter()},
               {"rho", InputParameter()},
               {"leff_pv", InputParameter()},
               {"Aann_pv", InputParameter()},
               {"Kvo_pv", InputParameter()},
               {"Kvc_pv", InputParameter()},
               {"dPopen_pv", InputParameter()},
               {"dPclose_pv", InputParameter()},
               {"leff_tc", InputParameter()},
               {"Aann_tc", InputParameter()},
               {"Kvo_tc", InputParameter()},
               {"Kvc_tc", InputParameter()},
               {"dPopen_tc", InputParameter()},
               {"dPclose_tc", InputParameter()},
               {"zeta_min", InputParameter(true)}}) {}

  void setup_dofs(DOFHandler& dofhandler) override;
  void update_solution(
      SparseSystem& system, std::vector<double>& parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) override;
  void post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) override;

  TripletsContributions num_triplets{0, 7, 18};

 private:
  static double clamp(double value, double low, double high);
  double activation(std::vector<double>& parameters) const;
};

#endif  // SVZERODSOLVER_MODEL_MATLAB_RIGHT_HEART_PULMONARY_HPP_
