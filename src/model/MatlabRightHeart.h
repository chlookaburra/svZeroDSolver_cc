// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/**
 * @file MatlabRightHeart.h
 * @brief Right-heart block patterned after USER2.f, with pulmonary loading supplied externally.
 */
#ifndef SVZERODSOLVER_MODEL_MATLAB_RIGHT_HEART_HPP_
#define SVZERODSOLVER_MODEL_MATLAB_RIGHT_HEART_HPP_

#include "Block.h"
#include "SparseSystem.h"

class MatlabRightHeart : public Block {
 public:
  enum ParamId {
    URA = 0,
    TS = 1,
    TR = 2,
    FD = 3,
    FS = 4,
    RV_SLOPE_D = 5,
    RV_OFFSET_D = 6,
    RV_SLOPE_S = 7,
    RV_OFFSET_S = 8,
    RHO = 9,
    LEFF_PV = 10,
    AANN_PV = 11,
    KVO_PV = 12,
    KVC_PV = 13,
    DPOPEN_PV = 14,
    DPCLOSE_PV = 15,
    LEFF_TC = 16,
    AANN_TC = 17,
    KVO_TC = 18,
    KVC_TC = 19,
    DPOPEN_TC = 20,
    DPCLOSE_TC = 21,
    PRESSURE_CONVERSION = 22,
    ZETA_MIN = 23,
    ZETA_CLOSED_THRESHOLD = 24
  };

  MatlabRightHeart(int id, Model* model)
      : Block(id, model, BlockType::matlab_right_heart, BlockClass::chamber,
              {{"Ura", InputParameter()},
               {"Ts", InputParameter()},
               {"Tr", InputParameter()},
               {"fD", InputParameter()},
               {"fS", InputParameter()},
               {"rv_slope_d", InputParameter()},
               {"rv_offset_d", InputParameter()},
               {"rv_slope_s", InputParameter()},
               {"rv_offset_s", InputParameter()},
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
               {"pressure_conversion", InputParameter(true)},
               {"zeta_min", InputParameter(true)},
               {"zeta_closed_threshold", InputParameter(true)}}) {}

  void setup_dofs(DOFHandler& dofhandler) override;
  void update_solution(
      SparseSystem& system, std::vector<double>& parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) override;
  void post_solve(Eigen::Matrix<double, Eigen::Dynamic, 1>& y) override;

  TripletsContributions num_triplets{6, 5, 10};

 private:
  static double clamp(double value, double low, double high);
  double activation(std::vector<double>& parameters) const;
};

#endif  // SVZERODSOLVER_MODEL_MATLAB_RIGHT_HEART_HPP_
