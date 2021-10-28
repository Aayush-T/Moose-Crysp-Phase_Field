//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeStressBase.h"

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticityStateVarRateComponent.h"

/**
 * FiniteStrainUObasedCP_PF_new_mod uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the
 * material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 *
 * Involves 4 different types of user objects that calculates:
 * State variables - update state variable (derive from CrystalPlasticityStateVariable)
 * State variable evolution compoment - individual component of state variable incremental rate
 * (derive from CrystalPlasticityStateVariableEvolutionRateComponent)
 * Slip resistance - calcuate slip resistances (derive from CrystalPlasticitySlipResistances)
 * Slip rates - calcuate flow direction and slip rates (derive from CrystalPlasticitySlipRates)
 */
class FiniteStrainUObasedCP_PF_new_mod : public ComputeStressBase
{
public:
  static InputParameters validParams();

  FiniteStrainUObasedCP_PF_new_mod(const InputParameters & parameters);

protected:
  /**
   * updates the stress at a quadrature point.
   */
  virtual void computeQpStress();

  /**
   * initializes the stateful properties such as
   * stress, plastic deformation gradient, slip system resistances, etc.
   */
  virtual void initQpStatefulProperties();

  /**
   * calls the residual and jacobian functions used in the
   * stress update algorithm.
   */
  virtual void calcResidJacob();

  /**
   * updates the slip system resistances and state variables.
   * override to modify slip system resistance and state variable evolution.
   */
  virtual void updateSlipSystemResistanceAndStateVariable();

  /**
   * set variables for stress and internal variable solve.
   */
  virtual void preSolveQp();

  /**
   * solve stress and internal variables.
   */
  virtual void solveQp();

  /**
   * update stress and internal variable after solve.
   */
  virtual void postSolveQp();

  /**
   * set variables for internal variable solve.
   */
  virtual void preSolveStatevar();

  /**
   * solve internal variables.
   */
  virtual void solveStatevar();

  /**
   * update internal variable after solve.
   */
  virtual void postSolveStatevar();

  /**
   * set variables for stress solve.
   */
  virtual void preSolveStress();

  /**
   * solves for stress, updates plastic deformation gradient.
   */
  virtual void solveStress();

  /**
   * update stress and plastic deformation gradient after solve.
   */
  virtual void postSolveStress();

  /**
   * calculate stress residual.
   */
  virtual void calcResidual();

  /**
   * calculate jacobian.
   */
  virtual void calcJacobian();

  /**
   * updates the slip rates.
   */
  virtual void getSlipRates();

  /**
   * calculate the tangent moduli for preconditioner.
   * Default is the elastic stiffness matrix.
   * Exact jacobian is currently implemented.
   * tan_mod_type can be modified to exact in .i file to turn it on.
   */
  virtual void calcTangentModuli();

  /**
   * calculate the elastic tangent moduli for preconditioner.
   */
  virtual void elasticTangentModuli();

  /**
   * calculate the exact tangent moduli for preconditioner.
   */
  virtual void elastoPlasticTangentModuli();

  /**
   * performs the line search update
   */
  bool lineSearchUpdate(const Real rnorm_prev, const RankTwoTensor);

  /**
   * evaluates convergence of state variables.
   */
  virtual bool isStateVariablesConverged();

  /// User objects that define the slip rate
  std::vector<const CrystalPlasticitySlipRate *> _uo_slip_rates;

  /// User objects that define the slip resistance
  std::vector<const CrystalPlasticitySlipResistance *> _uo_slip_resistances;

  /// User objects that define the state variable
  std::vector<const CrystalPlasticityStateVariable *> _uo_state_vars;

  /// User objects that define the state variable evolution rate component
  std::vector<const CrystalPlasticityStateVarRateComponent *> _uo_state_var_evol_rate_comps;

  /// Slip rates material property
  std::vector<MaterialProperty<std::vector<Real>> *> _mat_prop_slip_rates;

  /// Slip resistance material property
  std::vector<MaterialProperty<std::vector<Real>> *> _mat_prop_slip_resistances;

  /// State variable material property
  std::vector<MaterialProperty<std::vector<Real>> *> _mat_prop_state_vars;

  /// Old state variable material property
  std::vector<const MaterialProperty<std::vector<Real>> *> _mat_prop_state_vars_old;

  /// State variable evolution rate component material property
  std::vector<MaterialProperty<std::vector<Real>> *> _mat_prop_state_var_evol_rate_comps;

  /// Number of slip rate user objects
  unsigned int _num_uo_slip_rates;

  /// Number of slip resistance user objects
  unsigned int _num_uo_slip_resistances;

  /// Number of state variable user objects
  unsigned int _num_uo_state_vars;

  /// Number of state variable evolution rate component user objects
  unsigned int _num_uo_state_var_evol_rate_comps;

  /// Local state variable
  std::vector<std::vector<Real>> _state_vars_old;

  /// Local stored state variable (for sub-stepping)
  std::vector<std::vector<Real>> _state_vars_old_stored;

  /// Local old state variable
  std::vector<std::vector<Real>> _state_vars_prev;

  /// Stress residual equation relative tolerance
  Real _rtol;
  /// Stress residual equation absolute tolerance
  Real _abs_tol;
  /// Internal variable update equation tolerance
  Real _stol;
  /// Residual tolerance when variable value is zero. Default 1e-12.
  Real _zero_tol;

  /// Residual tensor
  RankTwoTensor _resid;

  /// Jacobian tensor
  RankFourTensor _jac;

  /// Maximum number of iterations for stress update
  unsigned int _maxiter;
  /// Maximum number of iterations for internal variable update
  unsigned int _maxiterg;

  /// Type of tangent moduli calculation
  MooseEnum _tan_mod_type;

  /// Maximum number of substep iterations
  unsigned int _max_substep_iter;

  /// Flag to activate line serach
  bool _use_line_search;

  /// Minimum line search step size
  Real _min_lsrch_step;

  /// Line search bisection method tolerance
  Real _lsrch_tol;

  /// Line search bisection method maximum iteration number
  unsigned int _lsrch_max_iter;

  /// Line search method
  MooseEnum _lsrch_method;

  MaterialProperty<RankTwoTensor> & _fp;
  const MaterialProperty<RankTwoTensor> & _fp_old;
  MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<RankTwoTensor> & _pk2_old;
  MaterialProperty<RankTwoTensor> & _lag_e;
  MaterialProperty<RankTwoTensor> & _update_rot;
  const MaterialProperty<RankTwoTensor> & _update_rot_old;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;

  /// Crystal rotation
  const MaterialProperty<RankTwoTensor> & _crysrot;

  RankTwoTensor _dfgrd_tmp;
  RankTwoTensor _fe, _fp_old_inv, _fp_inv;
  DenseVector<Real> _tau;
  std::vector<MaterialProperty<std::vector<RankTwoTensor>> *> _flow_direction;

  /// Flag to check whether convergence is achieved
  bool _err_tol;

  /// Used for substepping; Uniformly divides the increment in deformation gradient
  RankTwoTensor _delta_dfgrd, _dfgrd_tmp_old;
  /// Scales the substepping increment to obtain deformation gradient at a substep iteration
  Real _dfgrd_scale_factor;
  /// Current substep size
  Real _substep_dt;

  /// Coupled order parameter defining the crack
  const VariableValue & _c;

  /// Material property defining crack width, declared elsewhere
  const MaterialProperty<Real> & _l;

  /// Material property defining gc parameter, declared elsewhere
  const MaterialProperty<Real> & _gc;

  /// Material property defining pressure, declared elsewhere
  const MaterialProperty<Real> & _pressure;

  /// Use current value of history variable
  bool _use_current_hist;

  /// Use PETSc's VI (Reduced space active set solvers for variational inequalities based on Newton's method) solver
  bool _use_snes_vi_solver;

  /// History variable that prevents crack healing, declared in this material
  MaterialProperty<Real> & _H;

  /// Old value of history variable
  const MaterialProperty<Real> & _H_old;

  /// material property for fracture energy barrier
  const MaterialProperty<Real> & _barrier;

  /// Material property for elastic energy
  MaterialProperty<Real> & _E;

  /// Derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _dEdc;

  /// Second-order derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _d2Ed2c;

  /// Derivative of stress w.r.t damage variable
  MaterialProperty<RankTwoTensor> & _dstress_dc;

  /// Second-order derivative of elastic energy w.r.t damage variable and strain
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;

  /// Material property for energetic degradation function
  const MaterialProperty<Real> & _D;

  /// Derivative of degradation function w.r.t damage variable
  const MaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const MaterialProperty<Real> & _d2Dd2c;

  /// Material property for damage indicator function
  const MaterialProperty<Real> & _I;

  /// Derivative of damage indicator function w.r.t damage variable
  const MaterialProperty<Real> & _dIdc;

  /// Second-order derivative of damage indicator function w.r.t damage variable
  const MaterialProperty<Real> & _d2Id2c;

  ///elastic strain
  MaterialProperty<RankTwoTensor> & _ee;

  MaterialProperty<RankTwoTensor> & _pk2_new;

  ///plastic strain
  //MaterialProperty<RankTwoTensor> & _ep;

  /**
   * Method to split elastic energy based on stress spectral decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStrainSpectral(Real & F_pos, Real & F_neg, RankTwoTensor & ce);

  const Real _bulk_modulus_ref; // reference bulk modulus for vol/non-vol decomposition
};