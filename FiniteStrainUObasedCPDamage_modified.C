// Nicolò Grilli
// University of Bristol
// 26 Settembre 2021

#include "FiniteStrainUObasedCPDamage_modified.h"

// Are these needed?
#include "petscblaslapack.h"
#include "MooseException.h"
#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticityStateVarRateComponent.h"
#include "MathUtils.h"

registerMooseObject("TensorMechanicsApp", FiniteStrainUObasedCPDamage_modified);

InputParameters
FiniteStrainUObasedCPDamage_modified::validParams()
{
  InputParameters params = FiniteStrainUObasedCP::validParams();
  params.addClassDescription("UserObject based Crystal Plasticity system with damage.");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  return params;
}

FiniteStrainUObasedCPDamage_modified::FiniteStrainUObasedCPDamage_modified(const InputParameters & parameters)
  : FiniteStrainUObasedCP(parameters),
    _c(coupledValue("c")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _H(declareProperty<Real>("hist")), // History variable to avoid damage decrease
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"),
                                          getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
	_D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>(
        "D_name", getVar("c", 0)->name(), getVar("c", 0)->name()))
{
  _err_tol = false;

  _delta_dfgrd.zero();

  // resize the material properties for each userobject
  _mat_prop_slip_rates.resize(_num_uo_slip_rates);
  _mat_prop_slip_resistances.resize(_num_uo_slip_resistances);
  _mat_prop_state_vars.resize(_num_uo_state_vars);
  _mat_prop_state_vars_old.resize(_num_uo_state_vars);
  _mat_prop_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // resize the flow direction
  _flow_direction.resize(_num_uo_slip_rates);

  // resize local state variables
  _state_vars_old.resize(_num_uo_state_vars);
  _state_vars_old_stored.resize(_num_uo_state_vars);
  _state_vars_prev.resize(_num_uo_state_vars);

  // resize user objects
  _uo_slip_rates.resize(_num_uo_slip_rates);
  _uo_slip_resistances.resize(_num_uo_slip_resistances);
  _uo_state_vars.resize(_num_uo_state_vars);
  _uo_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // assign the user objects
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
  {
    _uo_slip_rates[i] = &getUserObjectByName<CrystalPlasticitySlipRate>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_rates")[i]);
    _mat_prop_slip_rates[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_rates")[i]);
    _flow_direction[i] = &declareProperty<std::vector<RankTwoTensor>>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_rates")[i] + "_flow_direction");
  }

  for (unsigned int i = 0; i < _num_uo_slip_resistances; ++i)
  {
    _uo_slip_resistances[i] = &getUserObjectByName<CrystalPlasticitySlipResistance>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_resistances")[i]);
    _mat_prop_slip_resistances[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_slip_resistances")[i]);
  }

  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
  {
    _uo_state_vars[i] = &getUserObjectByName<CrystalPlasticityStateVariable>(
        parameters.get<std::vector<UserObjectName>>("uo_state_vars")[i]);
    _mat_prop_state_vars[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_state_vars")[i]);
    _mat_prop_state_vars_old[i] = &getMaterialPropertyOld<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_state_vars")[i]);
  }

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rate_comps; ++i)
  {
    _uo_state_var_evol_rate_comps[i] = &getUserObjectByName<CrystalPlasticityStateVarRateComponent>(
        parameters.get<std::vector<UserObjectName>>("uo_state_var_evol_rate_comps")[i]);
    _mat_prop_state_var_evol_rate_comps[i] = &declareProperty<std::vector<Real>>(
        parameters.get<std::vector<UserObjectName>>("uo_state_var_evol_rate_comps")[i]);
  }

  _substep_dt = 0.0;
}

void
FiniteStrainUObasedCPDamage_modified::calcResidual()
{
  RankTwoTensor iden(RankTwoTensor::initIdentity), ce, ee, ce_pk2, eqv_slip_incr, pk2_new;

  Real F_pos, F_neg; // tensile and compressive part of the elastic strain energy

  getSlipRates();
  if (_err_tol)
    return;

  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
    for (unsigned int j = 0; j < _uo_slip_rates[i]->variableSize(); ++j)
      eqv_slip_incr +=
          (*_flow_direction[i])[_qp][j] * (*_mat_prop_slip_rates[i])[_qp][j] * _substep_dt;

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  // Decompose ee into positive and negative eigenvalues
  // and calculate elastic energy and stress
  computeStrainSpectral(F_pos, F_neg, ee, pk2_new);

  // calculate history variable and
  // assign elastic free energy to _E
  // for the fracture model
  computeHistoryVariable(F_pos, F_neg);

  // Anisotropic undamaged
  // pk2_new = _elasticity_tensor[_qp] * ee;

  _resid = _pk2[_qp] - pk2_new;
}

void
FiniteStrainUObasedCPDamage_modified::computeStrainSpectral(Real & F_pos, Real & F_neg,
                                                   RankTwoTensor & ee, RankTwoTensor & pk2_new)
{
  //Anisotropic elasticity
  Real Kb = 0.0;    //Kb is reference bulk modulus
  Real Je, delta;   //Je is relative volume,
//  delta is the trace of volumetric part of elastic G-L strain
  RankTwoTensor ident(RankTwoTensor::initIdentity), eqv_slip_incr;
  Real ap_vol, an_vol, ap_cpl;
  // ap_vol is positive volumetric free Energy
  // an_vol is negative volumetric free Energy
  // ap_cpl is the diﬀerence between the speciﬁc free energy used
  //in linear thermoelasticity and the one of an isotropic linear
  //elastic material with bulk modulus
  
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      Kb += _elasticity_tensor[_qp](i, i, j, j);

  //checked, it is good either way I : C : I

  Kb = (1.0 / 9.0) * Kb;
  Je = _fe.det();
  //eta = 1 - Je;
  delta = 1.5 * (std::pow(Je,2.0/3.0) - 1);
  ap_cpl = 0.0;
  ap_vol = 0.0;
  an_vol = 0.0;
  //volumentric free Energy
  if (Je > 1)
    ap_vol = 0.5 * Kb * (std::pow(delta,2.0));
  else
    an_vol = 0.5 * Kb * (std::pow(delta,2.0));

  //volumetric coupling free Energy
  RankTwoTensor elastic_enrgy_tensor;
  elastic_enrgy_tensor = _elasticity_tensor[_qp] * ee;
  elastic_enrgy_tensor = 0.5 * ee * elastic_enrgy_tensor;
  // Doubt

  ap_cpl += elastic_enrgy_tensor.trace();
  ap_cpl = ap_cpl + ap_vol;
  //
  // // Assume isotropic elasticity to calculate elastic strain energy
  // Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1); // Lame parameter
  // Real mu = _elasticity_tensor[_qp](0, 1, 0, 1); // Shear modulus
  //
  // // Compute eigenvectors and eigenvalues of Green-Lagrange strain and projection tensor
  // RankTwoTensor eigvec;
  // std::vector<Real> eigval(LIBMESH_DIM);
  // RankFourTensor Ppos;
  //
  // Ppos = ee.positiveProjectionEigenDecomposition(eigval, eigvec);
  //
  // // Calculate array of tensors of outerproduct of eigen vectors
  // std::vector<RankTwoTensor> etens(LIBMESH_DIM);
  //
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   etens[i].vectorOuterProduct(eigvec.column(i), eigvec.column(i));
  //
  // // Separate out positive and negative eigen values
  // std::vector<Real> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   epos[i] = (std::abs(eigval[i]) + eigval[i]) / 2.0;
  //   eneg[i] = -(std::abs(eigval[i]) - eigval[i]) / 2.0;
  // }
  //
  // // Separate positive and negative sums of all eigenvalues
  // Real etr = 0.0;
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   etr += eigval[i];
  //
  // const Real etrpos = (std::abs(etr) + etr) / 2.0;
  // const Real etrneg = -(std::abs(etr) - etr) / 2.0;
  //
  // Calculate the tensile (positive) and compressive (negative) parts of Cauchy stress
  // RankTwoTensor stress0pos, stress0neg;
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   stress0pos += etens[i] * (lambda * etrpos + 2.0 * mu * epos[i]);
  //   stress0neg += etens[i] * (lambda * etrneg + 2.0 * mu * eneg[i]);
  // }

  // sum squares of epos and eneg
  // Real pval(0.0), nval(0.0);
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   pval += epos[i] * epos[i];
  //   nval += eneg[i] * eneg[i];
  // }

  // Positive part of the stress is degraded by _D[_qp]
  // pk2_new = (_D[_qp] * stress0pos) + stress0neg;
  // if (Je > 1)
    pk2_new = _D[_qp] * _elasticity_tensor[_qp] * ee;
  // else
  //   pk2_new = _D[_qp] *
  // Energy with positive principal strains
  // F_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
  // F_neg = -lambda * etrneg * etrneg / 2.0 + mu * nval;
  F_pos = ap_cpl;
  F_neg = an_vol;

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  //_dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  _dstress_dc[_qp] = _dDdc[_qp] * _elasticity_tensor[_qp] * ee;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history variable
  // if (_use_current_hist)
  //   _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];
  if (_use_current_hist)
     _d2Fdcdstrain[_qp] = _dDdc[_qp] * _elasticity_tensor[_qp] * ee;
  // _Jacobian_mult is already defined in the CP base class

}

// compute history variable and assign to _E
// which is used by the fracture model for damage growth
// Damage grows only because of the positive part of the elastic energy F_pos
void
FiniteStrainUObasedCPDamage_modified::computeHistoryVariable(Real & F_pos, Real & F_neg)
{
  // Assign history variable
  Real hist_variable = _H_old[_qp];

  // _use_snes_vi_solver option not implemented

  if (F_pos > _H_old[_qp])
    _H[_qp] = F_pos;
  else
    _H[_qp] = _H_old[_qp];

  if (_use_current_hist)
    hist_variable = _H[_qp];

  // _barrier not implemented

  // Elastic free energy density and derivatives
  _E[_qp] = hist_variable * _D[_qp] + F_neg;
  _dEdc[_qp] = hist_variable * _dDdc[_qp];
  _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];

}
