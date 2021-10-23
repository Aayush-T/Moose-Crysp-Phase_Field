[GlobalParams]
  displacements = 'disp_x disp_y'
#  outputs = 'csv exodus'
[]

# add damage variable
#[Variables]
#  [./c]
#    order = FIRST
#    family = LAGRANGE
#	[./InitialCondition]
#      type = ConstantIC
#      value = '0.5'
#    [../]
#  [../]
#[]

[Mesh]
  type = FileMesh
  file = crack_mesh.e
[]

[AuxVariables]
  [./pk2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./slip_increment]
   order = CONSTANT
   family = MONOMIAL
  [../]
#  [./el_str]
#   order = CONSTANT
#   family = MONOMIAL
#  [../]
#  [./pl_str]
#   order = CONSTANT
#   family = MONOMIAL
#  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./All]
        strain = FINITE
        add_variables = true
        generate_output = stress_yy
        #use_displaced_mesh = true
      [../]
    [../]
  [../]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = F
        kappa = kappa_op
        mobility = L
      [../]
    [../]
  [../]
[]

# damage is constant in this test
# dc/dt = 0
[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
    use_displaced_mesh = true
  [../]
  [./off_disp]
    type = AllenCahnElasticEnergyOffDiag
    variable = c
    displacements = 'disp_x disp_y'
    mob_name = L
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./pk2]
   type = RankTwoAux
   variable = pk2
   rank_two_tensor = pk2
   index_j = 1
   index_i = 1
   execute_on = timestep_end
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./slip_inc]
   type = MaterialStdVectorAux
   variable = slip_increment
   property = slip_rate_gss
   index = 0
   execute_on = timestep_end
  [../]
  [./gss]
    type = MaterialStdVectorAux
    variable = gss
    property = state_var_gss
    index = 0
    execute_on = timestep_end
  [../]
#  [./el_str]
#    type = RankTwoAux
#    variable = el_str
#    rank_two_tensor = ee
#    index_j = 1
#    index_i = 1
#    execute_on = timestep_end
#  [../]
#  [./pl_str]
#    type = RankTwoAux
#    variable = pl_str
#    rank_two_tensor = ee
#    index_j = 1
#    index_i = 1
#    execute_on = timestep_end
#  [../]
[]

[BCs]
[./ydisp]
  type = FunctionDirichletBC
  variable = disp_y
  boundary = 2
  #function = '0.0001*t'
  function = '0.01*t'
[../]
[./yfix]
  type = DirichletBC
  variable = disp_y
  boundary = 1
  value = 0
[../]
[./xfix]
  type = DirichletBC
  variable = disp_x
  boundary = '1 2'
  value = 0
[../]
[]

[UserObjects]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSS
    variable_size = 12
    slip_sys_file_name = input_slip_sys.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1'
    uo_state_var_name = state_var_gss
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 4 8 12'
    group_values = '60.8 60.8 60.8'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1.0 541.5 109.8 2.5'
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '1e-3 0.05 1e-6'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  [./crysp] # new class with damage
    #type = FiniteStrainUObasedCP_PF_new
    #type = FiniteStrainUObasedCP_PF_new_pk2new_global
    #type = FiniteStrainUObasedCP_PF_new_pk2new_global_strain
    type = FiniteStrainUObasedCPDamage_modified
	  c = c
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    E_name = 'elastic_energy'
    D_name = 'degradation'
  #  F_name = 'local_fracture_energy'
  #  decomposition_type = stress_spectral
    use_current_history_variable = true
    outputs = exodus
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    #C_ijkl = '120.0 80.0'
    #fill_method = symmetric_isotropic
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1.0e-6'
    derivative_order = 2
    outputs = exodus
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
    outputs = exodus
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
    outputs = exodus
  [../]
[]

[Postprocessors]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./pk2]
   type = ElementAverageValue
   variable = pk2
  [../]
  [./fp_yy]
    type = ElementAverageValue
    variable = fp_yy
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
  [../]
  [./gss]
    type = ElementAverageValue
    variable = gss
  [../]
  [./slip_increment]
   type = ElementAverageValue
   variable = slip_increment
  [../]
  #[./el_str]
  #  type = ElementAverageValue
  #  variable = el_str
  #[../]
  #[./pl_str]
  #  type = ElementAverageValue
  #  variable = pl_str
  #[../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.001
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  nl_abs_tol = 1e-10
  nl_rel_step_tol = 1e-10
  dtmax = 10.0
  nl_rel_tol = 1e-10
  dtmin = 0.001
  num_steps = 2000
  nl_abs_step_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
