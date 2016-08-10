from __future__ import division
import iotbx.pdb
import random
import time
import sys
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx import group_args
import mmtbx.utils
import restraints
import macro_cycle
import mmtbx.command_line
import cif_charges
# ase / mopac specific imports
from ase.calculators.mopac import Mopac
from ase.calculators.terachem import TeraChem


if 1:
  random.seed(1)
  flex.set_random_seed(1)

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string="""
sf_algorithm = *direct fft
  .type = choice(multi=False)
restraints = cctbx *qm
  .type = choice(multi=False)
number_of_macro_cycles=50
  .type = int
data_weight=None
  .type = float
file_out_prefix=None
  .type = str
use_convergence_test = True
  .type = bool
update_all_scales = True
  .type = bool
refinement_target_name = *ml ls_wunit_k1
  .type = choice
max_bond_rmsd = 0.03
  .type = float
""", enable_twin_law=False)

def create_fmodel(f_obs, r_free_flags, xray_structure, log, update_all_scales, target_name):
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xray_structure,
    target_name    = target_name)
  if(update_all_scales):
    fmodel.update_all_scales(remove_outliers=False)
    fmodel.show(show_header=False, show_approx=False)
  print >> log, "r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), fmodel.r_free())
  return fmodel

def get_inputs(pdb_file_name):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  return group_args(
    pdb_hierarchy  = pdb_inp.construct_hierarchy(),
    xray_structure = pdb_inp.xray_structure_simple())
def run(args, log = sys.stdout):
  cmdline = mmtbx.command_line.load_model_and_data(
    args          = args,
    master_phil   = get_master_phil(),
    create_fmodel = False,
    out           = log)
  params = cmdline.params
  fmodel = create_fmodel(
    f_obs            = cmdline.f_obs,
    r_free_flags      = cmdline.r_free_flags,
    xray_structure    = cmdline.xray_structure,
    log               = log,
    target_name       = params.refinement_target_name,
    update_all_scales = params.update_all_scales)
  data_weight=None
  if(params.restraints == "cctbx"):
    grm = restraints.get_grm(
      cif_objects  = cmdline.cif_objects,
      file_name    = cmdline.pdb_file_names[0]).restraints_manager
    restraints_manager = restraints.from_cctbx(
      geometry_restraints_manager = grm)
  else:
    assert params.restraints == "qm"
    ### check the total charge and spin
    restraints_manager = restraints.from_qm(
      pdb_hierarchy = cmdline.pdb_hierarchy.deep_copy(),
  #    calculator    = Mopac(label='./ase/ase',charge="0", spin=0,
  #                          job_type='1SCF GRADIENTS AUX(0,PRECISION=9)',
  #                          RELSCF=None,functional='PM7'))
      calculator     = TeraChem(label='./ase/ase')) 
  states = mmtbx.utils.states(
    pdb_hierarchy  = cmdline.pdb_hierarchy,
    xray_structure = cmdline.xray_structure)
  states.add(sites_cart = cmdline.xray_structure.sites_cart())
  print "params.data_weight:", params.data_weight
  print "params.restraints:", params.restraints
  cif_file=None
  if len(cmdline.cif_objects)!=0:
    cif_file=cmdline.cif_objects[0][0]		
  res_charges=cif_charges.res_charges(cmdline.pdb_hierarchy.deep_copy(),cif_file)	
  print "residue charges",res_charges
  #####
  geometry_rmsd_manager=restraints.get_grm(
    cif_objects  = cmdline.cif_objects,
    file_name    = cmdline.pdb_file_names[0]).restraints_manager
  macro_cycle.run(
    fmodel                 = fmodel,
    states                 = states,
    pdb_hierarchy          = cmdline.pdb_hierarchy,
    restraints_manager     = restraints_manager,
    number_of_macro_cycles = params.number_of_macro_cycles,
    file_label             = params.restraints,
    data_weight            = params.data_weight,
    file_out_prefix        = params.file_out_prefix,
    use_convergence_test   = params.use_convergence_test,
    max_bond_rmsd          = params.max_bond_rmsd,
    geometry_rmsd_manager  = geometry_rmsd_manager)


if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time: %6.4f"%(time.time()-t0)
