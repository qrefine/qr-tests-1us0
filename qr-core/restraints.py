from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex
import ase.units as ase_units
from ase.io import read as ase_io_read
import os
import mmtbx.restraints
from libtbx import group_args
import mmtbx.utils
from libtbx.utils import null_out

def get_grm(file_name, cif_objects):
  from mmtbx import monomer_library
  import mmtbx.monomer_library.server
  import mmtbx.monomer_library.pdb_interpretation
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
    pdb_interpretation_params = params,
    stop_for_unknowns         = True,
    log                       = null_out(),
    cif_objects               = cif_objects,
    use_neutron_distances     = True)
  processed_pdb_file, junk = processed_pdb_files_srv.\
    process_pdb_files(pdb_file_names = [file_name])
  xray_structure = processed_pdb_file.xray_structure()
  sctr_keys = \
    xray_structure.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    assume_hydrogens_all_missing = not has_hd,
    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(
    geometry = geometry, normalization = False)
  return group_args(
    restraints_manager = restraints_manager,
    pdb_hierarchy      = processed_pdb_file.all_chain_proxies.pdb_hierarchy)

class from_cctbx(object):
  def __init__(self, geometry_restraints_manager):
    self.geometry_restraints_manager = geometry_restraints_manager

  def target_and_gradients(self, sites_cart):
    es = self.geometry_restraints_manager.energies_sites(
      sites_cart=sites_cart, compute_gradients=True)
    return es.target, es.gradients

class from_qm(object):
  def __init__(self, pdb_hierarchy, file_name = "./ase/tmp_ase.pdb",calculator=None):
    folder_name="./ase"
    if(os.path.exists(folder_name)!=True):
      os.mkdir(folder_name)
    self.pdb_hierarchy = pdb_hierarchy
    self.file_name = file_name
    self.calculator=calculator

  def target_and_gradients(self, sites_cart):
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    self.pdb_hierarchy.write_pdb_file(file_name=self.file_name)
    atoms=ase_io_read(self.file_name)
    atoms.set_calculator(self.calculator)
    atoms.get_calculator().run()
    calculator=atoms.get_calculator()
    energy= calculator.energy_free
    gradients=(-1.0)*calculator.forces
    unit_convert=ase_units.mol/ase_units.kcal
    return energy*unit_convert, flex.vec3_double(gradients*unit_convert)
