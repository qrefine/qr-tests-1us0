from __future__ import division
import iotbx.pdb
import minimizer
import os
from scitbx.array_family import flex
import copy

class convergence(object):
  def __init__(self, fmodel, r_tolerance, rmsd_tolerance):
    self.r_start = fmodel.r_work()
    self.sites_cart_start = fmodel.xray_structure.sites_cart()
    self.r_tolerance=r_tolerance
    self.rmsd_tolerance=rmsd_tolerance
    self.number_of_convergence_occurances=0

  def is_converged(self, fmodel):
    r = fmodel.r_work()
    sites_cart = fmodel.xray_structure.sites_cart()
    r_diff=abs(self.r_start-r)
    rmsd_diff=self.sites_cart_start.rms_difference(sites_cart)
    self.sites_cart_start = sites_cart
    self.r_start=r
    if(r_diff<self.r_tolerance and rmsd_diff<self.rmsd_tolerance):
      self.number_of_convergence_occurances+=2
    if(self.number_of_convergence_occurances==2 or r<0.005):
      return True
    else: return False

class results(object):
  def __init__(self, r, b, xrs, max_bond_rmsd):
    self.max_bond_rmsd = max_bond_rmsd
    self.rs   = flex.double()
    self.bs   = flex.double()
    self.xrss = []
    self.rs.append(r)
    self.bs.append(b)
    self.xrss.append(xrs.deep_copy_scatterers())

  def update(self, r, b, xrs):
    self.rs.append(r)
    self.bs.append(b)
    self.xrss.append(xrs.deep_copy_scatterers())

  def choose_best(self):
    # This modifies in-place the content of class member, not a good thing to do
    # but fine as a quick fix.
    self.rs   = self.rs  [1:]
    self.bs   = self.bs  [1:]
    self.xrss = self.xrss[1:]
    sel = self.bs < self.max_bond_rmsd
    if(sel.size()==0): return None
    min_r = flex.min(self.rs.select(sel))
    i_best = None
    for i in xrange(self.rs.size()):
      if(abs(self.rs[i]-min_r)<1.e-5):
        i_best = i
        break
    return self.xrss[i_best].deep_copy_scatterers(), self.rs[i_best]

    # XXX old, not using free-r
    #sel_where = self.rs < 0.2
    #sel = self.rs.select(sel_where)
    #if(sel==None or sel.size()==0): return None
    #bs = self.bs.select(sel_where)
    #xrss = [self.xrss[i].deep_copy_scatterers() for i in sel_where.iselection()]
    #sel = flex.sort_permutation(bs)
    #bs_all_list = list(self.bs)
    #bs_list = list(bs)
    #print "best mc:   ",bs_all_list.index(min(bs_list))-1
    #return xrss[sel[0]]

def get_initial_restraints_weight_scale(
      initial_fmodel,
      restraints_manager,
      sites,
      data_weight,
      max_bond_rmsd,
      geometry_rmsd_manager):
  initial_restraints_weight_scale = 32
  counter=0
  while(counter<10):
    counter+=1
    fmodel = initial_fmodel.deep_copy()
    minimizer.run(
      fmodel                  = fmodel,
      restraints_manager      = restraints_manager,
      sites                   = True,
      data_weight             = data_weight,
      restraints_weight_scale = initial_restraints_weight_scale)
    cctbx_rm_bonds_rmsd =  geometry_rmsd_manager.energies_sites(
      sites_cart        = fmodel.xray_structure.sites_cart(),
      compute_gradients = False).geometry.bond_deviations()[2]
    print "Initial weight optimization counter:", \
      counter, initial_restraints_weight_scale, cctbx_rm_bonds_rmsd
    if(cctbx_rm_bonds_rmsd<=max_bond_rmsd):
      initial_restraints_weight_scale /= 2.
    else:
      return initial_restraints_weight_scale*2.0
  return 1

def run(fmodel,
        restraints_manager,
        states,
        pdb_hierarchy,
        number_of_macro_cycles,
        max_bond_rmsd,
        file_label,
        data_weight           = None,
        file_out_prefix       = None,
        r_tolerance           = 0.001,
        rmsd_tolerance        = 0.01,
        use_convergence_test  = False,
        geometry_rmsd_manager = None):
  cctbx_rm_bonds_rmsd =  geometry_rmsd_manager.energies_sites(
    sites_cart        = fmodel.xray_structure.sites_cart(),
    compute_gradients = False).geometry.bond_deviations()[2]
  print "rmsd(b), start: %7.4f"%cctbx_rm_bonds_rmsd
  conv_test = convergence(
    fmodel=fmodel, r_tolerance=r_tolerance, rmsd_tolerance=rmsd_tolerance)
  restraints_weight_scale = 1.0
  if(data_weight is None or abs(data_weight)>1.e-9):
    # looks like we don't need it?
    #
    #restraints_weight_scale = get_initial_restraints_weight_scale(
    #  initial_fmodel          = fmodel,
    #  restraints_manager      = restraints_manager,
    #  sites                   = True,
    #  data_weight             = data_weight,
    #  max_bond_rmsd           = max_bond_rmsd,
    #  geometry_rmsd_manager   = geometry_rmsd_manager)
    print "initial_restraints_weight_scale", restraints_weight_scale
    res = results(r=fmodel.r_free(), b=cctbx_rm_bonds_rmsd,
      xrs=fmodel.xray_structure, max_bond_rmsd=max_bond_rmsd)
  #
  fmodel_copy = fmodel.deep_copy()
  rws = [restraints_weight_scale]
  for macro_cycle in xrange(number_of_macro_cycles):
    # refine coordinates
    fmodel = fmodel_copy.deep_copy()
    conv_test = convergence(
      fmodel=fmodel, r_tolerance=r_tolerance, rmsd_tolerance=rmsd_tolerance)
    for i in xrange(15):
      minimized = minimizer.run(
        fmodel                  = fmodel,
        restraints_manager      = restraints_manager,
        sites                   = True,
        data_weight             = data_weight,
        restraints_weight_scale = restraints_weight_scale)
      fmodel.update_all_scales(remove_outliers=False)
      if(use_convergence_test and conv_test.is_converged(fmodel=fmodel)):
        print "  Convergence reached at i:", i
        break
      cctbx_rm_bonds_rmsd =  geometry_rmsd_manager.energies_sites(
        sites_cart        = fmodel.xray_structure.sites_cart(),
        compute_gradients = False).geometry.bond_deviations()[2]
      if(cctbx_rm_bonds_rmsd>max_bond_rmsd*10.0):
        print "  Distorted at i:", i
        break
    states.add(sites_cart=fmodel.xray_structure.sites_cart())
    rw = fmodel.r_work()
    rf = fmodel.r_free()
    #
    cctbx_rm_bonds_rmsd =  geometry_rmsd_manager.energies_sites(
      sites_cart        = fmodel.xray_structure.sites_cart(),
      compute_gradients = False).geometry.bond_deviations()[2]
    if(data_weight is None or abs(data_weight)>1.e-9):
      if(cctbx_rm_bonds_rmsd<0.01):
        restraints_weight_scale /= 2.
      if(cctbx_rm_bonds_rmsd>max_bond_rmsd or (rf>rw and abs(rf-rw)*100.>5.)):
        restraints_weight_scale *= 2.
    ###
    if(data_weight is None or abs(data_weight)>1.e-9):
      res.update(r=fmodel.r_free(), b=cctbx_rm_bonds_rmsd, xrs=fmodel.xray_structure)
    fmt="mc: %3d Rw: %6.4f Rf: %6.4f rmsd(b): %7.4f restraint_scale: %4.1f min_steps: %d"
    print fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free(), cctbx_rm_bonds_rmsd,
      restraints_weight_scale, minimized.minimizer.nfun())
    if(restraints_weight_scale in rws):
      print " Convergence at macro_cycle:", macro_cycle
      break
    rws.append(restraints_weight_scale)
  ###
  if(data_weight is None or abs(data_weight)>1.e-9):
    xrs_best, r_best = res.choose_best()
    if(xrs_best is not None):
      fmodel.update_xray_structure(xray_structure = xrs_best,
        update_f_calc=True, update_f_mask=True)
      fmodel.update_all_scales(remove_outliers=False)
      pdb_hierarchy.adopt_xray_structure(xrs_best)
      print "Best r_work: %6.4f r_free: %6.4f"%(fmodel.r_work(),fmodel.r_free())
    else:
      print " r_factor (best): None"
      print " take the last structure"
      print "Best r_work: %6.4f r_free: %6.4f"%(fmodel.r_work(),fmodel.r_free())
  ###
  if(1):
    pdb_hierarchy.adopt_xray_structure(fmodel.xray_structure)
    if(file_out_prefix==None):
      file_out_prefix="./pdb/"
      if( os.path.exists(file_out_prefix)) !=True:
        os.mkdir(file_out_prefix)
      job_type= "_refined" if (data_weight==None)  else "_opt"
      file_out_prefix=file_out_prefix+file_label+ job_type
    pdb_hierarchy.write_pdb_file(file_name=file_out_prefix+".pdb",
      crystal_symmetry = fmodel.xray_structure.crystal_symmetry())
    states.write(file_name=file_out_prefix+"_all_states.pdb")
    print "see the result structure: ", file_out_prefix+".pdb"
