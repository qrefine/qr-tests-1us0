from __future__ import division
from scitbx.array_family import flex
from cctbx import xray
import scitbx.lbfgs
import random

def show_histogram(data, n_slots):
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  tmp = None
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    #print "%10.3f - %-10.3f : %d" % (lc_1, hc_1, n_1)
    lc_1 = hc_1
    if(tmp is None): tmp = hc_1
  return tmp

def compute_weight(fmodel, rm):
  random.seed(1)
  flex.set_random_seed(1)
  #
  fmodel_dc = fmodel.deep_copy()
  xrs = fmodel_dc.xray_structure.deep_copy_scatterers()
  xrs.shake_sites_in_place(mean_distance=0.2)
  fmodel_dc.update_xray_structure(xray_structure=xrs, update_f_calc=True)
  x_target_functor = fmodel_dc.target_functor()
  tgx = x_target_functor(compute_gradients=True)
  gx = flex.vec3_double(tgx.\
          gradients_wrt_atomic_parameters(site=True).packed())
  #
  tc, gc = rm.target_and_gradients(sites_cart=xrs.sites_cart())
  x = gc.norm()
  y = gx.norm()
#  ################
#  #print "-"*10
#  #print list(flex.sqrt(gx.dot()))
#  #print list(flex.sqrt(gc.dot()))
#  gx_d = flex.sqrt(gx.dot())
#  vx = show_histogram(data=gx_d, n_slots=10)
#  sel = gx_d<vx*10
#  #y = gx.select(sel).norm()
#  #print
#  gc_d = flex.sqrt(gc.dot())
#  vc = show_histogram(data=gc_d, n_slots=10)
#  sel = gc_d<vc*10
#  x = gc.select(sel).norm()
#  #print sel.size(), sel.count(True)
#  #print "-"*10
#  #x = gc.select(sel).norm()
#  #y = gx.select(sel).norm()
#  ################

  if(y != 0.0): return x/y
  else:         return 1.0 # ad hoc default fallback

class run(object):
  def __init__(self,
        fmodel,
        restraints_manager,
        max_iterations=50,
        sites = False,
        u_iso = False,
        restraints_weight=1,
        data_weight=None,
        restraints_weight_scale=1):
    self.cntr = 0
    self.rws = restraints_weight_scale
    self.fmodel = fmodel
    self.restraints_weight = restraints_weight
    self.data_weight = data_weight
    self.restraints_manager = restraints_manager
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.sites = sites
    self.u_iso = u_iso
    if(self.sites):
      self.x = self.fmodel.xray_structure.sites_cart().as_double()
    if(self.u_iso):
      assert self.fmodel.xray_structure.scatterers().size() == \
        self.fmodel.xray_structure.use_u_iso().count(True)
      self.x = self.fmodel.xray_structure.extract_u_iso_or_u_equiv()
    if(self.sites):
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)
    if(self.u_iso):
      sel = flex.bool(
        self.fmodel.xray_structure.scatterers().size(), True).iselection()
      self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
        iselection = sel)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    self.fmodel.xray_structure.tidy_us()
    self.fmodel.xray_structure.apply_symmetry_sites()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True,
      update_f_mask  = True)

  def compute_functional_and_gradients(self):
    sites_cart = flex.vec3_double(self.x)
    if(self.sites):
      self.fmodel.xray_structure.set_sites_cart(
        sites_cart = sites_cart)
    if(self.u_iso):
      self.fmodel.xray_structure.set_u_iso(values = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    tgx = self.x_target_functor(compute_gradients=True)
    if(self.sites):
      self.cntr+=1
      tx,gx = None,None
      if(self.data_weight is None or abs(self.data_weight)>1.e-9):
        # data (x-ray) target and grads
        tx = tgx.target_work()
        gx = flex.vec3_double(tgx.\
          gradients_wrt_atomic_parameters(site=True).packed())
      # restraints target and grads
      tc, gc = self.restraints_manager.target_and_gradients(
        sites_cart=sites_cart)
      # compute weight
      if(self.cntr == 1 and self.data_weight is None):
        self.data_weight = compute_weight(
          fmodel=self.fmodel,
          rm     = self.restraints_manager)
      # compute target and grads
      if(abs(self.data_weight)>1.e-9):
        ## refinement
        # add them together to get final total target
        f = tx*self.data_weight + self.restraints_weight*tc*self.rws
        g = gx*self.data_weight + self.restraints_weight*gc*self.rws
      else:
        ##optimization
        f = self.restraints_weight*tc
        g = self.restraints_weight*gc
    if(self.u_iso):
      tx = tgx.target_work()
      gx = tgx.gradients_wrt_atomic_parameters(u_iso=True)
      f = tx
      g = gx
    return f, g.as_double()
