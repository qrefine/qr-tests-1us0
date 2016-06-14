from __future__ import division
import math
from libtbx import easy_run
import os
import iotbx.pdb
from mmtbx import model_statistics
from scitbx.array_family import flex
from libtbx import group_args
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation

def get_grm(file_name):
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl=False
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = monomer_library.server.server(),
    ener_lib       = monomer_library.server.ener_lib(
      use_neutron_distances=True),
    file_name      = file_name,
    params         = params,
    force_symmetry = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    assume_hydrogens_all_missing = False,
    plain_pairs_radius = 5.0)
  return group_args(
    geometry      = geometry,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)

def get_model_stat(file_name):
  grm = get_grm(file_name = file_name)
  es = grm.geometry.energies_sites(
    sites_cart=grm.pdb_hierarchy.atoms().extract_xyz())
  return es.bond_deviations()[2]

def dist(site1, site2):
  return math.sqrt(
    (site1[0]-site2[0])**2 +
    (site1[1]-site2[1])**2 +
    (site1[2]-site2[2])**2)

def run(file_name="a87_99_h.pdb"):
  #
  prefix="results_terachem"
  #sub_path = "%s/refine_cctbx/"%prefix
  rmsd_dirs = ["0.3/","0.6/","0.9/","1.2/","1.5/"]
  file_names = ["0.pdb","1.pdb","2.pdb","3.pdb","4.pdb","5.pdb","6.pdb","7.pdb",
                "8.pdb","9.pdb"]
  # Identify i_seqs of O-H pairs involved into H bonds
  h_bonds_i_seqs = []
  h = iotbx.pdb.input(file_name="a87_99_h.pdb").construct_hierarchy()
  atoms = h.atoms()
  for i, a_i in enumerate(list(atoms)):
    for j, a_j in enumerate(list(atoms)):
      if(i<j):
        if(a_i.name.strip() in ["H","O"] and a_j.name.strip() in ["O","H"] and
          a_i.name.strip() != a_j.name.strip()):
          d = dist(a_i.xyz, a_j.xyz)
          if(d<2.2 and d>1.7):
            #print a_i.name.strip(), a_j.name.strip(), a_j.i_seq, a_i.i_seq, \
            #  a_j.i_seq-a_i.i_seq, d
            h_bonds_i_seqs.append([a_i.i_seq, a_j.i_seq])
  print "Total number of H bonds in one file:", len(h_bonds_i_seqs)
  # loop over refinements
  for sub_path in [#"opt_cctbx/",
                   #"opt_qm/",
                   "refine_cctbx/",
                   "refine_qm/"
                  ]:
    print sub_path
    for rmsd_dir in rmsd_dirs:
      h_bonds = flex.double()
      cntr = 0
      for fn in file_names:
        file_name = prefix+"/"+sub_path+rmsd_dir+fn
        if(os.path.exists(file_name)==False): continue
        #file_name = "/Users/pafonine/qr/qrefine/helix/a87_99_h.pdb"
        bond_rmsd_mean = get_model_stat(file_name=file_name)
        if(bond_rmsd_mean>0.05): continue
        sites_cart = iotbx.pdb.input(file_name=file_name).atoms().extract_xyz()
        for pair in h_bonds_i_seqs:
          d = dist(sites_cart[pair[0]], sites_cart[pair[1]])
          h_bonds.append(d)
        cntr+=1
      sel  = h_bonds < 2.2
      sel &= h_bonds > 1.7
      if(h_bonds.size()>0):
        print rmsd_dir, "%8.3f %8.3f %8.3f"%h_bonds.min_max_mean().as_tuple(), \
          "%7.2f"%(sel.count(True)*100./(len(h_bonds_i_seqs)*cntr))
      else:
        print rmsd_dir, "N/A"

if __name__ == "__main__":
  run()
