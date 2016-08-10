from __future__ import division
import iotbx.pdb
import libtbx.load_env
import os

def res_charges(pdb_hierarchy, cif_file=None):
  aa=['GLU', 'ASP','ARG', 'LYS','HIS', 'SER','THR','ASN','GLN','CYS','SEC','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']	
  charge_dict = dict.fromkeys(['GLU', 'ASP'], -1)
  charge_dict.update(dict.fromkeys(['ARG', 'LYS'], 1))
  charge_dict.update(dict.fromkeys(['HIS', 'SER','THR','ASN','GLN','CYS','SEC','GLY','PRO'], 0))
  charge_dict.update(dict.fromkeys(['ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP'],0))
  # read charge_dict from given cif file	
  cif_file_charge(cif_file,charge_dict)	
  charges=[]
  for chain in pdb_hierarchy.chains():
    chain_residue_groups=list(chain.residue_groups())
    size_chain_residue_groups=len(chain_residue_groups)  
    for i in range(size_chain_residue_groups):
      residue_group=chain_residue_groups[i]  	
      resname=list(residue_group.unique_resnames())[0]
      if resname not in charge_dict:
	# read charge_dict from cif file in cctbx library
        charge_dict.update(dict.fromkeys([resname], 0.0))
        file_name="chem_data/geostd/"+resname[0].lower()+"/data_"+resname+".cif"
        file_path = libtbx.env.find_in_repositories(
    		relative_path=file_name,
    		test=os.path.isfile)
        if(file_path==None):
          file_name="chem_data/mon_lib"+resname[0].lower()+"/"+resname+".cif"
          file_path = libtbx.env.find_in_repositories(
                  relative_path=file_name,
                  test=os.path.isfile)
          if(file_path==None):
	    print resname ,"please put its parameters in the ligand cif file, not found in library"	
        cif_file_charge(file_path,charge_dict)
      atoms=residue_group.atoms()
      atom_names=[]
      for atom in atoms:
	atom_names.append(atom.name.strip())
      # charge: NTerminal+1, CTerminal-1	
      if resname in aa and "H1" in atom_names:
	charges.append(charge_dict[resname]+1.0)
      elif resname in aa and "OXT" in atom_names: 
	charges.append(charge_dict[resname]-1.0)
      else:    
        charges.append(charge_dict[resname])	
#      if i==0:
#        charges.append(charge_dict[resname]+1.0)
#      elif i==size_chain_residue_groups-1:
#        charges.append(charge_dict[resname]-1.0) 
#      else:    
#        charges.append(charge_dict[resname])
  return charges

def cif_file_charge(cif_file=None,charge_dict=None):
  if cif_file != None and cif_file.endswith(".cif"):
    print "cif file :",cif_file
    charge_file = open(cif_file,"r").readlines()
    charges={}
    for line in charge_file:
      if 'data_comp_' in line and 'data_comp_list' not in line:	
	resname=line.split()[0][10:]
	if resname not in charges.keys():	
           	charges[resname]=0.0
    for line in charge_file:
      columns=line.split()	
      if(len(columns)==9):
        letter_in_column4=False
	for  c in  columns[4]:
	  if c.isalpha():
	     letter_in_column4=True
	     break		
	if columns[0] in charges.keys() and ( not letter_in_column4):
          charges[columns[0]]=charges[columns[0]]+float(columns[4])
    for resname in charges.keys():		 
      charge_dict.update(dict.fromkeys([resname], charges[resname]))  


if __name__ == "__main__":
  run()
