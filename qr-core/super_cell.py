import iotbx.pdb

class expand(object):
  def __init__(self, pdb_hierarchy, crystal_symmetry):
    self.pdb_hierarchy = pdb_hierarchy
    self.crystal_symmetry = crystal_symmetry
    self.ph_p1 = self.pdb_hierarchy.expand_to_p1(
      crystal_symmetry=self.crystal_symmetry)
    # super-cell
    self.root = iotbx.pdb.hierarchy.root()
    sites_cart = self.ph_p1.atoms().extract_xyz()
    a,b,c = self.crystal_symmetry.unit_cell().parameters()[:3]
    cntr=0
    for x in [-a, 0, a]:
      for y in [-b, 0, b]:
        for z in [-c, 0, c]:
          ph_ = self.ph_p1.deep_copy()
          ph_.atoms().set_xyz(sites_cart+[x,y,z])
          models = ph_.models()
          md = models[0].detached_copy()
          md.id = str(cntr)
          self.root.append_model(md)
          cntr+=1

  def write_super_cell(self, file_name="super.pdb"):
    self.root.write_pdb_file(file_name=file_name)

  def write_p1(self, file_name="p1.pdb"):
    self.ph_p1.write_pdb_file(file_name=file_name,
      crystal_symmetry=self.crystal_symmetry)

def exercise():
  pdb_inp = iotbx.pdb.input(file_name="1yjp.pdb")
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()

  o = expand(pdb_hierarchy = ph, crystal_symmetry = cs)
  o.write_super_cell()
  o.write_p1()

if __name__ == "__main__":
  exercise()
