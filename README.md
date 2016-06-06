

1) Compute simulated diffraction data:
  
  phenix.fmodel a87_99_h.pdb add_random_error_to_amplitudes_percent=5 type=real high_res=4 low_res=6 r_free=0.1
  
2) Run Q|R refinement, sample command:
  
  python ../qr-core/qrefine.py a87_99_h.pdb.mtz perturbed/1.5/4.pdb restraints=cctbx update_all_scales=False
