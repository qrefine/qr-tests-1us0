import os
from os import *
import time
import datetime
from libtbx import easy_run
import iotbx.pdb

def run():
    """This runner script is used to carry out a large set of SE refinements using different strategies."""
# print "Executing qRefine scripts: " datetime.now()
    total_start_time = time.time()

    # These are the 4 strategies
    #               Quantum REFINE  Classical Refine   Quantum Optimization             Classical Optimization
    strategies = {"refine_cctbx":"restraints=cctbx ","refine_qm":"restraints=qm ","opt_qm":"restraints=qm data_weight=0","opt_cctbx":"restraints=cctbx data_weight=0 "}
    print "Strategies selected are: "

    # This gives us the ability to direct the output PDBs of the script
    results_prefix = "./results_b_rmsd_0.03_terachem_selectrmsd_max32/"#"../p26/04_results/" 
    perturbed_prefix="./perturbed/"#"../p26/03_perturbed/"
    mtz_prefix="./"#"../p26/02_mtz/"
    #p26_pdb_prefix="../p26/01c_fixOXT/"
    #
    folder_name=results_prefix
    make_folder(folder_name)

    #Loop over refinement or optimization strategies
    strategy_start_time = 0

    for strategy_name,strategy in strategies.iteritems():
      if(strategy_name=="refine_cctbx"  ):  
        print strategy
        strategy_start_time = time.time()
#folder_name="../p26/04_results/"+strategy_name
	folder_name=results_prefix+strategy_name
        make_folder(folder_name)
	file_out_prefix = [results_prefix]
#	file_out_prefix = ["../p26/04_results"]
	file_out_prefix.append(strategy_name)
        #Loop over data files
        data_files = os.listdir(mtz_prefix)
        for data_file in data_files:
            data_file_start_time = time.time()
            if(data_file.endswith(".mtz")):
                data_file_name = "."#data_file[:-4] # drop the .mtz extension
                #p26_pdb_file = p26_pdb_prefix+data_file_name+".pdb" 
                #p26_pdb_xray = iotbx.pdb.input(file_name = p26_pdb_file).xray_structure_simple()
		folder_name=results_prefix+strategy_name+"/"+data_file_name
		make_folder(folder_name)	      
                file_out_prefix.append(data_file_name)
		pertubations = os.listdir(perturbed_prefix+data_file_name+"/")
		for pertubation in pertubations:
             #    if(pertubation=="0.3"):
		  if(os.path.isdir(perturbed_prefix+"/"+data_file_name+"/"+pertubation)):
		    file_out_prefix.append(pertubation)
                    make_folder(results_prefix+strategy_name+"/"+data_file_name+"/"+pertubation)
                #Loop over Snapshots
                    snapshotdir_name=perturbed_prefix+data_file_name+"/"+pertubation+"/"
                    snapshots = os.listdir(snapshotdir_name)
                    for snapshot in snapshots:
                        if(snapshot.endswith("pdb")):
                           snapshot_file_name = snapshot[:-4]
#   assert snapshot_file_name == data_file_name # make sure these two files correspond
                           file_out_prefix.append(snapshot_file_name)
#output = "/".join(file_out_prefix)
	                   output=results_prefix+strategy_name+"/"+data_file_name+"/"+pertubation+"/"+snapshot_file_name
			   result_pdb=output+".pdb"
			   result_pdb_exists=os.path.exists(result_pdb) 
                           if( result_pdb_exists != True):
			   #Execute job
                             cmd =  " ".join(["phenix.python ./qr-core/qrefine.py ",snapshotdir_name+snapshot, mtz_prefix+data_file, strategy, "update_all_scales=False"," file_out_prefix="+output, "  >  ", output+ ".log" ])
			     print cmd
                             easy_run.call(cmd)
                           #result_pdb_xray = iotbx.pdb.input(file_name = result_pdb).xray_structure_simple()
                           #rmsd_p26_result = p26_pdb_xray.sites_cart().rms_difference(result_pdb_xray.sites_cart())
                           #print result_pdb, "    rmsd to p26:  " ,rmsd_p26_result  
                           
                # for logging timings, if needed.
                if(0):
                    data_file_time = time.time() - data_file_start_time
                    print "Time taken for " ,data_file_name, "was " , data_file_time
        strategy_time = time.time() - strategy_start_time
        print "Time taken for " ,strategy ,"was ", strategy_time


    total_time = time.time() - total_start_time
    print "Time taken for entire batch was:", total_time

def  make_folder(folder_name):
   if(os.path.exists(folder_name)!=True):
       os.mkdir(folder_name)	

if (__name__ == "__main__"):
  run()
