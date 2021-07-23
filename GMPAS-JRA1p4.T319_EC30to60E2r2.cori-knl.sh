#####################################################################
# Template for E3SM-ARRM simulations: G-cases.
# Written for NERSC, directory modifications needed for other LCFs.
# Author: Darin Comeau, LANL.
#
# Things to check:
# - change USER
# - change REPO_ROOT, CASE_ROOT directories if desired
# - check COMPSET
# - check GRID
# - check PE_LAYOUT
# - check RUN_OPTIONS, default is 5 day test
# 
#####################################################################

#####################################################################
# Define user paths
#####################################################################

export USER=ameza9
export PROJECT=e3sm
export MACHINE=cori-knl
export REPO_ROOT=/global/cfs/cdirs/$PROJECT/$USER/e3sm_repos
export CASE_ROOT=/global/cfs/cdirs/$PROJECT/$USER/e3sm_cases
#where the model output will be set cannot and should not be changes
export RUN_ROOT=/global/cscratch1/sd/$USER/e3sm_scratch/$MACHINE

if [ ! -d "$REPO_ROOT" ]; then
	mkdir $REPO_ROOT
fi
if [ ! -d "$CASE_ROOT" ]; then
	mkdir $CASE_ROOT
fi

#####################################################################
# Download repository
#####################################################################

#get project output from master repository 
#export REPO_NAME=master.`date +"%Y%m%d"`
export REPO_NAME="E3SM_playground"
export E3SM_ROOT=$REPO_ROOT/$REPO_NAME

if [ ! -d "$REPO_ROOT/$REPO_NAME" ]; then
	mkdir $REPO_ROOT/$REPO_NAME
	cd $REPO_ROOT/$REPO_NAME
	git clone git@github.com:E3SM-Project/E3SM.git
	cd E3SM
	git submodule update --init --recursive
fi

#####################################################################
# Define case: Different simulations are called cases depending on configuration
#####################################################################
# G-gcase(active ocean and sea ice, prescribed atmosphere and land)
# JRA - forcing data set (newer than IAS) 
# GRID - TL319 (where atmosphere lives) 
#     - EC30to60 30km resolution at high latitudes, 60km at mid-latitudes 
# TAG - case name for simulation 


#export COMPSET=GMPAS-IAF
#export PRECISION=single
export COMPSET=GMPAS-JRA1p4
export GRID=TL319_EC30to60E2r2
export COMPILER=intel
export TAG=3monthtest_single
export E3SM_CASE=`date +"%Y%m%d"`.${COMPSET}.${GRID}.${TAG}.${MACHINE}

#####################################################################
# Create case (won't have to change) 
#####################################################################

cd $E3SM_ROOT/cime/scripts
./create_newcase \
-case $CASE_ROOT/$E3SM_CASE \
-compiler $COMPILER \
-mach $MACHINE \
-project $PROJECT \
-compset $COMPSET \
-res $GRID

cd $CASE_ROOT/$E3SM_CASE

#####################################################################
# PE_LAYOUT (Processor Layout, might not need to change after set)
#####################################################################

# Large layout, Cori-KNL
# Number of processors being used 
# Each node has 68 processors 

export NPROCS_CPL=640
export NPROCS_ATM=640
export NPROCS_LND=640
export NPROCS_ROF=640
export NPROCS_ICE=640
export NPROCS_OCN=960
export NPROCS_GLC=640
export NPROCS_WAV=640

export NTHRDS_CPL=1
export NTHRDS_ATM=1
export NTHRDS_LND=1
export NTHRDS_ROF=1
export NTHRDS_ICE=1
export NTHRDS_OCN=1
export NTHRDS_GLC=1
export NTHRDS_WAV=1

export ROOTPE_CPL=0
export ROOTPE_ATM=0
export ROOTPE_LND=0
export ROOTPE_ROF=0
export ROOTPE_ICE=0
export ROOTPE_OCN=680
export ROOTPE_GLC=0
export ROOTPE_WAV=0

#writes into an xml file, provides a record for all changes done in case directory 
./xmlchange -file env_mach_pes.xml  -id NTASKS_CPL  -val $NPROCS_CPL
./xmlchange -file env_mach_pes.xml  -id NTASKS_ATM  -val $NPROCS_ATM
./xmlchange -file env_mach_pes.xml  -id NTASKS_LND  -val $NPROCS_LND
./xmlchange -file env_mach_pes.xml  -id NTASKS_ROF  -val $NPROCS_ROF
./xmlchange -file env_mach_pes.xml  -id NTASKS_ICE  -val $NPROCS_ICE
./xmlchange -file env_mach_pes.xml  -id NTASKS_OCN  -val $NPROCS_OCN
./xmlchange -file env_mach_pes.xml  -id NTASKS_GLC  -val $NPROCS_GLC
./xmlchange -file env_mach_pes.xml  -id NTASKS_WAV  -val $NPROCS_WAV

./xmlchange -file env_mach_pes.xml  -id NTHRDS_CPL  -val $NTHRDS_CPL
./xmlchange -file env_mach_pes.xml  -id NTHRDS_ATM  -val $NTHRDS_ATM
./xmlchange -file env_mach_pes.xml  -id NTHRDS_LND  -val $NTHRDS_LND
./xmlchange -file env_mach_pes.xml  -id NTHRDS_ROF  -val $NTHRDS_ROF
./xmlchange -file env_mach_pes.xml  -id NTHRDS_ICE  -val $NTHRDS_ICE
./xmlchange -file env_mach_pes.xml  -id NTHRDS_OCN  -val $NTHRDS_OCN
./xmlchange -file env_mach_pes.xml  -id NTHRDS_GLC  -val $NTHRDS_GLC
./xmlchange -file env_mach_pes.xml  -id NTHRDS_WAV  -val $NTHRDS_WAV

./xmlchange -file env_mach_pes.xml  -id ROOTPE_CPL  -val $ROOTPE_CPL
./xmlchange -file env_mach_pes.xml  -id ROOTPE_ATM  -val $ROOTPE_ATM
./xmlchange -file env_mach_pes.xml  -id ROOTPE_LND  -val $ROOTPE_LND
./xmlchange -file env_mach_pes.xml  -id ROOTPE_ROF  -val $ROOTPE_ROF
./xmlchange -file env_mach_pes.xml  -id ROOTPE_ICE  -val $ROOTPE_ICE
./xmlchange -file env_mach_pes.xml  -id ROOTPE_OCN  -val $ROOTPE_OCN
./xmlchange -file env_mach_pes.xml  -id ROOTPE_GLC  -val $ROOTPE_GLC
./xmlchange -file env_mach_pes.xml  -id ROOTPE_WAV  -val $ROOTPE_WAV

./xmlchange -file env_mach_pes.xml -id MAX_MPITASKS_PER_NODE -val 64
./xmlchange -file env_mach_pes.xml -id MAX_TASK:S_PER_NODE -val 64

#####################################################################
# Setup and build
#####################################################################

#./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS --append --val='-cosp'
# Sets up and builds the model (typically takes around 10-15 minutes) 
./case.setup
./case.build

#####################################################################
# Modify namelist files
#####################################################################
# After building the model, there are still some things that you can change 
# Namelist files set different parameter values 

#cat <<EOF >> user_nl_mpaso
# config_use_data_icebergs = false
#EOF

#cat <<EOF >> user_nl_mpassi
# config_use_data_icebergs = false
#EOF

#####################################################################
#Model runs for 5 days, restarts after a certain number of days 

# RUN_OPTIONS
#
# SMOKE TEST:
# ./xmlchange -file env_run.xml -id STOP_OPTION -val ndays
# ./xmlchange -file env_run.xml -id STOP_N -val 5
# ./xmlchange -file env_run.xml -id REST_OPTION -val ndays
# ./xmlchange -file env_run.xml -id REST_N -val 1
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_QUEUE -val debug
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "00:30:00"
#
# FULL RUN:
# ./xmlchange -file env_run.xml -id STOP_OPTION -val nyears
# ./xmlchange -file env_run.xml -id STOP_N -val 5
# ./xmlchange -file env_run.xml -id REST_OPTION -val nyears
# ./xmlchange -file env_run.xml -id REST_N -val 1
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_QUEUE -val regular
# ./xmlchange -file env_batch.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "10:00:00"
#
# AFTER FIRST FULL RUN:
# ./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE
# ./xmlchange -file env_run.xml -id RESUBMIT -val 8
#
#####################################################################
#Restarts lets you continue to run models for extended periods of times while abiding by nersc que rules 
./xmlchange -file env_run.xml -id DOUT_S_ROOT --val "/global/cscratch1/sd/$USER/e3sm_scratch/$MACHINE/$CASE/archive"

./xmlchange -file env_run.xml -id BUDGETS -val TRUE

./xmlchange -file env_run.xml -id STOP_OPTION -val nmonths
./xmlchange -file env_run.xml -id STOP_N -val 3
./xmlchange -file env_run.xml -id REST_OPTION -val nmonths
./xmlchange -file env_run.xml -id REST_N -val 1
#./xmlchange -file env_run.xml -id HIST_OPTION -val nyears
#./xmlchange -file env_run.xml -id HIST_N -val 5

#regular queue for 30 minutes for 5 days of simulation  i
#regular queue for 13 horus for 3 years of simulation
./xmlchange -file env_workflow.xml -subgroup case.run -id JOB_QUEUE -val regular
#./xmlchange -file env_workflow.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "13:00:02"
./xmlchange -file env_workflow.xml -subgroup case.run -id JOB_WALLCLOCK_TIME -val "03:30:00"



./case.submit

#####################################################################
# Move to run directory to check on job
#
# cd $RUN_ROOT/$E3SM_CASE/run
#
#####################################################################
