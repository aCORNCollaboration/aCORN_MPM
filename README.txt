###################################################################
# aCORN analysis toolkit
# by Michael P. Mendenhall (NIST 2014--2015)


------------------------------------
---------- Dependencies ------------

MPMUtils https://github.com/mpmendenhall/MPMUtils


----------------------------------
---------- Environment -----------

export MPMUTILS=${HOME}/Applications/MPMUtils/

export ACORN_DATA=${HOME}/Data/aCORN/
export ACORN_RAWDAT=${ACORN_DATA}/Raw_PIXIE/
export ACORN_PULSEDAT=${ACORN_DATA}/PulseDat/
export ACORN_REDUCED_DATA=${ACORN_DATA}/Reduced/
export ACORN_REDUCED_CALDATA=${ACORN_DATA}/Reduced_Calib/
export ACORN_REDUCED_ROOT=${ACORN_DATA}/ROOT/
export ACORN_DB=${ACORN_DATA}/aCORN_info.db
export ACORN_METADATA=${HOME}/Documents/aCORN/MetaData_06102014.csv
export ACORN_WISHBONE=${ACORN_DATA}/Wishbone/
export ACORN_SUMMARY=${HOME}/Documents/aCORN/Analyzed/
export ACORN_AUX=${HOME}/Applications/aCORN_MPM/Aux/
export ACORN_MCOUT=${ACORN_DATA}/GeantMC/

--------------------------
--------- Build ----------

make clean
make -j4
make standalone -j4


----------------------------
--------- Workflow ---------

#############
# convert .red, .rd2 files to ROOT format
cd Scripts/launchers
./ReducedToROOT.py

#############
# Analyze one wishbone series
./WishboneScanner <series>

#############
# Analyze all wishbone data
cd Scripts/launchers
./WishboneAnalysis.py --analyze

#############
# Re-draw wishbone analyzer plots without full re-analysis
cd Scripts/launchers
./WishboneAnalysis.py --replot

#############
# generate wishbone summary report
cd Scripts/report/
./ReportGen.py

#############
# example background-subtracted source set
make standalone -j4
./PMT_Gainmatcher
