# Setup DUNE environment & Pandora software tags
source envDUNE.sh
# Setup edep-sim
cd $MY_TEST_AREA/edep-sim
source setup.sh
export LD_LIBRARY_PATH=${EDEP_ROOT}/install/lib:$LD_LIBRARY_PATH
cd $MY_TEST_AREA
