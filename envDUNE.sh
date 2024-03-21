# Setup DUNE environment
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# Setup various package versions
# For libxxhash.so
setup dune_oslibs v1_0_0
# Compiler
setup gcc v9_3_0
# Git
setup git v2_20_1
# ROOT
setup root v6_22_08d -q e20:p392:prof
# GEANT4 for edep-sim
setup geant4 v4_10_6_p01c -q e19:prof
# CMAKE
setup cmake v3_24_1
# Clang formatting
setup clang v7_0_0
# LibTorch for Deep Learning vertexing
setup libtorch v1_6_0d -q e20

# For Fermigrid jobs
setup ifdhc
setup jobsub_client

# Set main working directory by optional run argument
dirName=$1
testArea=$PWD

if [[ $dirName ]]; then
    testArea=$dirName
fi
echo "MY_TEST_AREA is $testArea"
export MY_TEST_AREA=${testArea}

# Set Pandora package versions
# Use "git fetch --tags" to update tag list for a package
export PANDORA_PFA_VERSION=v04-07-00
export PANDORA_SDK_VERSION=v03-04-01
export PANDORA_MONITORING_VERSION=v03-05-00
export PANDORA_LAR_CONTENT_VERSION=v04_07_00
export PANDORA_LC_CONTENT_VERSION=v03-01-06
export PANDORA_EXAMPLE_CONTENT_VERSION=v03-01-00
export PANDORA_LAR_RECO_VERSION=v04-07-00
export PANDORA_LC_RECO_VERSION=v03-01-05
export PANDORA_LAR_RECO_ND_VERSION=v01-00-00

# Set FW_SEARCH_PATH for Pandora xml run files & machine learning data etc.
export PATH=$MY_TEST_AREA/LArRecoND/bin:$PATH
export FW_SEARCH_PATH=$MY_TEST_AREA/LArRecoND/settings
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData:$FW_SEARCH_PATH

# DLVtx torch & protobuf libraries (should be added by earlier package "setup" steps)
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib
