# Build LArRecoND with Pandora and other required packages.
# MY_TEST_AREA specifies the working base directory (default = $PWD)

# PandoraPFA (Basic cmake setup)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/PandoraPFA.git
cd PandoraPFA
git checkout $PANDORA_PFA_VERSION


# PandoraSDK (Abstract interface and design)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/PandoraSDK.git
cd PandoraSDK
git checkout $PANDORA_SDK_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH=$MY_TEST_AREA/PandoraPFA/cmakemodules ..
make -j4 install


# PandoraMonitoring
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/PandoraMonitoring.git
cd PandoraMonitoring
git checkout $PANDORA_MONITORING_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK ..
make -j4 install


# Eigen
cd $MY_TEST_AREA
git clone https://gitlab.com/libeigen/eigen.git Eigen3
cd Eigen3
git checkout 3.4.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$MY_TEST_AREA/Eigen3/ ..
make -j4 install


# LArContent (Algorithms)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/LArContent.git
cd LArContent
git checkout $PANDORA_LAR_CONTENT_VERSION
mkdir build
cd build
# Need to include cvmfs versions of compatible libtorch & protobuf for DL vertexing.
# If you need to specify the cmake locations of these, then add the following cmake flag:
#-DCMAKE_PREFIX_PATH="/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake;/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib/cmake/protobuf"

cmake -DCMAKE_CXX_STANDARD=17 \
-DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring \
-DPANDORA_LIBTORCH=ON \
-DEigen3_DIR=$MY_TEST_AREA/Eigen3/share/eigen3/cmake/ ..
make -j4 install


# Build edep-sim (needs Geant4 & ROOT) for LArRecoND input files
cd $MY_TEST_AREA
git clone https://github.com/ClarkMcGrew/edep-sim
cd edep-sim
source setup.sh
export CMAKE_PREFIX_PATH=${EDEP_ROOT}/${EDEP_TARGET}
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${EDEP_ROOT}/install \
-DGeant4_DIR=/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1 ..
make -j4 install


# LArRecoND
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArRecoND.git
cd LArRecoND
git checkout master
mkdir build
cd build
# Need to specify LArContent, LArDLContent, Geant4, ROOT & edep-sim. May also need librotch & protobuf:
#-DCMAKE_PREFIX_PATH="/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake;/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib/cmake/protobuf"

cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK/ \
-DPANDORA_LIBTORCH=ON -DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring/ \
-DLArContent_DIR=$MY_TEST_AREA/LArContent/ -DLArDLContent_DIR=$MY_TEST_AREA/LArContent/ \
-DGeant4_DIR=/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1 \
-DEDepSim_DIR=$MY_TEST_AREA/edep-sim/build ..
make -j4 install


# LArMachineLearningData (for BDT/DLVtx data files etc)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/LArMachineLearningData.git
cd LArMachineLearningData
git checkout $PANDORA_LAR_RECO_VERSION
# Download training data files from google drive
#source download.sh sbnd
#source download.sh dune
#source download.sh dunend

cd $MY_TEST_AREA
