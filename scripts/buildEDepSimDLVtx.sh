# Build edep-sim
cd $MY_TEST_AREA
git clone https://github.com/ClarkMcGrew/edep-sim
cd edep-sim
source setup.sh
export CMAKE_PREFIX_PATH=${EDEP_ROOT}/${EDEP_TARGET}
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${EDEP_ROOT}/install \
-DGeant4_DIR="/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1" ..
make -j4 install

# PandoraPFA (cmake setup and .clang-format file)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraPFA.git
cd PandoraPFA
git checkout $PANDORA_PFA_VERSION

# PandoraSDK (Abstract interface and software development kit)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraSDK.git
cd PandoraSDK
git checkout $PANDORA_SDK_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH=$MY_TEST_AREA/PandoraPFA/cmakemodules ..
make -j4 install

# PandoraMonitoring
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraMonitoring.git
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
git checkout $EIGEN_VERSION
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$MY_TEST_AREA/Eigen3/ ..
make -j4 install

# LArContent (Algorithms) with LibTorch (DLVtx)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArContent.git
cd LArContent
git checkout $PANDORA_LAR_CONTENT_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring \
-DPANDORA_LIBTORCH=ON \
-DCMAKE_PREFIX_PATH="/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake;/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib/cmake/protobuf" \
-DEigen3_DIR=$MY_TEST_AREA/Eigen3/share/eigen3/cmake/ ..
make -j4 install

# LArRecoND with edep-sim (Geant4) and LibTorch (DLVtx)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArRecoND.git
cd LArRecoND
git checkout $PANDORA_LAR_RECO_ND_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK/ \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring/ \
-DLArContent_DIR=$MY_TEST_AREA/LArContent/ \
-DPANDORA_LIBTORCH=ON \
-DCMAKE_PREFIX_PATH="/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake;/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib/cmake/protobuf" \
-DLArDLContent_DIR=$MY_TEST_AREA/LArContent/ \
-DGeant4_DIR="/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1" \
-DUSE_EDEPSIM=ON -DEDepSim_DIR=$MY_TEST_AREA/edep-sim/build ..
make -j4 install

# LArMachineLearningData (for BDT files etc)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArMachineLearningData.git
cd LArMachineLearningData
git checkout $PANDORA_LAR_MLDATA_VERSION
# Download training files from CERNBox
#. download.sh dunend
#. download.sh dune lbl
#. download.sh uboone

cd $MY_TEST_AREA
