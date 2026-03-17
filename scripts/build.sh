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

# LArContent (Algorithms) without LibTorch (no DLVtx)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArContent.git
cd LArContent
git checkout $PANDORA_LAR_CONTENT_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring \
-DEigen3_DIR=$MY_TEST_AREA/Eigen3/share/eigen3/cmake/ ..
make -j4 install

# LArRecoND
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArRecoND.git
cd LArRecoND
git checkout $PANDORA_LAR_RECO_ND_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK/ \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring/ \
-DLArContent_DIR=$MY_TEST_AREA/LArContent ..
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
