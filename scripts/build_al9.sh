#!/bin/bash

RunWithError() {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "Error running command: $@"
        echo "Exiting with status $status"
        exit $status
    fi
    return $status
}

# Has the Pandora project dir been set?
if [ -z "$PANDORA_PROJECT_DIR" ]; then
    echo "Please set the PANDORA_PROJECT_DIR environment variable to the Pandora project directory."
    exit 1
fi

# Set default versions if not already set
PANDORA_PFA_VERSION=${PANDORA_PFA_VERSION:-"v05-00-00"}
PANDORA_SDK_VERSION=${PANDORA_SDK_VERSION:-"v05-00-00"}
PANDORA_MONITORING_VERSION=${PANDORA_MONITORING_VERSION:-"v05-00-00"}
PANDORA_LAR_CONTENT_VERSION=${PANDORA_LAR_CONTENT_VERSION:-"v05_02_01"}
PANDORA_LAR_MLDATA_VERSION=${PANDORA_LAR_MLDATA_VERSION:-"v05-02-01"}

# PandoraPFA (cmake setup and .clang-format file)
RunWithError cd $PANDORA_PROJECT_DIR
RunWithError git clone https://github.com/PandoraPFA/PandoraPFA.git
RunWithError cd PandoraPFA
RunWithError git checkout $PANDORA_PFA_VERSION

# PandoraSDK (Abstract interface and software development kit)
RunWithError cd $PANDORA_PROJECT_DIR
RunWithError git clone https://github.com/PandoraPFA/PandoraSDK.git
RunWithError cd PandoraSDK
RunWithError git checkout $PANDORA_SDK_VERSION
RunWithError mkdir build
RunWithError cd build
RunWithError cmake -DCMAKE_MODULE_PATH=$PANDORA_PROJECT_DIR/PandoraPFA/cmakemodules ..
RunWithError make -j4 install

# PandoraMonitoring
RunWithError cd $PANDORA_PROJECT_DIR
RunWithError git clone https://github.com/PandoraPFA/PandoraMonitoring.git
RunWithError cd PandoraMonitoring
RunWithError git checkout $PANDORA_MONITORING_VERSION
RunWithError mkdir build
RunWithError cd build
RunWithError cmake \
    -DCMAKE_MODULE_PATH="$PANDORA_PROJECT_DIR/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
    -DPandoraSDK_DIR=$PANDORA_PROJECT_DIR/PandoraSDK/build/install/lib64/cmake/PandoraSDK \
    ..

RunWithError make -j4 install

# LArContent (Algorithms) with LibTorch
RunWithError cd $PANDORA_PROJECT_DIR
RunWithError git clone https://github.com/PandoraPFA/LArContent.git
RunWithError cd LArContent
RunWithError git checkout $PANDORA_LAR_CONTENT_VERSION
RunWithError mkdir build
RunWithError cd build

# Find LibTorch from either FNAL Spack or a preinstalled path (e.g. CI).
if [ -z "$TORCH_PREFIX" ]; then
    if command -v spack >/dev/null 2>&1; then
        TORCH_PREFIX=$(spack location -i py-torch@1.13.1/dicc4vw 2>/dev/null)
    fi
fi

if [ -z "$TORCH_PREFIX" ] && [ -d "/pandora/libtorch" ]; then
    TORCH_PREFIX="/pandora/libtorch"
fi

if [ -z "$TORCH_PREFIX" ]; then
    echo "Unable to locate LibTorch. Set TORCH_PREFIX or provide Spack py-torch@1.13.1/dicc4vw."
    exit 1
fi

if [ -d "${TORCH_PREFIX}/share/cmake/Torch" ]; then
    TORCH_CMAKE_DIR="${TORCH_PREFIX}/share/cmake/Torch"
elif [ -d "${TORCH_PREFIX}/lib/python3.9/site-packages/torch/share/cmake/Torch" ]; then
    TORCH_CMAKE_DIR="${TORCH_PREFIX}/lib/python3.9/site-packages/torch/share/cmake/Torch"
else
    echo "Unable to locate Torch CMake config under ${TORCH_PREFIX}."
    exit 1
fi

RunWithError cmake \
    -DPANDORA_MONITORING=ON \
    -DPANDORA_LIBTORCH=ON \
    -DPandoraSDK_DIR=$PANDORA_PROJECT_DIR/PandoraSDK/build/install/lib64/cmake/PandoraSDK \
    -DPandoraMonitoring_DIR=$PANDORA_PROJECT_DIR/PandoraMonitoring/build/install/lib64/cmake/PandoraMonitoring \
    -DTorch_DIR=${TORCH_CMAKE_DIR} \
    ..
RunWithError make -j4 install

# LArRecoND
RunWithError cd $PANDORA_PROJECT_DIR
RunWithError cd LArRecoND
RunWithError mkdir build
RunWithError cd build
RunWithError cmake \
    -DPANDORA_MONITORING=ON \
    -DPANDORA_LIBTORCH=ON \
    -DPandoraSDK_DIR=$PANDORA_PROJECT_DIR/PandoraSDK/build/install/lib64/cmake/PandoraSDK \
    -DPandoraMonitoring_DIR=$PANDORA_PROJECT_DIR/PandoraMonitoring/build/install/lib64/cmake/PandoraMonitoring \
    -DLArContent_DIR=$PANDORA_PROJECT_DIR/LArContent/build/install/lib64/cmake/LArContent \
    -DLArDLContent_DIR=$PANDORA_PROJECT_DIR/LArContent/build/install/lib64/cmake/LArDLContent \
    -DTorch_DIR=${TORCH_CMAKE_DIR} ..
RunWithError make -j4 install

# LArMachineLearningData (for BDT files etc)
RunWithError cd $PANDORA_PROJECT_DIR
RunWithError git clone https://github.com/PandoraPFA/LArMachineLearningData.git
RunWithError cd LArMachineLearningData
RunWithError git checkout $PANDORA_LAR_MLDATA_VERSION
# Download training files: only do this once to avoid google drive's access restrictions (up to 24 hrs wait)
#. download.sh sbnd
#. download.sh dune
#. download.sh dunend

cd $PANDORA_PROJECT_DIR
