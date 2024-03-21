# LArRecoND
<!-- [![Build Status](https://travis-ci.org/PandoraPFA/LArRecoND.svg?branch=master)](https://travis-ci.org/PandoraPFA/LArRecoND) -->
<!-- [![Coverity Scan Build Status](https://scan.coverity.com/projects/13060/badge.svg)](https://scan.coverity.com/projects/pandorapfa-larrecond) -->

Standalone Pandora application for developing and testing DUNE ND reconstruction.

**Environment setup**

Make sure your computer can use [CVMFS](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html).

Setup the environment by sourcing the [envDUNE.sh](envDUNE.sh) script, which accepts an optional argument
to set the $MY_TEST_AREA environment variable that specifies the working directory where all of the
packages will be placed (defaults to the current directory):

```Shell
source envDUNE.sh MyTestAreaDirPath
```

If there are problems doing this, then try using a container such as
[apptainer](https://apptainer.org/docs/admin/main/installation.html) before running the script:

```Shell
apptainer shell -B /cvmfs/ /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest/
source envDUNE.sh MyTestAreaDirPath
```

**Building Pandora and LArRecoND software**

The [doBuild.sh](doBuild.sh) script contains a recipe for building LArRecoND and Pandora,
along with [edep-sim](https://github.com/ClarkMcGrew/edep-sim). This is based on instructions from the
[PandoraPFA/Documentation](https://github.com/PandoraPFA/Documentation#2-using-cmake-for-each-individual-package)
area. Note that extra cmake variables for LArContent and LArRecoND are needed to enable the Deep Learning (DL)
software. If different tags, branches or git forks are used, then rebuild the appropriate code by deleting and
remaking the `build` directories of the affected packages.

Various neutrino algorithms need to use MicroBooNE/SBND and DUNE training files from the
[LArMachineLearningData](https://github.com/PandoraPFA/LArMachineLearningData) package.
These need to be downloaded from google drive using the
[download.sh](https://github.com/PandoraPFA/LArMachineLearningData/blob/master/download.sh) script:

```Shell
source download.sh sbnd
source download.sh dune
source download.sh dunend
```

This should only be done once (for each data set), as repeated attempts to download them will eventually fail,
owing to download bandwidth restrictions imposed by google drive. If download warnings do occur, then waiting
around 1 day before trying again should reinstate download access.

**Running LArRecoND**

First, make sure the DUNE, Pandora and edep-sim environment is setup by sourcing the [setup.sh](setup.sh)
script for each new terminal/interactive session.

The LArRecoND software is run by using the `PandoraInterface` executable, which is created from the
[PandoraInterface.cxx](test/PandoraInterface.cxx) main program that is steered with xml files from the
[settings](settings) directory. This application uses the neutrino reconstruction methods defined in the
[LArContent](https://github.com/PandoraPFA/LArContent) package, as well as the 3D algorithms listed in
[LArNDContent.cc](src/LArNDContent.cc) along with the [MasterThreeDAlgorithm.cc](src/MasterThreeDAlgorithm.cc) class.

If everything has been built correctly, running

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -h
``` 
will list all available (required and optional) run options.

Here is an example of running the code using DL vertexing along with storing the MC and reco information
using hierarchy tools:

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_DLHierarchy.xml \
-r AllHitsNu -e LArEDepSim_numu_all_1.root -j LArTPC -N -n 10 -d ArgonCube
```

where the mandatory settings `-i`, `-r` and `-e` specify the xml steering file, reconstruction hit type
and the full name of the input file containing events (hit collections) in the default `EDepSim` format,
respectively. The option `-j` sets the hit projection method, `-N` prints out event info, `-n` sets the
number of events (in this case 10), while `-d` defines the Geant4 sensitive detector name containing
the hits we want to reconstruct.

You can tell if the DL vertexing is running if you see messages such as
`Loaded the TorchScript model PandoraNetworkDataFileName.pt` after the creation of the first event voxels.

Another example running both 3D and LArTPC reconstruction with DL vertexing and hierarchy tools on
2x2 MiniRun4 simulation data (5 events) is:

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml -j both -N -r AllHitsSliceNu \
-f SPMC -g data/MiniRun4/Merged2x2MINERvA_v3_withRock.root -k events -t Default -d volLArActive \
-v volArgonCubeDetector_PV -M -s 0 -n 5 -e data/MiniRun4/MiniRun4_1E19_RHC.flow.00001.FLOWTestMergedhits.root
```

Here, the hits are in SpacePoint MC [(SPMC)](include/LArSPMC.h) format (`-f`) from the `events` tree (`-k`)
inside the ROOT input file (`-e`). The geometry is set by the ROOT file specified by the `-g` option, along
with the physical placement geometry (`-v`) and sensitive (`-d`) detector names that belong to the
[TGeoManager](https://root.cern.ch/doc/master/classTGeoManager.html) object (`-t`) called `Default`.


LArRecoND is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## License and Copyright
Copyright (C), LArRecoND Authors

LArRecoND is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
