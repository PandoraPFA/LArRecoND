# LArRecoND

Standalone Pandora application for developing and running DUNE ND reconstruction.

## Building Pandora with LArRecoND

The [build.sh](scripts/build.sh) script contains a recipe for building LArRecoND with all of the required
[Pandora](https://github.com/PandoraPFA) packages, based on the instructions from
[PandoraPFA/Documentation](https://github.com/PandoraPFA/Documentation#2-using-cmake-for-each-individual-package),
using the versions defined in [tags.sh](scripts/tags.sh). This just requires the [ROOT](https://root.cern/install)
software to be installed on the system or available using an appropriate
[CVMFS](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html) repository.
The build script also sets up the [Eigen](https://gitlab.com/libeigen/eigen) header library for
[LArContent](https://github.com/PandoraPFA/LArContent).

Before building the software, the Pandora package versions need to be defined by sourcing the
[tags.sh](scripts/tags.sh) script, which also accepts an optional argument to set the
$MY_TEST_AREA environment variable, which specifies the working directory where all of the packages
will be placed (which defaults to the current directory if it is not given):

```Shell
source tags.sh MyTestAreaDirPath
source build.sh
```

### Alma9 environment at FNAL

The [Alma9_FNAL.sh](scripts/Alma9_FNAL.sh) script can be used to setup the cmake, gcc and ROOT environment at FNAL.
This also defines the Pandora package versions using the [tags.sh](scripts/tags.sh) script, with an optional argument
to set the $MY_TEST_AREA environment variable (which defaults to the current directory if left out). Then the
[build.sh](scripts/build.sh) script can be used to build Pandora along with LArRecoND. Note that building with
LibTorch and/or edep-sim (Geant4 with CLHEP) is not currently possible within this environment, since the required
software versions for these extra packages are not yet available or compatible.

```Shell
source Alma9_FNAL.sh MyTestAreaDirPath
source build.sh
```

### SL7 environment in FNAL container

The [ContainerSL7_FNAL.sh](scripts/ContainerSL7_FNAL.sh) script sets up the SL7 environment using an
[apptainer](https://apptainer.org/docs/admin/main/installation.html) container on the FNAL computers:

```Shell
source ContainerSL7_FNAL.sh
source SL7_FNAL.sh MyTestAreaDirPath
source build.sh
```

Note that you cannot mix the Alma9 and SL7 environments, i.e. sourcing [Alma9_FNAL.sh](scripts/Alma9_FNAL.sh)
followed by [SL7_FNAL.sh](scripts/SL7_FNAL.sh) will give compiler and other environment errors. It is best to
always start a fresh interactive terminal session for whatever build environment you need to use.

The SL7 environment can also be used to build LArRecoND with LibTorch (used for Deep Learning Vertexing)
and/or edep-sim enabled, as described below.

To use an FNAL-related SL7 CVMFS container for building the code on your laptop or work/home unix PC, do
(along with any other extra apptainer settings you need, such as adding more comma-separated
`-B` directory paths):

```Shell
apptainer shell -B /cvmfs/ /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest/
source SL7_FNAL.sh MyTestAreaDirPath
source build.sh
```

### Building with LibTorch for using Deep Learning Vertexing

The [buildDLVtx.sh](scripts/buildDLVtx.sh) script contains the recipe for building LArRecoND with LibTorch v1.6.0
that is needed for using the Deep Learning Vertexing. This also requires building LArContent
(which contains the vertexing algorithm) with LibTorch turned on. This recipe uses the LibTorch library that is
available on CVMFS using a container with the SL7 environment on Fermilab computers (replace the first script
appropriately for your own laptop/PC container setup):

```Shell
source ContainerSL7_FNAL.sh
source SL7_FNAL.sh MyTestAreaDirPath
source buildDLVtx.sh
```

### Building with edep-sim (and LibTorch)

The [buildEDepSimDLVtx.sh](scripts/buildEDepSimDLVtx.sh) script contains the recipe for building LArRecoND with
[edep-sim](https://github.com/ClarkMcGrew/edep-sim) enabled as well as LibTorch for the Deep Learning Vertexing.
This requires building edep-sim with
[Geant4](https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/)
(and [CLHEP](https://proj-clhep.web.cern.ch/proj-clhep/)).
The build script uses compatible libraries from CVMFS using a container with the SL7 environment on Fermilab
computers (replace the first script appropriately for your own laptop/PC container setup):

```Shell
source ContainerSL7_FNAL.sh
source SL7_FNAL.sh MyTestAreaDirPath
source buildEDepSimDLVtx.sh
```

### LArMachineLearningData

Various neutrino algorithms need to use MicroBooNE/SBND and DUNE training files from the
[LArMachineLearningData](https://github.com/PandoraPFA/LArMachineLearningData) package. These need to be
downloaded from Google Drive using the
[download.sh](https://github.com/PandoraPFA/LArMachineLearningData/blob/master/download.sh) script:

```Shell
cd $MY_TEST_AREA/LArMachineLearningData
source download.sh sbnd
source download.sh dune
source download.sh dunend
```

This should only be done once for each new data set, since repeated attempts to download them will eventually fail,
owing to automatic download bandwidth restrictions imposed by Google Drive. If download problems do occur, then you
will need to wait up to 1 day (12 to 24 hours) before trying again.


## Running LArRecoND

For each new terminal/interactive session, make sure the environment is setup by first running either the
[tags.sh](scripts/tags.sh), [Alma9_FNAL.sh](scripts/Alma9_FNAL.sh) or [SL7_FNAL.sh](scripts/SL7_FNAL.sh) scripts, where
the optional MyTestAreaDirPath parameter sets the $MY_TEST_AREA environment variable (which defaults to the current
working directory if this is not provided):

```Shell
source tags.sh MyTestAreaDirPath
```

```Shell
source Alma9_FNAL.sh MyTestAreaDirPath
```

```Shell
source ContainerSL7_FNAL.sh
source SL7_FNAL.sh MyTestAreaDirPath
```

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

If you get runtime warnings about missing parameter files, then make sure they are downloaded in the
[LArMachineLearningData](#larmachinelearningdata) directory and their relative locations are specified
by the $FW_SEARCH_PATH environment variable, which is set by the [tags.sh](scripts/tags.sh) script.

### Geometry files

The ND-LAr geometry needs to be provided as a ROOT file containing the
[TGeoManager](https://root.cern.ch/doc/master/classTGeoManager.html) object, specified by the `-g` run parameter.
These can be created from [GDML](https://gdml.web.cern.ch/GDML/) files using ROOT:

```Shell
root -l
TGeoManager::Import("GeometryFile.gdml")
gGeoManager->Export("GeometryFile.root")
.q
```

The GDML files for the 2x2 ArgonCube prototype geometry are available in the
[2x2_sim/geometry](https://github.com/DUNE/2x2_sim/tree/develop/geometry) repository.

### 2x2 data

The following example can be used to run LArRecoND 3D-clustering and 2D-projection neutrino reconstruction
algorithms (without deep learning vertexing) for the first 10 events from a 2x2 data file:

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD.xml \
-r AllHitsNu -e Input2x2Data.root -g Geometry2x2.root -n 10 -N
```

where the mandatory settings `-i`, `-r`, `-e` and `-g` specify the xml steering run file, reconstruction hit method,
the input data ROOT file containing the hits and the geometry ROOT file, respectively. The `-n` option sets the
number of events (in this case 10) while `-N` prints out event information.

The input hit data ROOT file uses the default [SpacePoint](include/LArSP.h) format, which needs to be previously
converted from the original HDF5 format by the [ndlarflow/h5_to_root_ndlarflow.py](ndlarflow/h5_to_root_ndlarflow.py) script:

```Shell
python ndlarflow/h5_to_root_ndlarflow.py inputHDF5File 1 outputDir
```

where the first argument is the input HDF5 file, the second integer specifies data (1) and the final argument
is the location of the output directory which will store the equivalent ROOT file.

To use deep learning vertexing (DLVtx), make sure LArRecoND and LArContent is first built with LibTorch enabled, then use
the [PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) settings file.
You can tell if the DL vertexing is running if you see messages such as
`Loaded the TorchScript model PandoraNetworkDataFileName.pt` when the first event is getting processed.

### 2x2 simulation

The following example can be used to run LArRecoND 3D-clustering and 2D-projection neutrino reconstruction
algorithms (without deep learning vertexing) for the first 10 events from a 2x2 Monte Carlo (MC) file:

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD.xml \
-r AllHitsNu -e Input2x2MC.root -g Geometry2x2.root -f SPMC -n 10 -N
```

where the mandatory settings `-i`, `-r`, `-e` and `-g` specify the xml steering run file, reconstruction hit method,
the input MC ROOT file containing the hits and the geometry ROOT file, respectively. The `-f SPMC` option sets the
input to use the [SpacePoint MC](include/LArSPMC.h) format, which stores all of the MC truth information;
this is not done by the default `-f SP` format option (the ROOT data structures are different).
The `-n` option sets the number of events (in this case 10) while `-N` prints out event information.

The input MC ROOT file needs to be previously converted from the original HDF5 format by the
[ndlarflow/h5_to_root_ndlarflow.py](ndlarflow/h5_to_root_ndlarflow.py) script:

```Shell
python ndlarflow/h5_to_root_ndlarflow.py inputHDF5File 0 outputDir
```

where the first argument is the input HDF5 file, the second integer specifies MC (0) and the final argument
is the location of the output directory which will store the equivalent ROOT file.

To use deep learning vertexing (DLVtx), make sure LArRecoND and LArContent is first built with LibTorch enabled, then use
the [PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) settings file.

### edep-sim

The following example can be used to run LArRecoND 3D-clustering and 2D-projection neutrino reconstruction
algorithms (without deep learning vertexing) for the first 10 events from an edep-sim MC file:

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_ThreeD.xml \
-r AllHitsNu -e EDepSimMC.root -g EDepSimMC.root -f EDepSim -n 10 -N
```

where the mandatory settings `-i`, `-r`, `-e` and `-g` specify the xml steering run file, reconstruction hit method,
the input ROOT file containing the hits and the geometry ROOT file, respectively. The `-f EDepSim` option sets the
input to use the [edep-sim format](https://github.com/ClarkMcGrew/edep-sim/tree/master/io), which also stores all of the
available MC truth information. The `-n` option sets the number of events (in this case 10) while `-N` prints out
event information. Usually, the TGeoManager geometry information is stored in the event input ROOT file, so the same
filename should be used for both the `-e` and `-g` options if this is indeed the case.

To use deep learning vertexing (DLVtx), make sure LArRecoND and LArContent is first built with LibTorch enabled, then use
the [PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) settings file.

### Event displays

Pandora uses ROOT's [TEve](https://root.cern/doc/master/group__TEve.html) module for event displays in monitoring algorithms such as
[LArVisualMonitoring](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.cc#L364).
This is enabled in the [PandoraSettings_LArRecoND_ThreeD.xml](settings/PandoraSettings_LArRecoND_ThreeD.xml) settings
file, for example. Calling the `LArVisualMonitoring` algorithm at specific locations in the xml file will run the event display
at that point in the reconstruction algorithm flow. To disable the event display (e.g. to run in batch jobs or if there are display
problems with ROOT), comment out or remove the visual monitoring calls in the xml run file, or set the global
`IsMonitoringEnabled` variable to false (which also disables the ROOT output from the hierarchy validation tools):

```xml
    <IsMonitoringEnabled>false</IsMonitoringEnabled>
```

### Hierarchy Tools validation and analysis output

The [HierarchyAnalysisAlgorithm.cc](src/HierarchyAnalysisAlgorithm.cc) class uses
[Hierarchy Tools](https://github.com/PandoraPFA/Documentation/blob/master/Hierarchy_Tools/Hierarchy_Tools_Overview.pdf)
to create an output ROOT file that contains summary information about the reconstructed Particle Flow Objects (PFOs)
and their best-matched MC particles. The hierarchy structure and logic is implemented by LArContent's
[LArHierarchyHelper](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArHelpers/LArHierarchyHelper.h)
class. The hierarchy analysis algorithm is enabled using the following example xml settings:

```xml
   <algorithm type = "LArHierarchyAnalysis">
        <CaloHitListName>CaloHitList2D</CaloHitListName>
        <PfoListName>RecreatedPfos</PfoListName>
        <AnalysisFileName>LArRecoND.root</AnalysisFileName>
        <AnalysisTreeName>LArRecoND</AnalysisTreeName>
        <FoldToPrimaries>true</FoldToPrimaries>
        <MinPurity>0.5</MinPurity>
        <MinCompleteness>0.1</MinCompleteness>
        <MinRecoHits>15</MinRecoHits>
        <MinRecoHitsPerView>5</MinRecoHitsPerView>
        <MinRecoGoodViews>2</MinRecoGoodViews>
        <RemoveRecoNeutrons>true</RemoveRecoNeutrons>
    </algorithm>

```

This creates the [TTree](https://root.cern.ch/doc/master/classTTree.html) `LArRecoND` in the output ROOT file `LArRecoND.root`
using the PFOs stored in Pandora's `RecreatedPfos` list along with the list of hits named `CaloHitList2D` (currently the
hierarchy tools can only use the 2D views). The hierarchy building and matching requires minimum quality selection criteria,
removes neutrons and folds all of the hierarchy to start from the initial neutrino primaries.

The hierarchy analysis algorithm sets the event number by incrementing the number of times the `Run()` function is called
(0 to N-1 for N events). If the `-e` input file contains event numbers that are not contiguous, then the following xml
parameter settings (which must be added to the previous ones) need to be included to set the event numbers correctly
(which include settings for also storing the run numbers and event timing information):

```xml
    <algorithm type = "LArHierarchyAnalysis">
        <EventFileName>EventFile.root</EventFileName>
        <EventTreeName>events</EventTreeName>
        <EventLeafName>event</EventLeafName>
	<RunLeafName>run</RunLeafName>
	<SubRunLeafName>subrun</SubRunLeafName>
	<UnixTimeLeafName>unix_ts</UnixTimeLeafName>
	<StartTimeLeafName>event_start_t</StartTimeLeafName>
	<EndTimeLeafName>event_end_t</EndTimeLeafName>
        <EventsToSkip>0</EventsToSkip>
    </algorithm>
```

Here, `EventFileName` needs to match the input file name specified by the `-e` run parameter, `EventTreeName` defines what TTree
contains the event numbers (which defaults to `events`) and `EventLeafName` defines the name of the event number variable
(which defaults to `event`). This workaround is needed since it is currently not possible to pass event (and run) information
between Pandora algorithms. By default, no events are skipped, but if the `-s` run option is used, then `EventsToSkip` must be equal
to this integer to ensure that the correct event numbers are found. The additional parameters `RunLeafName` and `SubRunLeafName`
define the run and sub-run numbers, `UnixTimeLeafName` defines the trigger unix time variable name (units in seconds), while
`StartTimeLeafName` and `EndTimeLeafName` define the decimal start and end time variable names (units in ticks = 0.1 microseconds),
respectively. If these are not provided, then they will all be set to zero.

The xml settings files [PandoraSettings_LArRecoND_ThreeD.xml](settings/PandoraSettings_LArRecoND_ThreeD.xml) and
[PandoraSettings_LArRecoND_ThreeD_DLVtx.xml](settings/PandoraSettings_LArRecoND_ThreeD_DLVtx.xml) contain
(commented out) examples of using LArContent's
[LArHierarchyValidation](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h)
algorithm, which can be used to create event and particle-level ROOT output files that contain much more detail of the hierarchy than the
above LArRecoND analysis algorithm, such as complete lists of all possible reco-MC matches.
The xml files also contain (commented out) examples of using the MicroBooNE validation algorithm
[LArNeutrinoEventValidation](https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArMonitoring/NeutrinoEventValidationAlgorithm.h),
which only works for events containing single neutrino interactions (with cosmic rays).


## Fermigrid jobs

The template python script [createFNALJobs.py](scripts/createFNALJobs.py) can be used to submit LArRecoND jobs
(in SL7 containers) on the [Fermigrid](https://dune.github.io/computing-basics/07-grid-job-submission/index.html)
batch system.

It has example settings for submitting MiniRun4 or edep-sim reconstruction jobs, and it should be
relatively straightforward to extend or modify it to deal with other event samples. It uses objects
to define the setup, geometry and reconstruction parameters, which change depending on the number and
format of the sample input files. The required Pandora packages and xml steering files are stored in a tarball
by the script, which is then used by each reconstruction job, along with a copy of the input data file.
The job directories use the `/pnfs/dune/scratch/users/$USER` area, which is visible to all of the batch nodes.
Each job copies the tarball and input file to its own temporary directory area, extracts the tarball, then runs
the LArRecoND executable and copies the output files to the appropriate `/pnfs/dune/scratch/users/$USER` area.
It is recommended that the completed job output files are first copied from the scratch area to a directory in
the `/exp/dune/data/users/$USER` data area before they are used for analysis.

As an example, the following will submit reconstruction jobs for a sample of MiniRun4 input files:

```Shell
python createFNALJobs.py --option MiniRun4
source runJobs_MiniRun4.sh
```

The job run file created by the python script depends on the sample option and the number of input files.
Separate job run files are also made for each sample, which can be sourced individually to split up the
job submission process.


## License and Copyright

LArRecoND is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

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
