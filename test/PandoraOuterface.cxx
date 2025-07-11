/**
 *  @file   LArRecoND/test/PandoraOuterface.cc
 *
 *  @brief  Implementation of the Post-Pandora high-level reco for DUNE ND
 *
 *  $Log: $
 */

#include "TFile.h"
#include "TTree.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

#include "Api/PandoraApi.h"
#include "Geometry/LArTPC.h"
#include "Helpers/XmlHelper.h"
#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

#include "LArNDContent.h"
#include "LArNDGeomSimple.h"
#include "LArRay.h"
#include "PandoraOuterface.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace pandora;
using namespace lar_nd_postreco;

int main(int argc, char *argv[])
{

  int errorNo(0);

  try
  {
    ParameterStruct pset;

    if (!ParseCommandLine(argc, argv, pset))
      return 1;

    ProcessPostReco(pset);
  }
  catch (const StatusCodeException &statusCodeException)
  {
    std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
    errorNo = 1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception: " << std::endl;
    errorNo = 1;
  }

  return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_postreco
{

void ProcessPostReco(const ParameterStruct &parameters)
{
  TFile *fileSource = TFile::Open(parameters.fileName.c_str(), "READ");
  if (!fileSource)
  {
      std::cout << "Error in ProcessPostReco(): can't open file " << parameters.fileName << std::endl;
      return;
  }

  TTree *recoTree = dynamic_cast<TTree *>(fileSource->Get("LArRecoND"));
  if (!recoTree)
  {
      std::cout << "Could not find the event tree LArRecoND" << std::endl;
      fileSource->Close();
      return;
  }

  std::unique_ptr<LArRecoNDFormat> pandoraIn = std::make_unique<LArRecoNDFormat>(recoTree);

  long nEntries = recoTree->GetEntries();

  std::cout << "Runninng ProcessEvents on " << nEntries << " entries with pixel pitch " << parameters.pixelPitch << " and track/shower separation score of " << parameters.trackScoreCut << std::endl;

  // Loop events
  for ( long entryIdx = 0; entryIdx < nEntries; ++entryIdx ) {
    pandoraIn->GetEntry(entryIdx);

    // Loop particles in the event
    unsigned int nParticles = pandoraIn->m_clusterID->size();
    for ( unsigned int particleIdx=0; particleIdx < nParticles; ++particleIdx ) {
      float trackScore = pandoraIn->m_trackScore->at(particleIdx);

      // Fitting it as a track
      if ( parameters.trackScoreCut < 0. || trackScore >= parameters.trackScoreCut ) {
	int sliceID = pandoraIn->m_sliceID->at(particleIdx);
	int clusterID = pandoraIn->m_clusterID->at(particleIdx);
	// Vertex
	CartesianVector vertexVector(pandoraIn->m_nuVtxX->at(particleIdx), pandoraIn->m_nuVtxY->at(particleIdx), pandoraIn->m_nuVtxZ->at(particleIdx));
	// Make the Cartesian Points
	CartesianPointVector pointVector;
	for ( unsigned int idxHits = 0; idxHits < pandoraIn->m_recoHitId->size(); ++idxHits ) {
	  if ( pandoraIn->m_recoHitSliceId->at(idxHits)==sliceID && pandoraIn->m_recoHitClusterId->at(idxHits)==clusterID ) {
	    CartesianVector thisHit(pandoraIn->m_recoHitX->at(idxHits), pandoraIn->m_recoHitY->at(idxHits), pandoraIn->m_recoHitZ->at(idxHits));
	    pointVector.push_back(thisHit);
	  }
	} // loop hits

	// Run the track fit info:
	// TODO: Make the MinTrajectoryPoints(default=2) and SlidingFitHalfWindow(20) configurable
	int minTrajectoryPoints = 2;
	float slidingFitHalfWindow = 20;

	lar_content::LArTrackStateVector trackStateVector;
	bool trackStateSuccess=false;
	try {
	  lar_content::LArPfoHelper::GetSlidingFitTrajectory( pointVector, vertexVector, slidingFitHalfWindow, parameters.pixelPitch, trackStateVector );
	  trackStateSuccess=true;
	}
	catch (const pandora::StatusCodeException&) {
	  trackStateSuccess=false;
	  std::cout << "Unable to extract sliding fit trajectory" << std::endl;
	}

	// Extract the track fit info
	if (!trackStateSuccess || trackStateVector.size() < minTrajectoryPoints) {
	  std::cout << " ------> TRACK FIT FAILURE" << std::endl;
	}
	else {
	  const lar_content::LArTrackState& trackStateStart = trackStateVector.front();
	  std::cout << "Track start point (" << trackStateStart.GetPosition().GetX() << ", " << trackStateStart.GetPosition().GetY() << ", " << trackStateStart.GetPosition().GetZ() << ")"
		    << ", and direction (" << trackStateStart.GetDirection().GetX() << ", " << trackStateStart.GetDirection().GetY() << ", " << trackStateStart.GetDirection().GetZ() << ")." << std::endl;
	  const lar_content::LArTrackState& trackStateEnd = trackStateVector.back();
	  std::cout << "Track end point (" << trackStateEnd.GetPosition().GetX() << ", " << trackStateEnd.GetPosition().GetY() << ", " << trackStateEnd.GetPosition().GetZ() << ")"
                    << ", and direction (" << trackStateEnd.GetDirection().GetX() << ", " << trackStateEnd.GetDirection().GetY() << ", " << trackStateEnd.GetDirection().GetZ() << ")." << std::endl;
	  float trklength = 0.;
	  // TODO: factor in invalid points. Maybe make use of TrackTrajectoryPoints?
	  for (unsigned int idxPt=0; idxPt < trackStateVector.size()-1; ++idxPt) {
	    const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	    const lar_content::LArTrackState& trackStateNext = trackStateVector.at(idxPt+1);
	    trklength+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
	  }
	  std::cout << "The track length is: " << trklength << std::endl;
	}
      } // TRACK FIT

      if ( parameters.trackScoreCut < 0. || trackScore < parameters.trackScoreCut ) {
	std::cout << "I would have fit this as a shower..." << std::endl;
      } // SHOWER FIT

      // ADD SHOWER FIT HERE

    } // loop particles
  }

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ParameterStruct &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int cOpt(0);

    parameters.pixelPitch = 4.0;

    while ((cOpt = getopt(argc, argv, "p:t:f:h")) != -1)
    {
      switch (cOpt)
      {
        case 'p':
	  parameters.pixelPitch = atof(optarg);
	  break;
        case 't':
	  parameters.trackScoreCut = atof(optarg);
	  break;
        case 'f':
	  parameters.fileName = optarg;
	  break;
        case 'h':
        default:
	  return PrintOptions();
      }
    }

    bool passed=true;
    if(!passed)
    {
        return PrintOptions();
    }
    return passed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl
              << "./bin/PandoraOuterface -f [file name] -p [pixel pitch] -t [track score cut]" << std::endl
              << std::endl;

    return false;
}

} // namespace lar_nd_postreco
