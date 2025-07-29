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

  std::cout << "Runninng ProcessEvents on " << nEntries << " entries with pixel pitch " 
	    << parameters.pixelPitch << " and track/shower separation score of " << parameters.trackScoreCut << std::endl;

  // Create the class where we'll store the output info
  NDRecoOutputData fOut( "LArRecoND_test.root", parameters.runTrackFit );

  // Loop events
  for ( long entryIdx = 0; entryIdx < nEntries; ++entryIdx ) {
    pandoraIn->GetEntry(entryIdx);

    // Fill up the branches of basic output
    fOut.FillBasicBranches(pandoraIn);

    // Track fit vectors of importance
    std::vector<float> trkStartX, trkStartY, trkStartZ, trkEndX, trkEndY, trkEndZ;
    std::vector<float> trkStartDirX, trkStartDirY, trkStartDirZ, trkEndDirX, trkEndDirY, trkEndDirZ;
    std::vector<float> trkLen;

    std::vector<int> trackFitSliceId, trackFitPfoId;
    std::vector<float> trackFitQ, trackFitRR, trackFitdx, trackFitdQdx;

    // Loop particles in the event
    unsigned int nParticles = pandoraIn->m_clusterID->size();
    for ( unsigned int particleIdx=0; particleIdx < nParticles; ++particleIdx ) {
      float trackScore = pandoraIn->m_trackScore->at(particleIdx);

      if ( !parameters.runTrackFit && !parameters.runShowerFit ) continue;

      int hitCounter(0);

      // Read in the vertex and point vector that will be the input to the track and shower fits
      int sliceID = pandoraIn->m_sliceID->at(particleIdx);
      int clusterID = pandoraIn->m_clusterID->at(particleIdx);
      CartesianVector vertexVector(pandoraIn->m_nuVtxX->at(particleIdx), pandoraIn->m_nuVtxY->at(particleIdx), pandoraIn->m_nuVtxZ->at(particleIdx));
      CaloHitList caloHitList;
      for ( unsigned int idxHits = 0; idxHits < pandoraIn->m_recoHitId->size(); ++idxHits ) {
	if ( pandoraIn->m_recoHitSliceId->at(idxHits)==sliceID && pandoraIn->m_recoHitClusterId->at(idxHits)==clusterID ) {
	  CartesianVector thisHit(pandoraIn->m_recoHitX->at(idxHits), pandoraIn->m_recoHitY->at(idxHits), pandoraIn->m_recoHitZ->at(idxHits));
	  lar_content::LArCaloHitParameters chParams;
	  chParams.m_positionVector = thisHit;
	  chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
	  chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
	  chParams.m_cellGeometry = pandora::RECTANGULAR;
	  chParams.m_cellSize0 = parameters.pixelPitch;
	  chParams.m_cellSize1 = parameters.pixelPitch;
	  chParams.m_cellThickness = parameters.pixelPitch;
	  chParams.m_nCellRadiationLengths = 1.f;
	  chParams.m_nCellInteractionLengths = 1.f;
	  chParams.m_time = 0.f;
	  chParams.m_inputEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_mipEquivalentEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_electromagneticEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_hadronicEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_isDigital = false;
	  chParams.m_hitType = pandora::TPC_3D;
	  chParams.m_hitRegion = pandora::SINGLE_REGION;
	  chParams.m_layer = 0;
	  chParams.m_isInOuterSamplingLayer = false;
	  chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter));
	  chParams.m_larTPCVolumeId = 0;
	  chParams.m_daughterVolumeId = 0;
	  // push back the calo hit
	  lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
	  caloHitList.push_back( ch );
	}
      } // loop hits

      // Fit it as a track?
      if ( parameters.runTrackFit && (parameters.trackScoreCut < 0. || trackScore >= parameters.trackScoreCut) ) {
	// Run the track fit info:
	// TODO: Make the MinTrajectoryPoints(default=2) and SlidingFitHalfWindow(20) configurable
	int minTrajectoryPoints = 2;
	float slidingFitHalfWindow = 20;

	lar_content::LArTrackStateVector trackStateVector;
	std::vector<int> indexVector;
	bool trackStateSuccess=false;
	try {
	  lar_content::LArPfoHelper::GetSlidingFitTrajectory( &caloHitList, vertexVector, slidingFitHalfWindow, parameters.pixelPitch, trackStateVector, &indexVector, true );
	  trackStateSuccess=true;
	}
	catch (const pandora::StatusCodeException&) {
	  trackStateSuccess=false;
	}

	// Extract the track fit info
	if (!trackStateSuccess || trackStateVector.size() < minTrajectoryPoints) {
	  trkStartX.push_back(-9999.);
	  trkStartY.push_back(-9999.);
	  trkStartZ.push_back(-9999.);
	  trkStartDirX.push_back(1.);
	  trkStartDirY.push_back(0.);
	  trkStartDirZ.push_back(0.);
	  trkEndX.push_back(-9999.);
          trkEndY.push_back(-9999.);
          trkEndZ.push_back(-9999.);
          trkEndDirX.push_back(1.);
          trkEndDirY.push_back(0.);
          trkEndDirZ.push_back(0.);
	  trkLen.push_back(0.);
	}
	else {
	  const lar_content::LArTrackState& trackStateStart = trackStateVector.front();
	  trkStartX.push_back(trackStateStart.GetPosition().GetX());
	  trkStartY.push_back(trackStateStart.GetPosition().GetY());
	  trkStartZ.push_back(trackStateStart.GetPosition().GetZ());
	  trkStartDirX.push_back(trackStateStart.GetDirection().GetX());
          trkStartDirY.push_back(trackStateStart.GetDirection().GetY());
          trkStartDirZ.push_back(trackStateStart.GetDirection().GetZ());
	  const lar_content::LArTrackState& trackStateEnd = trackStateVector.back();
	  trkEndX.push_back(trackStateEnd.GetPosition().GetX());
          trkEndY.push_back(trackStateEnd.GetPosition().GetY());
          trkEndZ.push_back(trackStateEnd.GetPosition().GetZ());
          trkEndDirX.push_back(trackStateEnd.GetDirection().GetX());
          trkEndDirY.push_back(trackStateEnd.GetDirection().GetY());
          trkEndDirZ.push_back(trackStateEnd.GetDirection().GetZ());
	  float trklength = 0.;
	  // TODO: may need to factor in invalid points? Maybe make use of TrackTrajectoryPoints?
	  for (unsigned int idxPt=0; idxPt < trackStateVector.size()-1; ++idxPt) {
	    const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	    const lar_content::LArTrackState& trackStateNext = trackStateVector.at(idxPt+1);
	    trklength+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
	  }
	  trkLen.push_back(trklength);

	  // Track calorimetry --> very rough first pass basically reimplemented from other test branch:
	  float lengthSoFar = 0.;
	  for (unsigned int idxPt=0; idxPt < trackStateVector.size()-1; ++idxPt ) {
	    const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	    const lar_content::LArTrackState& trackStateNext = trackStateVector.at(idxPt+1);
	    lengthSoFar+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
	    // Does not do any lifetime, spacecharge, diffusion, etc. corrections at least yet
	    if ( idxPt > 0 ){
	      const lar_content::LArTrackState& trackStatePrev = trackStateVector.at(idxPt-1);
	      float hitQ = trackState.GetCaloHit()->GetInputEnergy();
	      float hitRR = trklength - lengthSoFar;
	      float hitdx = std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
	      trackFitSliceId.push_back(sliceID);
	      trackFitPfoId.push_back(clusterID);
	      trackFitQ.push_back(hitQ);
	      trackFitRR.push_back(hitRR);
	      trackFitdx.push_back(hitdx);
	      float hitdQdx = hitdx > 0. ? hitQ/hitdx : -5.f;
	      trackFitdQdx.push_back(hitdQdx);
	    }
	  }
	}
      } // TRACK FIT

      if ( parameters.runShowerFit && (parameters.trackScoreCut < 0. || trackScore < parameters.trackScoreCut) ) {
	std::cout << "I would have fit this as a shower..." << std::endl;
      } // SHOWER FIT

    } // loop particles

    if ( parameters.runTrackFit ) {
      fOut.FillTrackBranches(trkStartX,trkStartY,trkStartZ,trkStartDirX,trkStartDirY,trkStartDirZ,trkEndX,trkEndY,trkEndZ,trkEndDirX,trkEndDirY,trkEndDirZ,trkLen);
      fOut.FillTrackCaloBranches(trackFitSliceId,trackFitPfoId,trackFitQ,trackFitRR,trackFitdx,trackFitdQdx);
    }

    // Write our branches to the output tree
    fOut.WriteToFile();
  } // loop entries

  // Close our output file
  fOut.CloseFile();

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ParameterStruct &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int cOpt(0);

    parameters.runTrackFit = false;
    parameters.runShowerFit = false;
    parameters.pixelPitch = 4.0;

    while ((cOpt = getopt(argc, argv, "p:c:f:tsh")) != -1)
    {
      switch (cOpt)
      {
        case 'p':
	  parameters.pixelPitch = atof(optarg);
	  break;
        case 't':
	  parameters.runTrackFit = true;
	  break;
        case 's':
	  parameters.runShowerFit = true;
	  break;
        case 'c':
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
              << "./bin/PandoraOuterface -f [file name] -p [pixel pitch] -c [track score cut] -t -s" << std::endl
	      << "    -t = run track fit" << std::endl
	      << "    -s = run shower fit"
              << std::endl;

    return false;
}

} // namespace lar_nd_postreco
