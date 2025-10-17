/**
 *  @file   LArRecoND/test/PandoraOuterface.cc
 *
 *  @brief  Implementation of the Post-Pandora high-level reco for DUNE ND
 *
 *  $Log: $
 */

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

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
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

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
#include <map>
#include <functional>

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

    if (!ReadSettings(pset))
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

  std::cout << "Hello! Runninng ProcessEvents on " << nEntries << " entries with pixel pitch " 
	    << parameters.pixelPitch << " and track/shower separation score of " << parameters.trackScoreCut << std::endl;

  // Create the class where we'll store the output info
  NDRecoOutputData fOut( parameters.outfileName, parameters.runTrackFit, parameters.runShowerFit );

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
    std::vector<float> trackFitX, trackFitY, trackFitZ;
    std::vector<float> trackFitQ, trackFitRR, trackFitdx, trackFitdQdx;
    
    //Shower fit vectors of importance
    std::vector<float> shwrCentroidX, shwrCentroidY, shwrCentroidZ, shwrStartX, shwrStartY, shwrStartZ;
    std::vector<float> shwrDirX, shwrDirY, shwrDirZ;
    std::vector<float> shwrLen;

    std::vector<int> shwrSliceId, shwrClusterId;

    std::vector<float> startTrkSlidingFitDirX,  startTrkSlidingFitDirY,  startTrkSlidingFitDirZ;  
    std::vector<float> startTrkPCADirX, startTrkPCADirY, startTrkPCADirZ;  
    std::vector<float> shwrStartPointsX, shwrStartPointsY, shwrStartPointsZ;
    std::vector<double> shwrdEdx;
    std::vector<float> shwrTotalE;

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
	  // Skip hit if it fails the threshold
	  if ( parameters.applyThreshold && pandoraIn->m_recoHitE->at(idxHits) < parameters.thresholdVal )
	    continue;
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

	// If user has set the Voxelize Z function, then rerun the track fit, starting from the output of the first fit
	lar_content::LArTrackStateVector trackStateVector_v2;
	std::vector<int> indexVector_v2;
        bool trackStateSuccess_v2=false;
	if ( parameters.voxelizeZ && trackStateSuccess ) {
	  if( parameters.verbosity >= 1 ) {
	    std::cout << "    INFO: Since voxelization is turned on, we will take the output of the track fit and try to voxelize now." << std::endl;
	    std::cout << "    ----> Input track has " << trackStateVector.size() << " track points." << std::endl;
	  }
	  try
	  {
	    int hitCounter_v1p5(0);
	    int hitCounter_v2(0);

	    // Initial calohit vector
	    std::vector< lar_content::LArCaloHit* > caloHitVect_v1;
	    for (unsigned int idxPt=0; idxPt < trackStateVector.size(); ++idxPt ) {
	      const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	      lar_content::LArCaloHitParameters chParams;
	      chParams.m_positionVector = trackState.GetCaloHit()->GetPositionVector();
	      chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
	      chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
	      chParams.m_cellGeometry = pandora::RECTANGULAR;
	      chParams.m_cellSize0 = parameters.pixelPitch;
	      chParams.m_cellSize1 = parameters.pixelPitch;
	      chParams.m_cellThickness = parameters.pixelPitch;
	      chParams.m_nCellRadiationLengths = 1.f;
	      chParams.m_nCellInteractionLengths = 1.f;
	      chParams.m_time = 0.f;
	      chParams.m_inputEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_mipEquivalentEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_electromagneticEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_hadronicEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_isDigital = false;
	      chParams.m_hitType = pandora::TPC_3D;
	      chParams.m_hitRegion = pandora::SINGLE_REGION;
	      chParams.m_layer = 0;
	      chParams.m_isInOuterSamplingLayer = false;
	      chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter_v1p5));
	      chParams.m_larTPCVolumeId = 0;
	      chParams.m_daughterVolumeId = 0;
	      lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
	      caloHitVect_v1.push_back( ch );
	    }
	    // Now let's construct the version that goes into the second pass track fit.
	    // 1. Loop through the vector and for each element, gather all the consecutive elements within epsilon of the z value
	    // 2. Within this subset, find the maximum Q hit, start here
	    //     a. Gather this hit and the ones within an x, y distance of the voxel setting
	    //     b. Make a new calohit that is the weighted mean of the (x, y, z) of these hits and the sum of the Q values
	    // 3. Repeat on the maximal Q value of the hits letf and continue repeating till all hits are swept up
	    // 4. Run track fit on this.
	    CaloHitList caloHitList_v2;
	    for ( unsigned int idxHit=0; idxHit < caloHitVect_v1.size(); ++idxHit ) {
	      if ( parameters.verbosity >= 2 ) std::cout << "      hit idx " << idxHit << " of " << caloHitVect_v1.size() << std::endl;
	      // Step 1
	      float thisZ = caloHitVect_v1.at(idxHit)->GetPositionVector().GetZ();
	      std::vector< lar_content::LArCaloHit* > caloHitVect_tmp;
	      caloHitVect_tmp.push_back( caloHitVect_v1.at(idxHit) );
	      bool stopLoop=false;
	      while ( !stopLoop && idxHit < caloHitVect_v1.size()-1 ) {
		if ( fabs(caloHitVect_v1.at(idxHit+1)->GetPositionVector().GetZ()-thisZ) < std::numeric_limits<float>::epsilon() ) {
		  caloHitVect_tmp.push_back( caloHitVect_v1.at(idxHit+1) );
		  idxHit+=1;
		}
		else stopLoop=true;
	      } // found all the hits that we need to check
	      if( parameters.verbosity >= 2 ) std::cout << "      --> At this stage of the voxelization, we have " << caloHitVect_tmp.size() << " hits to possibly merge." << std::endl;
	      // Step 2-3
	      if ( caloHitVect_tmp.size() == 1 ) {
		lar_content::LArCaloHitParameters chParams;
		chParams.m_positionVector = caloHitVect_tmp.at(0)->GetPositionVector();
		chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
		chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
		chParams.m_cellGeometry = pandora::RECTANGULAR;
		chParams.m_cellSize0 = parameters.pixelPitch;
		chParams.m_cellSize1 = parameters.pixelPitch;
		chParams.m_cellThickness = parameters.pixelPitch;
		chParams.m_nCellRadiationLengths = 1.f;
		chParams.m_nCellInteractionLengths = 1.f;
		chParams.m_time = 0.f;
		chParams.m_inputEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_mipEquivalentEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_electromagneticEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_hadronicEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_isDigital = false;
		chParams.m_hitType = pandora::TPC_3D;
		chParams.m_hitRegion = pandora::SINGLE_REGION;
		chParams.m_layer = 0;
		chParams.m_isInOuterSamplingLayer = false;
		chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter_v2));
		chParams.m_larTPCVolumeId = 0;
		chParams.m_daughterVolumeId = 0;
		lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
		caloHitList_v2.push_back( ch );
	      }
	      else {
		while ( caloHitVect_tmp.size() > 0 ) {
		  float maxQ = 0.;
		  float maxQ_X = 0.;
		  float maxQ_Y = 0.;
		  for ( unsigned int idxHit_inner=0; idxHit_inner < caloHitVect_tmp.size(); ++idxHit_inner ) {
		    if ( caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy() > maxQ ) {
		      maxQ = caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy();
		      maxQ_X = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetX();
		      maxQ_Y = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetY();
		    }
		  }
		  if ( parameters.verbosity >= 2 ) std::cout << "      --> Max Hit X = " << maxQ_X << ", Y = " << maxQ_Y << ", Q = " << maxQ << std::endl;
		  std::vector<float> xs, ys, zs, qs;
		  std::vector<unsigned int> toDelete;
		  for ( unsigned int idxHit_inner=0; idxHit_inner < caloHitVect_tmp.size(); ++idxHit_inner ) {
		    float thisX_inner = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetX();
		    float thisY_inner = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetY();
		    if ( parameters.verbosity >= 2 ) std::cout << "      --> This Hit X = " << thisX_inner << ", Y = " << thisY_inner << ", Q = " << caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy() <<std::endl;
                    if ( std::sqrt( std::pow(thisX_inner-maxQ_X,2) + std::pow(thisY_inner-maxQ_Y,2) ) < parameters.voxelZHW ) {
		      float thisZ_inner = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetZ();
		      float thisQ_inner = caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy();
		      xs.push_back(thisX_inner);
		      ys.push_back(thisY_inner);
		      zs.push_back(thisZ_inner);
		      qs.push_back(thisQ_inner);
		      toDelete.push_back(idxHit_inner);
                    }
                  }
		  if ( parameters.verbosity >= 2 ) std::cout << "      --> Making a new hit from " << xs.size() << " hit(s) and deleting " << toDelete.size() << " hits." << std::endl;
		  std::vector< lar_content::LArCaloHit* > caloHitVect_tmp_prev = caloHitVect_tmp;
		  caloHitVect_tmp.clear();
		  for ( unsigned int idxCopy=0; idxCopy < caloHitVect_tmp_prev.size(); ++idxCopy ) {
		    bool skipCopy = false;
		    for ( unsigned int checkIdx=0; checkIdx < toDelete.size(); ++checkIdx ) {
		      if ( idxCopy == toDelete[checkIdx] ) {
			skipCopy = true;
			break;
		      }
		    }
		    if ( skipCopy ) continue;
		    caloHitVect_tmp.push_back( caloHitVect_tmp_prev.at(idxCopy) );
		  }
		  // Make new hit:
		  float newHitX(0.), newHitY(0.), newHitZ(0.), newHitQ(0.);
		  for ( unsigned int idxUse=0; idxUse<xs.size(); ++idxUse ) {
		    newHitX+=xs[idxUse]*qs[idxUse];
		    newHitY+=ys[idxUse]*qs[idxUse];
		    newHitZ+=zs[idxUse]*qs[idxUse];
		    newHitQ+=qs[idxUse];
		  }
		  if ( newHitQ > 0. ){
		    newHitX/=newHitQ;
		    newHitY/=newHitQ;
		    newHitZ/=newHitQ;
		  }
		  lar_content::LArCaloHitParameters chParams;
		  chParams.m_positionVector = {newHitX, newHitY, newHitZ};
		  chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
		  chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
		  chParams.m_cellGeometry = pandora::RECTANGULAR;
		  chParams.m_cellSize0 = parameters.pixelPitch;
		  chParams.m_cellSize1 = parameters.pixelPitch;
		  chParams.m_cellThickness = parameters.pixelPitch;
		  chParams.m_nCellRadiationLengths = 1.f;
		  chParams.m_nCellInteractionLengths = 1.f;
		  chParams.m_time = 0.f;
		  chParams.m_inputEnergy = newHitQ;
		  chParams.m_mipEquivalentEnergy = newHitQ;
		  chParams.m_electromagneticEnergy = newHitQ;
		  chParams.m_hadronicEnergy = newHitQ;
		  chParams.m_isDigital = false;
		  chParams.m_hitType = pandora::TPC_3D;
		  chParams.m_hitRegion = pandora::SINGLE_REGION;
		  chParams.m_layer = 0;
		  chParams.m_isInOuterSamplingLayer = false;
		  chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter_v2));
		  chParams.m_larTPCVolumeId = 0;
		  chParams.m_daughterVolumeId = 0;
		  lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
		  caloHitList_v2.push_back( ch );
		  if( parameters.verbosity >= 2 )
		    std::cout << "      --> After this particular voxelization, we have " << caloHitVect_tmp.size() <<" hits remaining to possibly merge.\n"
			      << "          and caloHitList_v2 size is " << caloHitList_v2.size() << std::endl;
		}
	      } // Step 2-3
	    } // Steps 1-3
	    // Step 4
	    if( parameters.verbosity >= 1 ) std::cout << "    ----> DONE with the merging. Now running the new track fit." << std::endl;
	    lar_content::LArPfoHelper::GetSlidingFitTrajectory( &caloHitList_v2, vertexVector, slidingFitHalfWindow, parameters.pixelPitch, trackStateVector_v2, &indexVector_v2, true );

	    trackStateSuccess_v2=true;
	  }
	  catch(const pandora::StatusCodeException&) {
	    trackStateSuccess_v2=false;
	  }
	}

	lar_content::LArTrackStateVector trackStateVector_out = (parameters.voxelizeZ && trackStateSuccess_v2 ) ? trackStateVector_v2 : trackStateVector;
	if( parameters.verbosity >= 1 ) std::cout << "    INFO: The track state vector we are using for calorimetry analysis has " << trackStateVector_out.size() << " points." << std::endl;

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
	  const lar_content::LArTrackState& trackStateStart =
	    ( parameters.useVoxelizedStartStop && (trackStateSuccess_v2 && trackStateVector_out.size() >= minTrajectoryPoints) ) ? 
	    trackStateVector_out.front() : 
	    trackStateVector.front();
	  trkStartX.push_back(trackStateStart.GetPosition().GetX());
	  trkStartY.push_back(trackStateStart.GetPosition().GetY());
	  trkStartZ.push_back(trackStateStart.GetPosition().GetZ());
	  trkStartDirX.push_back(trackStateStart.GetDirection().GetX());
          trkStartDirY.push_back(trackStateStart.GetDirection().GetY());
          trkStartDirZ.push_back(trackStateStart.GetDirection().GetZ());
	  const lar_content::LArTrackState& trackStateEnd =
	    ( parameters.useVoxelizedStartStop && (trackStateSuccess_v2 && trackStateVector_out.size() >= minTrajectoryPoints) ) ?
            trackStateVector_out.back() :
            trackStateVector.back();
	  trkEndX.push_back(trackStateEnd.GetPosition().GetX());
          trkEndY.push_back(trackStateEnd.GetPosition().GetY());
          trkEndZ.push_back(trackStateEnd.GetPosition().GetZ());
          trkEndDirX.push_back(trackStateEnd.GetDirection().GetX());
          trkEndDirY.push_back(trackStateEnd.GetDirection().GetY());
          trkEndDirZ.push_back(trackStateEnd.GetDirection().GetZ());
	  float trklength = 0.;
	  // Get the length going point to point
	  if ( parameters.useVoxelizedStartStop && (trackStateSuccess_v2 && trackStateVector_out.size() >= minTrajectoryPoints) ) {
	    for (unsigned int idxPt=0; idxPt < trackStateVector_out.size()-1; ++idxPt) {
              const lar_content::LArTrackState& trackState = trackStateVector_out.at(idxPt);
              const lar_content::LArTrackState& trackStateNext = trackStateVector_out.at(idxPt+1);
              trklength+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
            }
	  }
	  else {
	    for (unsigned int idxPt=0; idxPt < trackStateVector.size()-1; ++idxPt) {
	      const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	      const lar_content::LArTrackState& trackStateNext = trackStateVector.at(idxPt+1);
	      trklength+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
	    }
	  }
	  trkLen.push_back(trklength);

	  // Track calorimetry --> very rough first pass basically reimplemented from other test branch:
	  // ! Consider the first and last points, but here we only have one side of dx
	  // ! Does not do any lifetime, spacecharge, diffusion, etc. corrections at least yet
	  if ( trackStateVector_out.size() >= minTrajectoryPoints ) {
	    float lengthSoFar = 0.;
	    for (unsigned int idxPt=0; idxPt < trackStateVector_out.size(); ++idxPt ) {
	      const lar_content::LArTrackState& trackState = trackStateVector_out.at(idxPt);
	      // First point
	      if ( idxPt == 0 ){
		float hitQ = trackState.GetCaloHit()->GetInputEnergy();
		float hitRR = trklength;
		float hitdx = 0.;
		if ( idxPt < trackStateVector_out.size()-1 ) {
		  const lar_content::LArTrackState& trackStateNext = trackStateVector_out.at(idxPt+1);
		  hitdx = std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
		}
		trackFitSliceId.push_back(sliceID);
		trackFitPfoId.push_back(clusterID);
		trackFitX.push_back(trackState.GetPosition().GetX());
		trackFitY.push_back(trackState.GetPosition().GetY());
		trackFitZ.push_back(trackState.GetPosition().GetZ());
		trackFitQ.push_back(hitQ);
		trackFitRR.push_back(hitRR);
		trackFitdx.push_back(hitdx);
		float hitdQdx = hitdx > 0. ? hitQ/hitdx : -5.f;
		trackFitdQdx.push_back(hitdQdx);
	      }
	      // Rest of points
	      if ( idxPt > 0 ){
		const lar_content::LArTrackState& trackStatePrev = trackStateVector_out.at(idxPt-1);
		lengthSoFar+=std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackState.GetPosition() ) );
		float hitQ = trackState.GetCaloHit()->GetInputEnergy();
		float hitRR = trklength - lengthSoFar;
		float hitdx = 0.;
		// Middle points
		if ( idxPt < trackStateVector_out.size()-1 ) {
		  const lar_content::LArTrackState& trackStateNext = trackStateVector_out.at(idxPt+1);
		  hitdx = std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
		}
		// Last point
		else {
		  hitdx = std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackState.GetPosition() ) );
		}
		trackFitSliceId.push_back(sliceID);
		trackFitPfoId.push_back(clusterID);
		trackFitX.push_back(trackState.GetPosition().GetX());
		trackFitY.push_back(trackState.GetPosition().GetY());
		trackFitZ.push_back(trackState.GetPosition().GetZ());
		trackFitQ.push_back(hitQ);
		trackFitRR.push_back(hitRR);
		trackFitdx.push_back(hitdx);
		float hitdQdx = hitdx > 0. ? hitQ/hitdx : -5.f;
		trackFitdQdx.push_back(hitdQdx);
	      }
	    }
	  }
	}
      } // TRACK FIT

      if ( parameters.runShowerFit && (parameters.trackScoreCut < 0. || trackScore < parameters.trackScoreCut) ) {
	//std::cout << "I would have fit this as a shower..." << std::endl;
     
	//Begin Defining Shower Direction Through a PCA      
	CartesianVector centroid(0.f, 0.f, 0.f);
        lar_content::LArPcaHelper::EigenVectors eigenVecs;
        lar_content::LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        lar_content::LArPcaHelper::RunPca(caloHitList, centroid, eigenValues, eigenVecs);
	
	//Define directions to be positive
	const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

	//Define the shower start position
	//loop over the caloHitList 
	float dotProduct, projection;
	std::map<float, CartesianVector> ProjectionMap;

	CartesianVector startShower(0.f, 0.f, 0.f);
	//Find projections for each hit along the primary axis and save them into a map from least to greatest	
	

	for (const CaloHit *const pCaloHit3D : caloHitList){
		
		dotProduct = axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid);
		projection = dotProduct/axisDirection.GetMagnitude();
		startShower = pCaloHit3D->GetPositionVector();
		ProjectionMap.insert({projection, startShower});	
	}

	// constants for looping through projection
	//Define a proximity radius and proximity threshold
	//
	CartesianVector currentShower(0.f, 0.f, 0.f);		
	CartesianVector minProjectionVector(0.f, 0.f, 0.f);

	float currentX ,currentY, currentZ;
	
	float minProjectionValue;

	float HitProximityMinX, HitProximityMinY, HitProximityMinZ;
	float HitProximityMaxX, HitProximityMaxY, HitProximityMaxZ;

	int HitProximityRadius = 2;
	int ProximityHitsCounter = 0;
	int ProximityHitsThreshold = 15;

	//loop through minProjections then loop through hits
	//if proxthreshold is met, break the loop
	
	for (auto projection1 : ProjectionMap){
	
		minProjectionVector = projection1.second;
		minProjectionValue = projection1.first;
//	std::cout << "Projection Value: " << projection1.first << std::endl;	
		for (auto projection2 : ProjectionMap){

			currentShower = projection2.second;
			currentX = currentShower.GetX();
			currentY = currentShower.GetY();
			currentZ = currentShower.GetZ();
	
			HitProximityMinX = (minProjectionVector.GetX() - HitProximityRadius);
			HitProximityMaxX = (minProjectionVector.GetX() + HitProximityRadius);

			HitProximityMinY = (minProjectionVector.GetY() - HitProximityRadius);
                	HitProximityMaxY = (minProjectionVector.GetY() + HitProximityRadius);
		
			HitProximityMinZ = (minProjectionVector.GetZ() - HitProximityRadius);
                	HitProximityMaxZ = (minProjectionVector.GetZ() + HitProximityRadius);

			if(HitProximityMinX < currentX && currentX < HitProximityMaxX && HitProximityMinY < currentY && currentY < HitProximityMaxY && HitProximityMinZ < currentZ && currentZ < HitProximityMaxZ ){
			ProximityHitsCounter ++;
				}	
		}
        	if(ProximityHitsCounter >= ProximityHitsThreshold){
                	break;
                	}
	}
        	

	//Define dE/dx of the shower
	//start 
	float ShowerStartX, ShowerStartY, ShowerStartZ;

	int showerStartLength = 5;
	int showerStartWidth = 4;

	//std::vector<float> ShowerStartingHits;

	ShowerStartX = minProjectionVector.GetX();
	ShowerStartY = minProjectionVector.GetY();
	ShowerStartZ = minProjectionVector.GetZ();
	
	float distanceFromShowerStart;	
	//float deltaHitX, deltaHitY, deltaHitZ;	
	
	float hitPCAOpeningAngle, hitPositionAlongAxis, hitPositionFromAxis;
	float totalQ = 0;

	CartesianVector showerStartCurrentHit(0.f, 0.f, 0.f);
	CartesianVector showerStartPCAProjection(0.f, 0.f, 0.f);
	CaloHitList showerStartCaloHitList;	
	CartesianVector hitProjectedPosition(0.f, 0.f, 0.f);
	


	for(const CaloHit *const pShowerStartCaloHit3D : caloHitList){
		showerStartPCAProjection = centroid + (axisDirection * minProjectionValue);
		showerStartCurrentHit = pShowerStartCaloHit3D->GetPositionVector();
if(entryIdx == 31){
	std::cout <<"Event: "<< entryIdx << "current Hit" << showerStartCurrentHit << std::endl;
}
		totalQ += pShowerStartCaloHit3D->GetInputEnergy();	
		
		if(showerStartPCAProjection == showerStartCurrentHit){
			continue;
			}

		hitPCAOpeningAngle = axisDirection.GetOpeningAngle(showerStartPCAProjection - showerStartCurrentHit);		
	if(entryIdx == 31){
	std::cout << "hitPCAOpeningAngle: " << hitPCAOpeningAngle << std::endl;
}
		distanceFromShowerStart = std::sqrt(showerStartCurrentHit.GetDistanceSquared(showerStartPCAProjection));
	if(entryIdx == 31){
std::cout << "distanceFromShowerStart" << distanceFromShowerStart <<std::endl;	

}		hitPositionAlongAxis = distanceFromShowerStart * std::abs(std::cos(hitPCAOpeningAngle));
		hitPositionFromAxis = distanceFromShowerStart * std::sin(hitPCAOpeningAngle);		
	if(entryIdx == 31){
	std::cout << "Along Axis: " << hitPositionAlongAxis << "From Axis: " << hitPositionFromAxis << std::endl;
}
		if(hitPositionAlongAxis < showerStartLength && hitPositionFromAxis < showerStartWidth){
if(entryIdx == 31){	
 std::cout<< "Event: " << entryIdx << "Along Axis: " << hitPositionAlongAxis << "From Axis: " << hitPositionFromAxis << std::endl;
}
			showerStartCaloHitList.push_back(pShowerStartCaloHit3D);
               

	       //Add positions of hits to a branch to look at later
			shwrStartPointsX.push_back(showerStartCurrentHit.GetX());
                        shwrStartPointsY.push_back(showerStartCurrentHit.GetY());
                        shwrStartPointsZ.push_back(showerStartCurrentHit.GetZ());
			} 
		else{
			continue;
			}	

	}	

   	shwrTotalE.push_back(totalQ);

      // std::cout << "Number of Shower Start Points: " << showerStartCaloHitList.size() << std::endl;


	if(showerStartCaloHitList.size() < 2){
		std::cout<< "Reco Hit List Too Small. Skipping fit. Event is: " << entryIdx << std::endl;
        
		startTrkSlidingFitDirX.push_back(-9999.);
        	startTrkSlidingFitDirY.push_back(-9999.);
        	startTrkSlidingFitDirZ.push_back(-9999.);
	
		startTrkPCADirX.push_back(-9999.);
		startTrkPCADirY.push_back(-9999.);
		startTrkPCADirZ.push_back(-9999.);
	
		continue;
		}

	else{
//track fit of beginning points
	float showerSlidingFitHalfWindow = 20;
	std::vector<int> showerFitIndexVector;

	lar_content::LArTrackStateVector showerFitTrackStateVector;
	
	lar_content::LArPfoHelper::GetSlidingFitTrajectory( &showerStartCaloHitList, vertexVector, showerSlidingFitHalfWindow, parameters.pixelPitch, showerFitTrackStateVector, &showerFitIndexVector, true );

if(showerFitTrackStateVector.size() < 1){

std::cout << "TRACK STATE FAILED " << " Event is: " << entryIdx << "Slice: "<< sliceID << "Cluster: " << clusterID << std::endl;

        startTrkSlidingFitDirX.push_back(-9999.);
        startTrkSlidingFitDirY.push_back(-9999.);
        startTrkSlidingFitDirZ.push_back(-9999.);
	
	shwrdEdx.push_back(-9999.);
	
	}

else{

	startTrkSlidingFitDirX.push_back(showerFitTrackStateVector.front().GetDirection().GetX());
	startTrkSlidingFitDirY.push_back(showerFitTrackStateVector.front().GetDirection().GetY());
	startTrkSlidingFitDirZ.push_back(showerFitTrackStateVector.front().GetDirection().GetZ());


	//dQ/dx
	double hitQ_Shower;
	float cosgamma, projectedShowerDir, pitch;	
	std::vector<double> Q_ShowerVector;

	for(unsigned int iTrackHit = 0; iTrackHit < showerFitTrackStateVector.size(); ++iTrackHit){
		const lar_content::LArTrackState& showerFitTrackState = showerFitTrackStateVector.at(iTrackHit);
		hitQ_Shower = showerFitTrackState.GetCaloHit()->GetInputEnergy();
		Q_ShowerVector.push_back(hitQ_Shower);
		}

	projectedShowerDir = std::sqrt(std::pow(showerFitTrackStateVector.front().GetDirection().GetY(), 2) + std::pow(showerFitTrackStateVector.front().GetDirection().GetZ(), 2));

if(projectedShowerDir == 0){
	std::cout << "Pitch is calculated to be inf. dEdx calculation will be skipped." << std::endl;
	shwrdEdx.push_back(-9999.);
	}
else{
	cosgamma = showerFitTrackStateVector.front().GetDirection().GetZ() / projectedShowerDir;

		if(cosgamma == 0){
			std::cout << "Pitch is calculated to be inf. dEdx calculation will be skipped." << std::endl;
			shwrdEdx.push_back(-9999.);
		}	
		else{
			pitch = parameters.pixelPitch/cosgamma;

			if(pitch == 0){
				std::cout << "Pitch is calculated to be 0. dEdx calculation will be skipped." << std::endl;
				shwrdEdx.push_back(-9999.);
			}
			else{
				double dQdx = TMath::Median(Q_ShowerVector.size(), &Q_ShowerVector[0]) / pitch;
				shwrdEdx.push_back(dQdx);
				}	
		}

	}
}	
//PCA of beginning points
//

        //Begin Defining Shower Direction Through a PCA
 CartesianVector centroidShowerStart(0.f, 0.f, 0.f);
 lar_content::LArPcaHelper::EigenVectors eigenVecsShowerStart;
 lar_content::LArPcaHelper::EigenValues eigenValuesShowerStart(0.f, 0.f, 0.f);
 lar_content::LArPcaHelper::RunPca(showerStartCaloHitList, centroidShowerStart, eigenValuesShowerStart, eigenVecsShowerStart);

//Define directions to be positive
 const CartesianVector axisDirectionShowerStart(eigenVecsShowerStart.at(0).GetZ() > 0.f ? eigenVecsShowerStart.at(0) : eigenVecsShowerStart.at(0) * -1.f);
 
 	startTrkPCADirX.push_back(axisDirectionShowerStart.GetX());
        startTrkPCADirY.push_back(axisDirectionShowerStart.GetY());
        startTrkPCADirZ.push_back(axisDirectionShowerStart.GetZ());


}

	//Fill shower branches
	shwrLen.push_back(minProjectionValue); //doing this to silence a warning
        shwrCentroidX.push_back(centroid.GetX());
        shwrCentroidY.push_back(centroid.GetY());
        shwrCentroidZ.push_back(centroid.GetZ());
        shwrStartX.push_back(ShowerStartX);
        shwrStartY.push_back(ShowerStartY);
        shwrStartZ.push_back(ShowerStartZ);
        shwrSliceId.push_back(sliceID);
	shwrClusterId.push_back(clusterID);
        shwrDirX.push_back(axisDirection.GetX());
        shwrDirY.push_back(axisDirection.GetY());
        shwrDirZ.push_back(axisDirection.GetZ());
		


	 } // SHOWER FIT


}

    if ( parameters.runTrackFit ) {
      fOut.FillTrackBranches(trkStartX,trkStartY,trkStartZ,trkStartDirX,trkStartDirY,trkStartDirZ,trkEndX,trkEndY,trkEndZ,trkEndDirX,trkEndDirY,trkEndDirZ,trkLen);
      fOut.FillTrackCaloBranches(trackFitSliceId,trackFitPfoId,trackFitX,trackFitY,trackFitZ,trackFitQ,trackFitRR,trackFitdx,trackFitdQdx);
    }
   if ( parameters.runShowerFit){
     fOut.FillShowerBranches(shwrCentroidX,shwrCentroidY,shwrCentroidZ,shwrStartX, shwrStartY, shwrStartZ, shwrDirX, shwrDirY, shwrDirZ, shwrLen, shwrSliceId, shwrClusterId, startTrkSlidingFitDirX, startTrkSlidingFitDirY, startTrkSlidingFitDirZ, startTrkPCADirX, startTrkPCADirY, startTrkPCADirZ, shwrStartPointsX, shwrStartPointsY, shwrStartPointsZ, shwrdEdx, shwrTotalE);
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

    bool hasInputFile=false;
    bool hasXmlFile=false;

    while ((cOpt = getopt(argc, argv, "x:f:o:h")) != -1)
    {
      switch (cOpt)
      {
        case 'x':
	  parameters.xmlName = optarg;
	  hasXmlFile = true;
	  break;
        case 'f':
	  parameters.fileName = optarg;
	  hasInputFile = true;
	  break;
        case 'o':
	  parameters.outfileName = optarg;
	  break;
        case 'h':
        default:
	  return PrintOptions();
      }
    }

    bool passed=hasXmlFile && hasInputFile;
    if(!passed)
    {
        return PrintOptions();
    }
    return passed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ReadSettings(ParameterStruct &parameters)
{
  TiXmlDocument xmlDocument(parameters.xmlName);

  if (!xmlDocument.LoadFile()) {
    std::cout << "XML document (" << parameters.xmlName << ") not loaded. Returning." << std::endl;
    return false;
  }

  const TiXmlHandle xmlDocumentHandle(&xmlDocument);
  const TiXmlHandle xmlHandle(TiXmlHandle(xmlDocumentHandle.FirstChildElement().Element()));

  try
  {
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldRunTrackFit", parameters.runTrackFit) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldRunShowerFit", parameters.runShowerFit) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackScoreCut", parameters.trackScoreCut) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PixelPitch", parameters.pixelPitch) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldApplyHitThreshold", parameters.applyThreshold) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChargeThreshold", parameters.thresholdVal) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldVoxelizeZ", parameters.voxelizeZ) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VoxelHalfWidthZ", parameters.voxelZHW) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseVoxelizedStartStop", parameters.useVoxelizedStartStop) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbosity", parameters.verbosity) );
  }
  catch (StatusCodeException &statusCodeException)
  {
    std::cout << "Failed to initialized parameters in the XML file. Status code " << statusCodeException.GetStatusCode() << ". Returning." << std::endl;
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl
              << "./bin/PandoraOuterface -f [file name] -o [out name] -p [pixel pitch] -c [track score cut] -t -s" << std::endl
	      << "    -t = run track fit" << std::endl
	      << "    -s = run shower fit"
              << std::endl;

    return false;
}

} // namespace lar_nd_postreco
