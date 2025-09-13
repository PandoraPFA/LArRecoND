/**
 *  @file   src/CheatingEventSlicingThreeDTool.cc
 *
 *  @brief  Implementation of the 3D event slicing tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
//new
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "Helpers/MCParticleHelper.h"

/*#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"*/

#include "CheatingEventSlicingThreeDTool.h"

//#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

CheatingEventSlicingThreeDTool::CheatingEventSlicingThreeDTool(){}
	
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingEventSlicingThreeDTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//New
////------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingThreeDTool::RunSlicing(const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    const HitTypeToNameMap & /*clusterListNames*/, Slice3DList &slice3DList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    MCParticleToSlice3DMap mcParticleToSliceMap;
    std::cout <<"debug00" << std::endl;
    this->InitializeMCParticleToSlice3DMap(pAlgorithm, caloHitListNames, mcParticleToSliceMap);
    std::cout <<"debug0" << std::endl;
    this->FillSlices(pAlgorithm, TPC_VIEW_U, caloHitListNames, mcParticleToSliceMap);
    std::cout <<"debug1" << std::endl;
    this->FillSlices(pAlgorithm, TPC_VIEW_V, caloHitListNames, mcParticleToSliceMap);
    std::cout <<"debug2" << std::endl;
    this->FillSlices(pAlgorithm, TPC_VIEW_W, caloHitListNames, mcParticleToSliceMap);
    std::cout <<"debug3" << std::endl;
    this->FillSlices(pAlgorithm, TPC_3D, caloHitListNames, mcParticleToSliceMap);
    std::cout <<"debug4" << std::endl;

    MCParticleVector mcParticleVector;
    for (const auto &mapEntry : mcParticleToSliceMap)
        mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    std::cout <<"debug6" << std::endl;
    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const Slice3D &slice(mcParticleToSliceMap.at(pMCParticle));

        if (!slice.m_caloHitListU.empty() || !slice.m_caloHitListV.empty() || !slice.m_caloHitListW.empty() || !slice.m_caloHitList3D.empty())
            slice3DList.push_back(slice);
    }
    std::cout <<"debug7" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingThreeDTool::InitializeMCParticleToSlice3DMap(
    const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames, MCParticleToSlice3DMap &mcParticleToSliceMap) const
{
    for (const auto &mapEntry : caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, mapEntry.second, pCaloHitList));
        std::cout << "mapEntry.second = " << mapEntry.second << std::endl;

        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            MCParticleVector mcParticleVector;
	    std::cout << "Cheated algo weight map size = " << pCaloHit->GetMCParticleWeightMap().size() << std::endl;
            for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap())
            {
                std::cout << "entry n. = " << mcParticleVector.size() << " weightMapEntry.second = " << weightMapEntry.second << std::endl;
                mcParticleVector.push_back(weightMapEntry.first);
            }
            std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);
            

	    const MCParticle* mainMCParticle(nullptr);
            try
            {
                mainMCParticle=MCParticleHelper::GetMainMCParticle(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
            }
	    if(!mainMCParticle || mcParticleToSliceMap.count(LArMCParticleHelper::GetParentMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit)))))
            {
		continue;
            }

            //if (mcParticleToSliceMap.count(LArMCParticleHelper::GetParentMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit)))))
            //{ 
            //    continue;
	    //}

            if (!mcParticleToSliceMap.insert(MCParticleToSlice3DMap::value_type(LArMCParticleHelper::GetParentMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))), Slice3D())).second)
	    {
                throw StatusCodeException(STATUS_CODE_FAILURE);
	    }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingThreeDTool::FillSlices(const Algorithm *const pAlgorithm, const HitType hitType,
    const HitTypeToNameMap &caloHitListNames, MCParticleToSlice3DMap &mcParticleToSliceMap) const
{
    if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType) && (TPC_3D != hitType))
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, caloHitListNames.at(hitType), pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            MCParticleToSlice3DMap::iterator mapIter = mcParticleToSliceMap.find(LArMCParticleHelper::GetParentMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))));

            if (mcParticleToSliceMap.end() == mapIter)
	    {
                throw StatusCodeException(STATUS_CODE_FAILURE);
	    }
            Slice3D &slice(mapIter->second);
            CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU
                    : (TPC_VIEW_V == hitType)                ? slice.m_caloHitListV
                    : (TPC_VIEW_W == hitType)                ? slice.m_caloHitListW
                                                             : slice.m_caloHitList3D);
            caloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------


} // namespace lar_content
