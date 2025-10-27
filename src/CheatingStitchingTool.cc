/**
 *  @file   src/CheatingStitching.cc
 *
 *  @brief  Implementation of the cheating stitching tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <larpandoracontent/LArHelpers/LArMCParticleHelper.h>
#include <larpandoracontent/LArHelpers/LArPfoHelper.h>

#include "CheatingStitchingTool.h"

using namespace pandora;

namespace lar_content
{
CheatingStitchingTool::CheatingStitchingTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingStitchingTool::Run(const MasterAlgorithm *const pAlgorithm, const PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap,
    PfoToFloatMap &stitchedPfosToX0Map)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::cout << "Cheating stitching tool: there are " << pMultiPfoList->size() << " PFOs to consider for stitching" << std::endl;
    std::cout << "There are " << pfoToLArTPCMap.size() << " PFOs with associated TPCs" << std::endl;

    if (this->GetPandora().GetGeometry()->GetLArTPCMap().size() < 2)
        return;

    if (pfoToLArTPCMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    std::cout << "Now performing cheating stitching based on MCParticle associations" << std::endl;

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));
    std::cout << "There are " << pMCParticleList->size() << " MCParticles in the event" << std::endl;
    const PfoList *pPfoList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pPfoList));
    std::cout << "There are " << pPfoList->size() << " PFOs in the event" << std::endl;

    std::map<const MCParticle *, std::map<unsigned int, std::vector<const ParticleFlowObject *>>> mcParticleToLArTPCToPfosMap;

    // Find all PFOs associated to each MCParticle
    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        CaloHitList allHits;
        LArPfoHelper::GetAllCaloHits(pPfo, allHits);
        std::map<const MCParticle *, unsigned int> mcParticleToHitCountMap;

        for (const CaloHit *const pCaloHit : allHits)
        {
            MCParticleVector mcParticleVector;
            for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap())
                mcParticleVector.push_back(weightMapEntry.first);
            std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

            // Consider only the primary MCParticle associated to this hit
            if (!mcParticleVector.empty())
                ++mcParticleToHitCountMap[mcParticleVector.front()];
        }

        // Find the MCParticle with the most hits associated to this PFO
        const MCParticle *pMCParticle{nullptr};
        size_t maxHits(0);

        for (const auto &mcParticleHitCountPair : mcParticleToHitCountMap)
        {
            if (mcParticleHitCountPair.second > maxHits)
            {
                maxHits = mcParticleHitCountPair.second;
                pMCParticle = mcParticleHitCountPair.first;
            }
        }

        // Add this PFO to the map for this MCParticle
        if (nullptr == pMCParticle)
            continue;

        const unsigned int tpcVolumeId = pfoToLArTPCMap.at(pPfo)->GetLArTPCVolumeId();
        mcParticleToLArTPCToPfosMap[pMCParticle][tpcVolumeId].push_back(pPfo);
    }

    std::cout << "There are " << mcParticleToLArTPCToPfosMap.size() << " MCParticles associated to PFOs" << std::endl;

    // Now, join together any PFOs from different TPCs associated to the same MCParticle
    for (const auto &mcParticlePfosPair : mcParticleToLArTPCToPfosMap)
    {
        const std::map<unsigned int, std::vector<const ParticleFlowObject *>> &tpcToPfosMap = mcParticlePfosPair.second;

        // Nothing to stitch if only one TPC involved
        if (tpcToPfosMap.size() < 2)
            continue;

        std::cout << "Stitching together PFOs from " << tpcToPfosMap.size() << " TPCs" << std::endl;

        // Sort the PFOs in each TPC by number of hits.
        // Since this is JUST cheated stitching, we will take the largest PFO in each TPC as the "main" one to stitch to.
        std::map<unsigned int, const ParticleFlowObject *> tpcToMainPfoMap;
        for (const auto &tpcPfosPair : tpcToPfosMap)
        {
            const unsigned int tpcVolumeId = tpcPfosPair.first;
            const std::vector<const ParticleFlowObject *> &pfosInTpc = tpcPfosPair.second;
            const ParticleFlowObject *pLargestPfo{nullptr};

            size_t maxHits(0);
            for (const ParticleFlowObject *const pPfo : pfosInTpc)
            {
                const size_t nHits = LArPfoHelper::GetNumberOfThreeDHits(pPfo) + LArPfoHelper::GetNumberOfTwoDHits(pPfo);
                if (nHits > maxHits)
                {
                    maxHits = nHits;
                    pLargestPfo = pPfo;
                }
            }

            tpcToMainPfoMap[tpcVolumeId] = pLargestPfo;
        }

        std::cout << "Stitching " << tpcToMainPfoMap.size() << " PFOs together" << std::endl;

        // Stitch all other PFOs in this TPC to the main one
        auto tpcToMainPfoIter = tpcToMainPfoMap.begin();
        const ParticleFlowObject *const pReferencePfo = tpcToMainPfoIter->second;
        ++tpcToMainPfoIter;

        std::cout << "Reference PFO has " << LArPfoHelper::GetNumberOfThreeDHits(pReferencePfo) << " 3D hits and "
                  << LArPfoHelper::GetNumberOfTwoDHits(pReferencePfo) << " 2D hits" << std::endl;

        for (; tpcToMainPfoIter != tpcToMainPfoMap.end(); ++tpcToMainPfoIter)
        {
            const ParticleFlowObject *const pCurrentPfo = tpcToMainPfoIter->second;

            // Stitch the current PFO to the reference PFO
            pAlgorithm->StitchPfos(pReferencePfo, pCurrentPfo, pfoToLArTPCMap);
        }

        std::cout << "Stitched PFO now has " << LArPfoHelper::GetNumberOfThreeDHits(pReferencePfo) << " 3D hits and "
                  << LArPfoHelper::GetNumberOfTwoDHits(pReferencePfo) << " 2D hits" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingStitchingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
