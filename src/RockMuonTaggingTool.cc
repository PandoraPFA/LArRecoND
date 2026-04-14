/**
 *  @file   src/RockMuonTaggingTool.cc
 *
 *  @brief  Implementation of the rock muon tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingTool.h"

#include "RockMuonTaggingTool.h"

using namespace pandora;

namespace lar_content
{

RockMuonTaggingTool::RockMuonTaggingTool() :
    m_tagRockMuons(false),
    m_marginX(5.f), // [cm]
    m_marginY(5.f),
    m_marginZ(5.f)
  {
  }

//------------------------------------------------------------------------------------------------------------------------------------------
bool RockMuonTaggingTool::IsOutsideBox(const float x, const float y, const float z) const 
{
  const float BoxXmin = m_face_Xa + m_marginX;
  const float BoxXmax = m_face_Xc - m_marginX;

  const float BoxYmin = m_face_Yb + m_marginY;
  const float BoxYmax = m_face_Yt - m_marginY;

  const float BoxZmin = m_face_Zu + m_marginZ;
  const float BoxZmax = m_face_Zd - m_marginZ;

  bool IsOutRangeX = (x < BoxXmin || x > BoxXmax);
  bool IsOutRangeY = (y < BoxYmin || y > BoxYmax);
  bool IsOutRangeZ = (z < BoxZmin || z > BoxZmax);

  return (IsOutRangeZ || IsOutRangeY || IsOutRangeX);
}

//------------------------------------------------------------------------------------------------------------------------------------------
void RockMuonTaggingTool::CheckIfThroughgoing(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsThroughgoingMap) const
{
    for (const CRCandidate &candidate : candidates)
    {
      bool isThroughgoing = (
          (candidate.m_endPoint1.GetX() != std::numeric_limits<float>::max()) && // for some reason some candidates happen to have "default" values set as inf 
          IsOutsideBox(candidate.m_endPoint1.GetX(), candidate.m_endPoint1.GetY(), candidate.m_endPoint1.GetZ()) &&
          IsOutsideBox(candidate.m_endPoint2.GetX(), candidate.m_endPoint2.GetY(), candidate.m_endPoint2.GetZ()));

      if (!pfoToIsThroughgoingMap.insert(PfoToBoolMap::value_type(candidate.m_pPfo, isThroughgoing)).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
} 

//------------------------------------------------------------------------------------------------------------------------------------------
void RockMuonTaggingTool::TagRockMuons(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsLikelyRockMuonMap, const PfoToBoolMap &pfoToIsThroughgoingMap) const
{
    int nof_tagged_rockmus = 0;
    for (const CRCandidate &candidate : candidates)
    {
	bool likelyRockMuon = false;

        if(m_tagRockMuons && (pfoToIsThroughgoingMap.at(candidate.m_pPfo)))
        {
          likelyRockMuon = true;
          nof_tagged_rockmus++;
        }

	if (!pfoToIsLikelyRockMuonMap.insert(PfoToBoolMap::value_type(candidate.m_pPfo, likelyRockMuon)).second)
          throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
    std::cout << "nof tagged rock muons : " << nof_tagged_rockmus << "\n";
}

//------------------------------------------------------------------------------------------------------------------------------------------
void RockMuonTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const /*pAlgorithm*/)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // TODO First time only, TODO Refactor with master algorithm
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    m_face_Xa = parentMinX;
    m_face_Xc = parentMaxX;
    m_face_Yb = parentMinY;
    m_face_Yt = parentMaxY;
    m_face_Zu = parentMinZ;
    m_face_Zd = parentMaxZ;
    
    if(m_tagRockMuons)
     std::cout << "Detector boundaries to tag rock muons : x [" 
               << m_face_Xa + m_marginX <<"," << m_face_Xc - m_marginX<< "] "
               << "y [" << m_face_Yb + m_marginY <<"," << m_face_Yt - m_marginY<< "] "
               << "z [" << m_face_Zu + m_marginZ <<"," << m_face_Zd - m_marginZ << "] \n";

    PfoToPfoListMap pfoAssociationMap;
    this->GetPfoAssociations(parentCosmicRayPfos, pfoAssociationMap);

    PfoToSliceIdMap pfoToSliceIdMap;
    this->SliceEvent(parentCosmicRayPfos, pfoAssociationMap, pfoToSliceIdMap);

    CRCandidateList candidates;
    this->GetCRCandidates(parentCosmicRayPfos, pfoToSliceIdMap, candidates);

    PfoToBoolMap pfoToInTimeMap;
    this->CheckIfInTime(candidates, pfoToInTimeMap);

    PfoToBoolMap pfoToIsContainedMap;
    this->CheckIfContained(candidates, pfoToIsContainedMap);

    PfoToBoolMap pfoToIsTopToBottomMap;
    this->CheckIfTopToBottom(candidates, pfoToIsTopToBottomMap);

    UIntSet neutrinoSliceSet;
    this->GetNeutrinoSlices(candidates, pfoToInTimeMap, pfoToIsContainedMap, neutrinoSliceSet);

    PfoToBoolMap pfoToIsLikelyCRMuonMap;
    this->TagCRMuons(candidates, pfoToInTimeMap, pfoToIsTopToBottomMap, neutrinoSliceSet, pfoToIsLikelyCRMuonMap);

    // (clear) rock muon tagging
    PfoToBoolMap pfoToIsLikelyRockMuonMap;
    if(m_tagRockMuons)
    {
      PfoToBoolMap pfoToIsThroughgoingMap;
      this->CheckIfThroughgoing(candidates, pfoToIsThroughgoingMap);

      this->TagRockMuons(candidates, pfoToIsLikelyRockMuonMap, pfoToIsThroughgoingMap);
    }
    // end rock muon tagging

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        const bool is_rock = pfoToIsLikelyRockMuonMap.at(pPfo); // dummy-proof (for me, basically) 
        const bool is_cosmic = pfoToIsLikelyCRMuonMap.at(pPfo); 

        if(m_tagRockMuons)
        {
          if (!(is_rock || is_cosmic))
            ambiguousPfos.push_back(pPfo);
        }
        else
        {
          if (!is_cosmic)
            ambiguousPfos.push_back(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode RockMuonTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CutMode", m_cutMode));
    std::transform(m_cutMode.begin(), m_cutMode.end(), m_cutMode.begin(), ::tolower);

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AngularUncertainty", m_angularUncertainty));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PositionalUncertainty", m_positionalUncertainty));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAssociationDist", m_maxAssociationDist));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitThreshold", m_minimumHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMargin", m_inTimeMargin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMaxX0", m_inTimeMaxX0));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TagRockMuons", m_tagRockMuons));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MarginX", m_marginX));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MarginY", m_marginY));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MarginZ", m_marginZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxNeutrinoCosTheta", m_maxNeutrinoCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosmicCosTheta", m_minCosmicCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCosmicCurvature", m_maxCosmicCurvature));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content
