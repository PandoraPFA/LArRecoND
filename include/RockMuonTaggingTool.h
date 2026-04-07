/**
 *  @file   include/RockMuonTaggingTool.h
 *
 *  @brief  Header file for the rock muon tagging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ROCK_MUON_TAGGING_TOOL_H
#define LAR_ROCK_MUON_TAGGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/CosmicRayTaggingTool.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  RockMuonTaggingAlgorithm class
 */
class RockMuonTaggingTool : public CosmicRayTaggingTool
{
  public:
    /**
     *  @brief  Default constructor
     */
    RockMuonTaggingTool();

    void FindAmbiguousPfos( const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm) override;

  private:
    /**
     *  @brief Check if a 3D point is inside the detector boundaies with margins 
     *
     *  @param  x point x coordinate
     *  @param  y point y coordinate
     *  @param  z point z coordinate
     */   
    bool IsOutsideBox(const float x, const float y, const float z) const;
  
    /**
     *  @brief  Check if each candidate is throughgoing (i.e emerging and exiting from any of the detector boundaies)
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsThroughgoingMap output mapping between candidates Pfos and if they are top to bottom
     */
    void CheckIfThroughgoing(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsThroughgoingMap) const;

    /**
     *  @brief  Tag Pfos which are likely to be a CR muon
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsLikelyRockMuonMap to receive the output mapping between Pfos and a boolean deciding if they are likely a rock muon
     *  @param  pfoToIsThroughgoingMap output mapping between candidates Pfos and if they are top to bottom
     */
    void TagRockMuons(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsLikelyRockMuonMap, const PfoToBoolMap &pfoToIsThroughgoingMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    bool m_tagRockMuons;  ///< bool to activate tagging of rock muons
    float m_marginX; ///< the minimum distance from the dector X-face to define a fiducial volume for tagging
    float m_marginY; ///< the minimum distance from the dector Y-face to define a fiducial volume for tagging
    float m_marginZ; ///< the minimum distance from the dector Z-face to define a fiducial volume for tagging
};
} // namespace lar_content

#endif // #ifndef LAR_ROCK_MUON_TAGGING_TOOL_H
