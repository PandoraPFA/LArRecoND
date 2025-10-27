/**
 *  @file   include/CheatingStitching.h
 *
 *  @brief  Header file for the cheating stitching tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_STITCHING_TOOL_H
#define LAR_CHEATING_STITCHING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/StitchingBaseTool.h"

namespace lar_content
{

/**
 *  @brief  CheatingStitchingTool class
 */
class CheatingStitchingTool : public StitchingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingStitchingTool();

    void Run(const MasterAlgorithm *const pAlgorithm, const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap,
        PfoToFloatMap &stitchedPfosToX0Map);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_STITCHING_TOOL_H
