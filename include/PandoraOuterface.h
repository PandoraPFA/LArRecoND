/**
 *  @file   LArRecoND/include/PandoraOuterface.h
 *
 *  @brief  Header file for PandoraOuterface.
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_OUTERFACE_H
#define PANDORA_ND_OUTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"

#ifdef USE_EDEPSIM
#include "TG4Event.h"
#endif

#include "TGeoManager.h"
#include "TGeoNode.h"

#include "LArGrid.h"
#include "LArHitInfo.h"
#include "LArRecoNDFormat.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

namespace pandora
{
class Pandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_postreco
{

struct ParameterStruct
{
  float pixelPitch = 0.4;
  float trackScoreCut = 0.5;
  std::string fileName = "";
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Read in the file and process the post-reconstruction
 *
 *  @param  parameters the input parameters controlling aspects of post-reco
 *
 */
void ProcessPostReco(const ParameterStruct &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 *
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], ParameterStruct &parameters);

/**
 *  @brief  Print the list of configurable options
 *
 *  @return false, to force abort
 */
bool PrintOptions();

} // namespace lar_nd_postreco

#endif // #ifndef PANDORA_ND_OUTERFACE_H
