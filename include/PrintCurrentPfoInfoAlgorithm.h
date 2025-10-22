/**
 *  @file   include/PrintCurrentPfoInfoAlgorithm.h
 *
 *  @brief  Header file for the 3D pre processing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PRINT_CURRENT_PFO_INFO_ALGORITHM_H 
#define LAR_PRINT_CURRENT_PFO_INFO_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{
  //------------------------------------------------------------------------------------------------------------------------------------------

  /**
 *  @brief  PrintCurrentPfoInfoAlgorithm class
 */
  class PrintCurrentPfoInfoAlgorithm : public pandora::Algorithm
  {
    public:
      /*
       * @brief Default constructor
       * */
      PrintCurrentPfoInfoAlgorithm();
      static std::map<std::string, int> AlgoExecutionCount;
      static std::ofstream pfoInfoOutputFile;
    
    private:

      pandora::StatusCode Run();
      void PrintCaloHitsInfo(const pandora::CaloHitList *& pCaloHitList, std::string HitsName, std::string STAGE);
      void PrintClusterListInfo(const pandora::ClusterList *& pClusterList, std::string clusterName, std::string STAGE);
      void PrintPfoInfo(const pandora::ParticleFlowObject *& pPfo, std::string STAGE, std::string LIST_NAME);
      pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

      std::string m_inputStageName;                 ///< Name of the stage at which the printing is called
      std::vector<std::string> m_inputPfoListNames; ///< Name of the input pfo list(s)
      std::vector<std::string> m_inputVertexListNames; 
      std::string m_inputCaloHitList3DName;         ///< Name of the input 3D calo hit list
      // hits
      std::string m_inputCaloHitListUName;
      std::string m_inputCaloHitListVName;
      std::string m_inputCaloHitListWName;
      std::string m_inputCaloHitList2DName;
      // clusters
      std::string m_inputClusterListName3D;
      std::string m_inputClusterListNameU;
      std::string m_inputClusterListNameV;
      std::string m_inputClusterListNameW;
  };                                    

} // namespace lar_content

#endif // LAR_PRINT_CURRENT_PFO_INFO_ALGORITHM_H 

