/**
 *  @file   src/PrintCurrentPfoInfoAlgorithm.cc
 *
 *  @brief  Implementation of the 3D list preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "PreProcessingThreeDAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include <fstream>
#include <iterator>

#include "PrintCurrentPfoInfoAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PrintCurrentPfoInfoAlgorithm::PrintCurrentPfoInfoAlgorithm() : 
  m_inputStageName{""},
  // hits info
  m_inputCaloHitList3DName{""},
  m_inputCaloHitListUName{""},
  m_inputCaloHitListVName{""},
  m_inputCaloHitListWName{""},
  m_inputCaloHitList2DName{""},
  // clusters info
  m_inputClusterListName3D{""},
  m_inputClusterListNameU{""},
  m_inputClusterListNameV{""},
  m_inputClusterListNameW{""}
{
}

  //------------------------------------------------------------------------------------------------------------------------------------------
std::ofstream PrintCurrentPfoInfoAlgorithm::pfoInfoOutputFile("pfos_info.txt");

std::map<std::string, int> PrintCurrentPfoInfoAlgorithm::AlgoExecutionCount;

StatusCode PrintCurrentPfoInfoAlgorithm::Run()
{

  std::cout << "Running PrintCurrentPfoInfoAlgorithm after STAGE : "<< m_inputStageName << "\n";

  PrintCurrentPfoInfoAlgorithm::AlgoExecutionCount[m_inputStageName]++;

  // print hits info ----------------------------------------------------------
  const CaloHitList *pCaloHitList3D{nullptr};
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputCaloHitList3DName, pCaloHitList3D));
  PrintCaloHitsInfo(pCaloHitList3D, m_inputCaloHitList3DName, m_inputStageName);

  const CaloHitList *pCaloHitListU{nullptr};
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListUName, pCaloHitListU));
  PrintCaloHitsInfo(pCaloHitListU, m_inputCaloHitListUName, m_inputStageName);

  const CaloHitList *pCaloHitListW{nullptr};
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListWName, pCaloHitListW));
  PrintCaloHitsInfo(pCaloHitListW, m_inputCaloHitListWName, m_inputStageName);

  const CaloHitList *pCaloHitListV{nullptr};
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListVName, pCaloHitListV));
  PrintCaloHitsInfo(pCaloHitListV, m_inputCaloHitListVName, m_inputStageName);

  const CaloHitList *pCaloHitList2D{nullptr};
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputCaloHitList2DName, pCaloHitList2D));
  PrintCaloHitsInfo(pCaloHitList2D, m_inputCaloHitList2DName, m_inputStageName);

  // print clusters info -----------------------------------------------------------
  const ClusterList *pClusterList3D{nullptr};

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListName3D, pClusterList3D));
  PrintClusterListInfo(pClusterList3D, m_inputClusterListName3D, m_inputStageName);

  const ClusterList *pClusterListU{nullptr};

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS,STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameU, pClusterListU));
  PrintClusterListInfo(pClusterListU, m_inputClusterListNameU, m_inputStageName);

  const ClusterList *pClusterListV{nullptr};

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameV, pClusterListV));
  PrintClusterListInfo(pClusterListV, m_inputClusterListNameV, m_inputStageName);

  const ClusterList *pClusterListW{nullptr};

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameW, pClusterListW));
  PrintClusterListInfo(pClusterListW, m_inputClusterListNameW, m_inputStageName);

  // print pfo info --------------------------------------------------------------------
  for (unsigned int i = 0; i < m_inputPfoListNames.size(); ++i)
  {
    const std::string pfoListName(m_inputPfoListNames.at(i));
    const PfoList *pPfoList{nullptr};

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
      continue;

    for (const ParticleFlowObject *pPfo : (*pPfoList))
      PrintPfoInfo(pPfo, m_inputStageName, pfoListName);

   }

  // pfoInfoOutputFile.close();
  
  // print CR Vertices info --------------------------------------------------------------------
  for (unsigned int i = 0; i < m_inputVertexListNames.size(); ++i)
  {
    const std::string VertexListName(m_inputVertexListNames.at(i));
    const VertexList *pVertexList{nullptr};
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, VertexListName, pVertexList));
    
    if (!pVertexList || pVertexList->empty())
      continue;

    for (const Vertex *vertex : (*pVertexList))
    {
      std::cout << "{ \"label\" : " << vertex->GetVertexLabel() 
                << ", \"type\" : " << vertex->GetVertexType()
                << ", \"X0\" : " << vertex->GetX0()
                << ", \"x\" : " << vertex->GetPosition().GetX() 
                << ", \"y\" : " << vertex->GetPosition().GetY() 
                << ", \"z\" : " << vertex->GetPosition().GetZ() 
                << "}\n";
    }

  }

  return STATUS_CODE_SUCCESS;
}

void PrintCurrentPfoInfoAlgorithm::PrintPfoInfo(const ParticleFlowObject *& pPfo, std::string STAGE, std::string LIST_NAME)
{
      ClusterList cluster3DList, cluster2DList;
      LArPfoHelper::GetThreeDClusterList(pPfo, cluster3DList);
      LArPfoHelper::GetTwoDClusterList(pPfo, cluster2DList);
    
      if(!pfoInfoOutputFile.is_open())
      {
        std::cout << "Warning: pfoInfoOutputFile not open \n";
        return;
      }

      for (const Cluster *pCluster : cluster3DList)
      {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *pCaloHit : caloHitList)
        {
          pfoInfoOutputFile   << "{ \"STAGE\" : " << "\"" << STAGE << "\""
                    << ", \"CallNumber\" : " << "\"" <<AlgoExecutionCount[m_inputStageName] << "\""
                    << ", \"pfoListName\" : " << "\"" << LIST_NAME << "\""
                    << ", \"pfo\" : " << "\"" << pPfo + AlgoExecutionCount[m_inputStageName] << "\""
                    << ", \"NofCluters\" : " << "\"" << cluster3DList.size()  << "\""
                    << ", \"Cluster\" : " << "\"" << pCluster + AlgoExecutionCount[m_inputStageName] << "\""
                    << ", \"CaloHit\" : " << "\"" << pCaloHit + AlgoExecutionCount[m_inputStageName] << "\""
                    << ", \"CaloHitType\" : " << "\""<< pCaloHit->GetHitType() << "\""
                    << ", \"CaloHitX\" : " << pCaloHit->GetPositionVector().GetX()
                    << ", \"CaloHitY\" : " << pCaloHit->GetPositionVector().GetY()
                    << ", \"CaloHitZ\" : " << pCaloHit->GetPositionVector().GetZ()
                    << "},\n";
        }
      } 

      for (const Cluster *pCluster : cluster2DList)
      {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *pCaloHit : caloHitList)
        {
          pfoInfoOutputFile   << "{ \"STAGE\" : " << "\"" << m_inputStageName << "\""
                    << ", \"CallNumber\" : " << "\"" <<AlgoExecutionCount[m_inputStageName] << "\""
                    << ", \"pfoListName\" : " << "\"" << LIST_NAME << "\""
                    << ", \"pfo\" : " << "\"" << pPfo << "\""
                    << ", \"NofCluters\" : " << "\"" << cluster2DList.size()  << "\""
                    << ", \"Cluster\" : " << "\"" << pCluster << "\""
                    << ", \"CaloHit\" : " << "\"" << pCaloHit << "\""
                    << ", \"CaloHitType\" : " << "\""<< pCaloHit->GetHitType() << "\""
                    << ", \"CaloHitX\" : " << pCaloHit->GetPositionVector().GetX()
                    << ", \"CaloHitY\" : " << pCaloHit->GetPositionVector().GetY()
                    << ", \"CaloHitZ\" : " << pCaloHit->GetPositionVector().GetZ()
                    << "},\n";
        }
      }
}

void PrintCurrentPfoInfoAlgorithm::PrintCaloHitsInfo(const CaloHitList *& pCaloHitList, std::string HitsName, std::string STAGE)
{
  if(!pCaloHitList)
  { 
    // std::cout << "pCaloHitList from" << HitsName << "is nullptr\n";
    return;
  }

  if(!pfoInfoOutputFile.is_open())
  {
    std::cout << "Warning: pfoInfoOutputFile not open \n";
    return;
  }
  
  std::cout << "PrintCurrentPfoInfoAlgorithm::PrintCaloHitsInfo start \n";

  for (const CaloHit *const pCaloHit : *pCaloHitList)
  {
    pfoInfoOutputFile << "{\"STAGE\" : " << "\"" << STAGE << "\"" 
            << ", \"CaloHitType\" : " << "\""<< HitsName << "\""
            << ", \"CaloHit\" : " << "\""<< pCaloHit << "\""
            << ", \"CaloHitX\" : " << pCaloHit->GetPositionVector().GetX() 
            << ", \"CaloHitY\" : " << pCaloHit->GetPositionVector().GetY()
            << ", \"CaloHitZ\" : " << pCaloHit->GetPositionVector().GetZ()
            << "}, \n";
  }

  std::cout << "PrintCurrentPfoInfoAlgorithm::PrintCaloHitsInfo end \n";
}

void PrintCurrentPfoInfoAlgorithm::PrintClusterListInfo(const ClusterList*& pClusterList, std::string clusterName, std::string STAGE)
{
  if(!pClusterList || pClusterList->empty())
  {
    // std::cout << __func__ << ": ClusterList " << clusterName << " is either nullptr or empty\n";
    return;
  }

  if(!pfoInfoOutputFile.is_open())
  {
    std::cout << "Warning: pfoInfoOutputFile not open \n";
    return;
  }

  std::cout << "PrintCurrentPfoInfoAlgorithm::PrintClusterListInfo start \n";

  for (const Cluster* pCluster : *pClusterList)
  {
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *pCaloHit : caloHitList)
    {
      pfoInfoOutputFile << "{\"STAGE\" : " << "\"" << STAGE << "\"" 
              << ", \"ClusterName\" : " << "\"" << clusterName << "\"" 
              << ", \"Cluster\" : " << "\"" << pCluster << "\"" 
              << ", \"CaloHit\" : " << "\"" << pCaloHit << "\"" 
              << ", \"CaloHitX\" : " << pCaloHit->GetPositionVector().GetX() 
              << ", \"CaloHitY\" : " << pCaloHit->GetPositionVector().GetY()
              << ", \"CaloHitZ\" : " << pCaloHit->GetPositionVector().GetZ()
              << "}, \n";
    }
  }

  std::cout << "PrintCurrentPfoInfoAlgorithm::PrintClusterListInfo end \n";
  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PrintCurrentPfoInfoAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputStageName", m_inputStageName));

  // hits
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList3DName", m_inputCaloHitList3DName));
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList2DName", m_inputCaloHitList2DName));
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListUName", m_inputCaloHitListUName));
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListVName", m_inputCaloHitListVName));
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListWNames", m_inputCaloHitListWName));

  // clusters
PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName3D", m_inputClusterListName3D));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

  // pfos
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));

  // CR Vertices
  PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputVertexListNames", m_inputVertexListNames));



  return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
