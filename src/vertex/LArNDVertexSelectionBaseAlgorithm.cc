/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the vertex selection base algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "LArNDVertexSelectionBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void DUNE_ND_VertexSelectionBaseAlgorithm::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexVector &filteredVertices) const
{
    for (const Vertex *const pVertex : *pInputVertexList)
    {
        unsigned int nAcceptableViews(0);

        if ((m_isEmptyViewAcceptable && kdTreeU.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeV.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeW.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW))
            ++nAcceptableViews;

        if (nAcceptableViews >= m_minVertexAcceptableViews)
            filteredVertices.push_back(pVertex);
    }

    std::sort(filteredVertices.begin(), filteredVertices.end(), SortByVertexZPosition);
}

float DUNE_ND_VertexSelectionBaseAlgorithm::GetVertexEnergy(const pandora::Vertex *const pVertex, const KDTreeMap &kdTreeMap) const
{
    float totalEnergy(0.f);

    if (this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeMap.at(TPC_VIEW_U)))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_U, kdTreeMap.at(TPC_VIEW_U));

    if (this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeMap.at(TPC_VIEW_V)))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_V, kdTreeMap.at(TPC_VIEW_V));

    if (this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeMap.at(TPC_VIEW_W)))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_W, kdTreeMap.at(TPC_VIEW_W));

    return totalEnergy;
}

} // namespace lar_content
