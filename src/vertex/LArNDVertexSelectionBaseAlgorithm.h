/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h
 *
 *  @brief  Header file for the vertex selection base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DUNE_ND_VERTEX_SELECTION_BASE_ALGORITHM_H
#define LAR_DUNE_ND_VERTEX_SELECTION_BASE_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  VertexSelectionBaseAlgorithm class
 */
class DUNE_ND_VertexSelectionBaseAlgorithm : public VertexSelectionBaseAlgorithm
{
protected:
    /**
     *  @brief  Filter the input list of vertices to obtain a reduced number of vertex candidates
     *
     *  @param  pInputVertexList the address of the input vertex list
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  filteredVertices to receive the filtered vertex list
     */
    void FilterVertexList(const pandora::VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
        HitKDTree2D &kdTreeW, pandora::VertexVector &filteredVertices) const override;

    /**
     *  @brief  Calculate the energy of a vertex candidate by summing values from all three planes
     *
     *  @param  pVertex the address of the vertex
     *  @param  kdTreeMap the map of 2D hit kd trees
     *
     *  @return the summed vertex energy
     */
    float GetVertexEnergy(const pandora::Vertex *const pVertex, const KDTreeMap &kdTreeMap) const override;
};

} // namespace lar_content

#endif // #ifndef LAR_DUNE_ND_VERTEX_SELECTION_BASE_ALGORITHM_H
