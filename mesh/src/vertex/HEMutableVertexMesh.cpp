/*
 * HEHEMutableVertexMesh.cpp
 *
 *  Created on: 30 Mar 2020
 *      Author: aydar
 */

#include "HEMutableVertexMesh.hpp"


template<unsigned SPACE_DIM>
HEMutableVertexMesh<SPACE_DIM>::HEMutableVertexMesh(std::vector<HENode<SPACE_DIM>*> HENodes,
                                                    std::vector<HEElement<SPACE_DIM>*> HEElements,
                                                    double cellRearrangementThreshold,
                                                    double t2Threshold,
                                                    double cellRearrangementRatio,
                                                    double protorosetteFormationProbability,
                                                    double protorosetteResolutionProbabilityPerTimestep,
                                                    double rosetteResolutionProbabilityPerTimestep)
                                                    :
    HEVertexMesh<SPACE_DIM>(HENodes, HEElements),
    AbstractMutableVertexMesh<SPACE_DIM, SPACE_DIM>(cellRearrangementThreshold,
                                                    cellRearrangementRatio,
                                                    t2Threshold,
                                                    protorosetteFormationProbability,
                                                    protorosetteResolutionProbabilityPerTimestep,
                                                    rosetteResolutionProbabilityPerTimestep)
{
    // Reset member variables and clear mHENodes and mElements
    Clear();

    // If in 3D, then also populate mFaces
    if (SPACE_DIM == 3)
    {
        EXCEPTION("3D MutableVertexMesh is not supported");
    }

    this->mMeshChangesDuringSimulation = true;
}

template<unsigned SPACE_DIM>
HEMutableVertexMesh<SPACE_DIM>::HEMutableVertexMesh()
{
    // Note that the member variables initialised above will be overwritten as soon as archiving is complete
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

template<unsigned SPACE_DIM>
HEMutableVertexMesh<SPACE_DIM>::~HEMutableVertexMesh()
{
    Clear();
}


template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::Clear()
{

}

template<unsigned SPACE_DIM>
bool HEMutableVertexMesh<SPACE_DIM>::CheckForSwapsFromShortEdges()
{
    return false;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::IdentifySwapType(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB)
{

}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::HandleHighOrderJunctions(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB)
{

}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{

}

// Explicit instantiation
template class HEMutableVertexMesh<1>;
template class HEMutableVertexMesh<2>;
template class HEMutableVertexMesh<3>;
