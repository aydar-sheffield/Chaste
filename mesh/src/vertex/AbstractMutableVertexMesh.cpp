/*
 * AbstractAbstractMutableVertexMesh.cpp
 *
 *  Created on: 30 Mar 2020
 *      Author: aydar
 */

#include "AbstractMutableVertexMesh.hpp"
#include "Exception.hpp"
#include <cassert>
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::AbstractMutableVertexMesh(double cellRearrangementThreshold,
                                                                             double t2Threshold,
                                                                             double cellRearrangementRatio,
                                                                             double protorosetteFormationProbability,
                                                                             double protorosetteResolutionProbabilityPerTimestep,
                                                                             double rosetteResolutionProbabilityPerTimestep)
        : mCellRearrangementThreshold(cellRearrangementThreshold),
          mCellRearrangementRatio(cellRearrangementRatio),
          mT2Threshold(t2Threshold),
          mProtorosetteFormationProbability(protorosetteFormationProbability),
          mProtorosetteResolutionProbabilityPerTimestep(protorosetteResolutionProbabilityPerTimestep),
          mRosetteResolutionProbabilityPerTimestep(rosetteResolutionProbabilityPerTimestep),
          mCheckForInternalIntersections(false),
          mDistanceForT3SwapChecking(5.0)
{
    // Threshold parameters must be strictly positive
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);
    assert(protorosetteFormationProbability >= 0.0);
    assert(protorosetteFormationProbability <= 1.0);
    assert(protorosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(protorosetteResolutionProbabilityPerTimestep <= 1.0);
    assert(rosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(rosetteResolutionProbabilityPerTimestep <= 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractMutableVertexMesh()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementThreshold() const
{
    return mCellRearrangementThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetT2Threshold() const
{
    return mT2Threshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementRatio() const
{
    return mCellRearrangementRatio;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteFormationProbability() const
{
    return this->mProtorosetteFormationProbability;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteResolutionProbabilityPerTimestep() const
{
    return this->mProtorosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetRosetteResolutionProbabilityPerTimestep() const
{
    return this->mRosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetDistanceForT3SwapChecking(double distanceForT3SwapChecking)
{
    mDistanceForT3SwapChecking = distanceForT3SwapChecking;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceForT3SwapChecking() const
{
    return mDistanceForT3SwapChecking;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCheckForInternalIntersections() const
{
    return mCheckForInternalIntersections;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementThreshold(double cellRearrangementThreshold)
{
    mCellRearrangementThreshold = cellRearrangementThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetT2Threshold(double t2Threshold)
{
    mT2Threshold = t2Threshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementRatio(double cellRearrangementRatio)
{
    mCellRearrangementRatio = cellRearrangementRatio;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteFormationProbability(double protorosetteFormationProbability)
{
    // Check that the new value is in [0,1]
    if (protorosetteFormationProbability < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteFormationProbability > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteFormationProbability = protorosetteFormationProbability;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (protorosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteResolutionProbabilityPerTimestep = protorosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (rosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (rosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mRosetteResolutionProbabilityPerTimestep = rosetteResolutionProbabilityPerTimestep;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCheckForInternalIntersections(bool checkForInternalIntersections)
{
    mCheckForInternalIntersections = checkForInternalIntersections;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMutableVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{
}
// Explicit instantiation
template class AbstractMutableVertexMesh<1, 1>;
template class AbstractMutableVertexMesh<1, 2>;
template class AbstractMutableVertexMesh<1, 3>;
template class AbstractMutableVertexMesh<2, 2>;
template class AbstractMutableVertexMesh<2, 3>;
template class AbstractMutableVertexMesh<3, 3>;


