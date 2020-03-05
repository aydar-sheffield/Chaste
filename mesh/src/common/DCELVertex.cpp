/*
 * DCELVertex.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "DCELVertex.hpp"

#include <cassert>

template<unsigned int SPACE_DIM>
void DCELVertex<SPACE_DIM>::CommonConstructor()
{
    DCELHalfEdge<SPACE_DIM>* mEdge = nullptr;
    DCELVertex<SPACE_DIM>* mNextVertex = nullptr;
    DCELVertex<SPACE_DIM>* mPreviousVertex = nullptr;
    mIsDeleted = false;
}

template<unsigned int SPACE_DIM>
DCELVertex<SPACE_DIM>::~DCELVertex()
{}


template<unsigned int SPACE_DIM>
DCELVertex<SPACE_DIM>::DCELVertex()
{
    for (unsigned int i=0; i<SPACE_DIM; ++i)
    {
        mLocation(i) = 0.0;
    }
    CommonConstructor();
}

template<unsigned int SPACE_DIM>
DCELVertex<SPACE_DIM>::DCELVertex(const c_vector<double, SPACE_DIM> location)
{
    mLocation = location;
    CommonConstructor();
}

template<unsigned int SPACE_DIM>
DCELVertex<SPACE_DIM>::DCELVertex(const double x, const double y, const double z)
{
    mLocation(0) = x;
    mLocation(1) = y;
    if (SPACE_DIM==3)
        mLocation(2) = z;
    CommonConstructor();
}

template <unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>* DCELVertex<SPACE_DIM>::GetOutgoingEdge() const
{
    return mEdge;
}

template <unsigned int SPACE_DIM>
void DCELVertex<SPACE_DIM>::SetOutgoingEdge(DCELHalfEdge<SPACE_DIM>* edge)
{
    mEdge = edge;
}

template<unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& DCELVertex<SPACE_DIM>::rGetLocation() const
{
    // This assert statement is a useful warning: when new nodes are created we overwrite previously deleted nodes if there are any.
    // This means that we can not use this method to interrogate deleted nodes about their position before deletion because we can't
    // guarantee that the node has not been overwritten already. Hence, when implementing new functionality we need to make sure
    // that this functionality does not rely on being able to interrogate deleted nodes for their location.
    // \todo #2401: make this an exception.
    assert(!mIsDeleted);
    return mLocation;
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM>& DCELVertex<SPACE_DIM>::rGetModifiableLocation()
{
    assert(!mIsDeleted);
    return mLocation;
}

template class DCELVertex<1>;
template class DCELVertex<2>;
template class DCELVertex<3>;



