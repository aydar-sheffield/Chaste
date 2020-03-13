/*
 * DCELEdge.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HalfEdge.hpp"

#include "HEElement.hpp"
#include "HEVertex.hpp"

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>::HalfEdge()
:
mTwin(nullptr), mNextEdge(nullptr), mPreviousEdge(nullptr),
mTargetVertex(nullptr), mElement(nullptr),mIndex(0)
{
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>::~HalfEdge()
{}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetTwinHalfEdge() const
{
    return mTwin;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetTwinHalfEdge(HalfEdge<SPACE_DIM>* edge)
{
    mTwin = edge;
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetNextHalfEdge() const
{
    return mNextEdge;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetNextHalfEdge(HalfEdge<SPACE_DIM>* edge)
{
    mNextEdge = edge;
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetPreviousHalfEdge() const
{
    return mPreviousEdge;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetPreviousHalfEdge(HalfEdge<SPACE_DIM>* edge)
{
    mPreviousEdge = edge;
}

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetOriginVertex() const
{
    return mTwin->GetTargetVertex();
}

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetTargetVertex() const
{
    return mTargetVertex;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetTargetVertex(HEVertex<SPACE_DIM>* vertex)
{
    mTargetVertex = vertex;
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetElement() const
{
    return mElement;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetElement(HEElement<SPACE_DIM>* element)
{
    mElement = element;
}

template<unsigned int SPACE_DIM>
unsigned int HalfEdge<SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetIndex(const unsigned new_index)
{
    mIndex = new_index;
}
template class HalfEdge<1>;
template class HalfEdge<2>;
template class HalfEdge<3>;

