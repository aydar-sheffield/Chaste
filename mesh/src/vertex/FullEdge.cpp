/*
 * FullEdge.cpp
 *
 *  Created on: 25 May 2020
 *      Author: aydar
 */

#include "FullEdge.hpp"
#include "HEElement.hpp"
#include <cassert>
template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::FullEdge()
:
mpFirst(nullptr),
mpSecond(nullptr),
mLength(0)
{}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::FullEdge(HalfEdge<SPACE_DIM>* edge0, HalfEdge<SPACE_DIM>* edge1)
:
mpFirst(edge0),
mpSecond(edge1),
mLength(0)
{}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::FullEdge(HalfEdge<SPACE_DIM>* edge0)
:
mpFirst(edge0),
mpSecond(edge0->GetTwinHalfEdge()),
mLength(0)
{}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::~FullEdge()
{}

template<unsigned int SPACE_DIM>
double FullEdge<SPACE_DIM>::GetLength() const
{
    return mLength;
}

template<unsigned int SPACE_DIM>
double FullEdge<SPACE_DIM>::ComputeLength()
{
    mLength = mpFirst->ComputeLength();
    mpSecond->UpdateLength(mLength);
    return mLength;
}

template<unsigned int SPACE_DIM>
std::set<unsigned int> FullEdge<SPACE_DIM>::GetContainingElementIndices() const
{
    std::set<unsigned int> indices;
    auto element_0 = mpFirst->GetElement();
    auto element_1 = mpSecond->GetElement();
    if (element_0)
        indices.insert(element_0->GetIndex());
    if (element_1)
        indices.insert(element_1->GetIndex());
    return indices;
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* FullEdge<SPACE_DIM>::operator()(unsigned int index)
{
    assert(index<=1);
    return index == 0 ? mpFirst : mpSecond;
}

template class FullEdge<1>;
template class FullEdge<2>;
template class FullEdge<3>;



