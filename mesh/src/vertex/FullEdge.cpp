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
FullEdge<SPACE_DIM>::FullEdge(const unsigned index)
:
mpFirst(nullptr),
mpSecond(nullptr),
mIndex(index)
{}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::FullEdge(HalfEdge<SPACE_DIM>* edge0, HalfEdge<SPACE_DIM>* edge1, const unsigned index)
:
mpFirst(edge0),
mpSecond(edge1),
mIndex(index)
{}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::FullEdge(HalfEdge<SPACE_DIM>* edge0, const unsigned index)
:
mpFirst(edge0),
mpSecond(edge0->GetTwinHalfEdge()),
mIndex(index)
{}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>::~FullEdge()
{}

template<unsigned int SPACE_DIM>
double FullEdge<SPACE_DIM>::GetLength() const
{
    assert(mpFirst->GetTwinHalfEdge()==mpSecond);
    assert(mpFirst->GetLength()==mpSecond->GetLength());
    return mpFirst->GetLength();
}

template<unsigned int SPACE_DIM>
double FullEdge<SPACE_DIM>::ComputeLength()
{
    return mpFirst->ComputeLength();
}

template<unsigned int SPACE_DIM>
unsigned int FullEdge<SPACE_DIM>::GetIndex() const
{
    return mIndex;
}

template<unsigned int SPACE_DIM>
void FullEdge<SPACE_DIM>::SetIndex(const unsigned int new_index)
{
    mIndex = new_index;
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

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* FullEdge<SPACE_DIM>::at(unsigned int index) const
{
    assert(index<=1);
    return index == 0 ? mpFirst : mpSecond;
}
template<unsigned int SPACE_DIM>
std::pair<HENode<SPACE_DIM>*, HENode<SPACE_DIM>* > FullEdge<SPACE_DIM>::GetNodes() const
{
    return std::pair<HENode<SPACE_DIM>*, HENode<SPACE_DIM>* >(mpFirst->GetTargetNode(), mpSecond->GetTargetNode());
}
template class FullEdge<1>;
template class FullEdge<2>;
template class FullEdge<3>;



