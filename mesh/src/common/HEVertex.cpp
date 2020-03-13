/*
 * HEVertex.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HEVertex.hpp"
#include "HEElement.hpp"
#include <cassert>

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>::HEVertex(unsigned index, std::vector<double> coords,
                              bool isBoundaryNode, HalfEdge<SPACE_DIM>* edge)
                              : Node<SPACE_DIM>(index, coords, isBoundaryNode),
                                mEdge(edge)
{
}

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>::HEVertex(unsigned index, c_vector<double, SPACE_DIM> location,
                              bool isBoundaryNode, HalfEdge<SPACE_DIM>* edge)
                              : Node<SPACE_DIM>(index, location, isBoundaryNode),
                               mEdge(edge)
{
}

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>::HEVertex(unsigned index, bool isBoundaryNode, double v1, double v2, double v3,
                              HalfEdge<SPACE_DIM>* edge)
                              : Node<SPACE_DIM>(index, isBoundaryNode, v1, v2, v3),
                               mEdge(edge)
{
}

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>::HEVertex(const Node<SPACE_DIM> &node, HalfEdge<SPACE_DIM>* edge)
:Node<SPACE_DIM>(node.GetIndex(), node.rGetLocation(), node.IsBoundaryNode()),
 mEdge(edge)
{
}

template<unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>::~HEVertex()
{}

template <unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEVertex<SPACE_DIM>::GetOutgoingEdge() const
{
    return mEdge;
}

template <unsigned int SPACE_DIM>
void HEVertex<SPACE_DIM>::SetOutgoingEdge(HalfEdge<SPACE_DIM>* edge)
{
    mEdge = edge;
}

template <unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEVertex<SPACE_DIM>::GetIncomingEdge() const
{
    return mEdge->GetTwinHalfEdge();
}

template <unsigned int SPACE_DIM>
HEVertex<SPACE_DIM>* HEVertex<SPACE_DIM>::GetNextVertex() const
{
    return mEdge->GetTargetVertex();
}

template <unsigned int SPACE_DIM>
void HEVertex<SPACE_DIM>::UpdateElementIndices()
{
    std::set<unsigned int> &element_indices = this->rGetContainingElementIndices();
    element_indices.clear();
    HalfEdge<SPACE_DIM>* edge = mEdge;
    HEElement<SPACE_DIM>* element;
    do
    {
        element = edge->GetElement();
        if (element)
            element_indices.insert(element->GetIndex());
        edge = edge->GetTwinHalfEdge()->GetNextHalfEdge();
    }
    while (edge != mEdge);
}
template class HEVertex<1>;
template class HEVertex<2>;
template class HEVertex<3>;



