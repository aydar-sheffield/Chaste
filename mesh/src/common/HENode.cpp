/*
 * HENode.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HENode.hpp"
#include "HEElement.hpp"
#include <cassert>

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>::HENode(unsigned index, std::vector<double> coords,
                              bool isBoundaryNode, HalfEdge<SPACE_DIM>* pEdge)
                              : Node<SPACE_DIM>(index, coords, isBoundaryNode),
                                mpEdge(pEdge)
{
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>::HENode(unsigned index, c_vector<double, SPACE_DIM> location,
                              bool isBoundaryNode, HalfEdge<SPACE_DIM>* pEdge)
                              : Node<SPACE_DIM>(index, location, isBoundaryNode),
                               mpEdge(pEdge)
{
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>::HENode(unsigned index, bool isBoundaryNode, double v1, double v2, double v3,
                              HalfEdge<SPACE_DIM>* pEdge)
                              : Node<SPACE_DIM>(index, isBoundaryNode, v1, v2, v3),
                               mpEdge(pEdge)
{
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>::HENode(const Node<SPACE_DIM> &rNode, HalfEdge<SPACE_DIM>* pEdge)
:Node<SPACE_DIM>(rNode.GetIndex(), rNode.rGetLocation(), rNode.IsBoundaryNode()),
 mpEdge(pEdge)
{
}

template<unsigned SPACE_DIM>
HENode<SPACE_DIM>::HENode(unsigned index, ChastePoint<SPACE_DIM> point, bool isBoundaryNode)
:Node<SPACE_DIM>(index, point, isBoundaryNode),
 mpEdge(nullptr)
{
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>::~HENode()
{
}

template <unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HENode<SPACE_DIM>::GetOutgoingEdge() const
{
    return mpEdge;
}

template <unsigned int SPACE_DIM>
void HENode<SPACE_DIM>::SetOutgoingEdge(HalfEdge<SPACE_DIM>* pEdge)
{
    mpEdge = pEdge;
}

template <unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HENode<SPACE_DIM>::GetIncomingEdge() const
{
    return mpEdge->GetTwinHalfEdge();
}

template <unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HENode<SPACE_DIM>::GetNextNode() const
{
    return mpEdge->GetTargetNode();
}

template <unsigned int SPACE_DIM>
std::set<unsigned int> HENode<SPACE_DIM>::GetContainingElementIndices() const
{
    std::set<unsigned int> element_indices;
    if (mpEdge)
    {
        assert(mpEdge!=nullptr);
        HalfEdge<SPACE_DIM>* edge = mpEdge;
        HEElement<SPACE_DIM>* element;
        //Iterate over outgoing edges...
        do
        {
            //.. and get the elements to which they belong
            element = edge->GetElement();
            if (element)
                element_indices.insert(element->GetIndex());
            edge = edge->GetTwinHalfEdge()->GetNextHalfEdge();
        }
        while (edge != mpEdge);
    }
    return element_indices;
}

template <unsigned int SPACE_DIM>
unsigned int HENode<SPACE_DIM>::GetNumContainingElements() const
{
    std::set<unsigned int> indices = GetContainingElementIndices();
    return indices.size();
}

template class HENode<1>;
template class HENode<2>;
template class HENode<3>;



