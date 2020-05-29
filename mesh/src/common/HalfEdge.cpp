/*
 * DCELEdge.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HalfEdge.hpp"

#include "HEElement.hpp"
#include "HENode.hpp"

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>::HalfEdge(HEElement<SPACE_DIM>* pElement)
:
mpTwin(nullptr), mpNextEdge(nullptr), mpPreviousEdge(nullptr),
mpTargetNode(nullptr), mpElement(pElement),mIndex(0), mIsDeleted(false), mLength(0)
{
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>::~HalfEdge()
{}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetTwinHalfEdge() const
{
    return mpTwin;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetTwinHalfEdge(HalfEdge<SPACE_DIM>* pEdge, const bool SetOtherEdgeToo)
{
    mpTwin = pEdge;
    if (SetOtherEdgeToo)
        pEdge->SetTwinHalfEdge(this);
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetNextHalfEdge() const
{
    return mpNextEdge;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetNextHalfEdge(HalfEdge<SPACE_DIM>* pEdge, const bool SetOtherEdgeToo)
{
    mpNextEdge = pEdge;
    if (SetOtherEdgeToo)
        pEdge->SetPreviousHalfEdge(this);
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetPreviousHalfEdge() const
{
    return mpPreviousEdge;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetPreviousHalfEdge(HalfEdge<SPACE_DIM>* pEdge, const bool SetOtherEdgeToo)
{
    mpPreviousEdge = pEdge;
    if (SetOtherEdgeToo)
        pEdge->SetNextHalfEdge(this);
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetOriginNode() const
{
    //return mpPreviousEdge->GetTargetNode();
    return mpTwin->GetTargetNode();
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetOriginNode(HENode<SPACE_DIM>* pNode)
{
    GetPreviousHalfEdge()->SetTargetNode(pNode);
    GetTwinHalfEdge()->SetTargetNode(pNode);
    pNode->SetOutgoingEdge(this);
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetTargetNode() const
{
    return mpTargetNode;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetTargetNode(HENode<SPACE_DIM>* pNode, const bool ModifyAdjacentEdges)
{
    mpTargetNode = pNode;
    if (ModifyAdjacentEdges)
        GetNextHalfEdge()->GetTwinHalfEdge()->SetTargetNode(pNode);
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>* HalfEdge<SPACE_DIM>::GetElement() const
{
    return mpElement;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetElement(HEElement<SPACE_DIM>* pElement)
{
    mpElement = pElement;
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

template<unsigned int SPACE_DIM>
bool HalfEdge<SPACE_DIM>::IsDeleted() const
{
    return mIsDeleted;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::SetDeletedStatus(const bool status, const bool SetTwin)
{
    mIsDeleted = status;
    if (SetTwin)
        mpTwin->SetDeletedStatus(status);
}

template<unsigned int SPACE_DIM>
bool HalfEdge<SPACE_DIM>::IsFullyInitialized() const
{
    return mpTwin&&mpNextEdge&&mpPreviousEdge&&mpTargetNode;
}

template<unsigned int SPACE_DIM>
bool HalfEdge<SPACE_DIM>::IsOnBoundary() const
{
    return !mpTwin->GetElement();
}

template<unsigned int SPACE_DIM>
double HalfEdge<SPACE_DIM>::ComputeLength() const
{
    c_vector<double, SPACE_DIM> edge_vector;
    edge_vector = GetOriginNode()->rGetLocation()-GetOriginNode()->rGetLocation();
    return norm_2(edge_vector);
}

template<unsigned int SPACE_DIM>
inline double HalfEdge<SPACE_DIM>::GetLength()
{
    return mLength;
}

template<unsigned int SPACE_DIM>
void HalfEdge<SPACE_DIM>::UpdateLength(const double new_length)
{
    mLength = new_length==0 ? ComputeLength() : new_length;
}
template class HalfEdge<1>;
template class HalfEdge<2>;
template class HalfEdge<3>;

