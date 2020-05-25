/*
 * HEElement.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HEElement.hpp"

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::CommonConstructor(const std::vector<HENode<SPACE_DIM>* > node_list)
{
    const unsigned int n_nodes = node_list.size();

    for (unsigned int i=0; i<n_nodes; ++i)
    {
        HalfEdge<SPACE_DIM>* out_edge = node_list[i]->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* edge;
        const unsigned int next_index = (i+1)%n_nodes;
        /**
         * If neighbouring halfedges are constructed, get the edge outgoing from this node
         */
        bool construct_edge = true;
        if (out_edge)
        {
            if (out_edge->GetPreviousHalfEdge()->GetTwinHalfEdge()->GetTargetNode()==node_list[next_index])
            {
                edge = out_edge->GetPreviousHalfEdge()->GetTwinHalfEdge();
                construct_edge = false;
            }
        }

        if (construct_edge)
        {
            edge = new HalfEdge<SPACE_DIM>();

            edge->SetTargetNode(node_list[next_index]);

            HalfEdge<SPACE_DIM>* twin_edge = new HalfEdge<SPACE_DIM>();
            twin_edge->SetTargetNode(node_list[i]);
            twin_edge->SetTwinHalfEdge(edge);
            edge->SetTwinHalfEdge(twin_edge);
        }
        node_list[i]->SetOutgoingEdge(edge);
        edge->SetElement(this);
    }


    for (unsigned int i=0; i<n_nodes; ++i)
    {
        //Link internal edges
        HalfEdge<SPACE_DIM>* out_edge = node_list[i]->GetOutgoingEdge();
        const unsigned int next_index = (i+1)%n_nodes;
        const unsigned int previous_index = (i-1+n_nodes)%n_nodes;
        HalfEdge<SPACE_DIM>* next_edge = node_list[next_index]->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* previous_edge = node_list[previous_index]->GetOutgoingEdge();
        out_edge->SetNextHalfEdge(next_edge);
        out_edge->SetPreviousHalfEdge(previous_edge);

        //Linking external (CW or boundary) edges
        HalfEdge<SPACE_DIM>* out_edge_twin = out_edge->GetTwinHalfEdge();
        HalfEdge<SPACE_DIM>* next_edge_twin = next_edge -> GetTwinHalfEdge();
        HalfEdge<SPACE_DIM>* previous_edge_twin = previous_edge -> GetTwinHalfEdge();

        //If on boundary
        if (!out_edge_twin->GetElement())
        {
            while(previous_edge_twin->GetElement())
                previous_edge_twin = previous_edge_twin->GetPreviousHalfEdge()->GetTwinHalfEdge();
            out_edge_twin->SetNextHalfEdge(previous_edge_twin);
            previous_edge_twin->SetPreviousHalfEdge(out_edge_twin);

            while(next_edge_twin->GetElement())
                next_edge_twin = next_edge_twin->GetNextHalfEdge()->GetTwinHalfEdge();
            out_edge_twin->SetPreviousHalfEdge(next_edge_twin);
            next_edge_twin->SetNextHalfEdge(out_edge_twin);
        }
    }

    mpHalfEdge = node_list[0]->GetOutgoingEdge();
    mNumNodes = n_nodes;
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index)
:AbstractElement<SPACE_DIM,SPACE_DIM>(index),
 mNumNodes(0)
{
    mpHalfEdge = nullptr;
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<HENode<SPACE_DIM>* > node_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index),
 mNumNodes(node_list.size())
{
    mpHalfEdge = nullptr;
    CommonConstructor(node_list);
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, HalfEdge<SPACE_DIM>* edge)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index),
mpHalfEdge(edge)
{
    RegisterWithHalfEdges();
}

template <unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<Node<SPACE_DIM>* > node_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index),
 mNumNodes(node_list.size())
 {
    std::vector<HENode<SPACE_DIM>* > henode_list;
    for(auto node:node_list)
    {
        henode_list.push_back(new HENode<SPACE_DIM>(*node));
    }
    mpHalfEdge = nullptr;
    CommonConstructor(henode_list);
 }

template <unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(const VertexElement<SPACE_DIM, SPACE_DIM> &rElement)
:AbstractElement<SPACE_DIM, SPACE_DIM>(rElement.GetIndex()),
 mNumNodes(rElement.GetNumNodes())
 {
    std::vector<HENode<SPACE_DIM>* > node_list;
    for (unsigned int i=0; i<rElement.GetNumNodes(); ++i)
    {
        node_list.push_back(new HENode<SPACE_DIM>(*rElement.GetNode(i)));
    }
    mpHalfEdge = nullptr;
    CommonConstructor(node_list);
 }

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::~HEElement()
{}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::GetHalfEdge(const unsigned int local_index) const
{
    assert(local_index<mNumNodes);
    //Outgoing internal edge from node local_index
    HalfEdge<SPACE_DIM>* out_edge = mpHalfEdge;
    for (unsigned int i=0; i<local_index; ++i)
    {
        out_edge = out_edge->GetNextHalfEdge();
    }
    return out_edge;
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::GetHalfEdge(HENode<SPACE_DIM>* pTarget) const
{
    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
    if (next_edge->GetTargetNode()!=pTarget)
    {
        do
        {
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge&&pTarget!=next_edge->GetTargetNode());
        if (next_edge == mpHalfEdge)
        {
            EXCEPTION("Halfedge search by node not found.");
        }
    }
    return next_edge;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::SetHalfEdge(HalfEdge<SPACE_DIM>* pEdge)
{
    mpHalfEdge = pEdge;
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HEElement<SPACE_DIM>::GetNode(unsigned int local_index) const
{
    HalfEdge<SPACE_DIM>* edge = GetHalfEdge(local_index);
    return edge->GetOriginNode();
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::UpdateNode(const unsigned& rIndex, HENode<SPACE_DIM>* pNode)
{
    assert(rIndex<mNumNodes);
    //Outgoing internal edge from node rIndex
    HalfEdge<SPACE_DIM>* out_edge = GetHalfEdge(rIndex);
    out_edge->SetOriginNode(pNode);
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::ReplaceNode(HENode<SPACE_DIM>* pOldNode, HENode<SPACE_DIM>* pNewNode)
{
    assert(pOldNode != pNewNode);

    HalfEdge<SPACE_DIM>* next_edge = GetHalfEdge(pOldNode)->GetNextHalfEdge();

    next_edge->SetOriginNode(pNewNode);
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::MarkAsDeleted()
{
    this->mIsDeleted = true;

    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;

    do
    {
        next_edge->SetElement(nullptr);
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge!=mpHalfEdge);

}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::RegisterWithHalfEdges()
{
    assert(mpHalfEdge);
    HalfEdge<SPACE_DIM>* edge = mpHalfEdge;
    mNumNodes = 0;
    do
    {
        mNumNodes++;
        edge->SetElement(this);
        edge = edge->GetNextHalfEdge();
    }while(edge!=mpHalfEdge);
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::AddNode(HalfEdge<SPACE_DIM>* pEdge, HENode<SPACE_DIM>* pNode)
{
    HalfEdge<SPACE_DIM>* prev_edge = pEdge->GetPreviousHalfEdge();

    //A new halfedge between prev_edge->OriginNode() and pNode is created with pNode as its target index.
    HalfEdge<SPACE_DIM>* new_edge = new HalfEdge<SPACE_DIM>(this);
    //Set the element to which the new twin edge belongs to
    HalfEdge<SPACE_DIM>* new_edge_twin = new HalfEdge<SPACE_DIM>(pEdge->GetTwinHalfEdge()->GetElement());

    //Set edge adjacency relations
    new_edge->SetTwinHalfEdge(new_edge_twin, true);
    new_edge->SetNextHalfEdge(pEdge, true);
    new_edge->SetPreviousHalfEdge(prev_edge, true);

    //External halfedges
    new_edge_twin->SetNextHalfEdge(prev_edge->GetTwinHalfEdge(), true);
    new_edge_twin->SetPreviousHalfEdge(pEdge->GetTwinHalfEdge(), true);

    //Set target Node relations
    pEdge->SetOriginNode(pNode);
    new_edge_twin->SetTargetNode(prev_edge->GetTargetNode());

    assert(new_edge->IsFullyInitialized());
    assert(new_edge_twin->IsFullyInitialized());

    if (pEdge==mpHalfEdge)
        mpHalfEdge = new_edge;
    mNumNodes++;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::DeleteNode(HENode<SPACE_DIM>* pNode)
{
    //Edge outgoing from pNode is deleted

    //Find incoming edge
    HalfEdge<SPACE_DIM>* in_edge = GetHalfEdge(pNode);
    HalfEdge<SPACE_DIM>* in_edge_twin = in_edge->GetTwinHalfEdge();

    //Outgoing edge from pNode
    HalfEdge<SPACE_DIM>* out_edge = in_edge->GetNextHalfEdge();
    HalfEdge<SPACE_DIM>* next_out_edge = out_edge->GetNextHalfEdge();

    in_edge->SetNextHalfEdge(next_out_edge, true);
    in_edge_twin->SetPreviousHalfEdge(next_out_edge->GetTwinHalfEdge(), true);

    in_edge->SetTargetNode(out_edge->GetTargetNode());

    //Make sure the edge outgoing from the node after pNode is not the deleted one
    if (out_edge->GetTargetNode()->GetOutgoingEdge()==out_edge->GetTwinHalfEdge())
        out_edge->GetTargetNode()->SetOutgoingEdge(in_edge_twin);

    out_edge->SetDeletedStatus(true,true);
}

template<unsigned int SPACE_DIM>
unsigned int HEElement<SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned int SPACE_DIM>
bool HEElement<SPACE_DIM>::IsElementOnBoundary() const
{
    bool is_element_on_boundary = false;
    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
    do
    {
        is_element_on_boundary = next_edge->GetTargetNode()->IsBoundaryNode();
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge!=mpHalfEdge &&!is_element_on_boundary);

    return is_element_on_boundary;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::GetVolume() const
{
    return mVolume;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::GetSurfaceArea() const
{
    return mSurfaceArea;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::ComputeVolume()
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    double element_volume = 0.0;
    if (SPACE_DIM == 2)
    {
        // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location;
        first_node_location = mpHalfEdge->GetOriginNode()->rGetLocation();
        c_vector<double, SPACE_DIM> pos_1;
        pos_1 = zero_vector<double>(SPACE_DIM);

        // Loop over vertices
        HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
        do
        {
            c_vector<double, SPACE_DIM> next_node_location = next_edge->GetTargetNode()->rGetLocation();
            c_vector<double, SPACE_DIM> pos_2 = this->GetVectorFromAtoB(first_node_location, next_node_location);

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            element_volume += 0.5 * (this_x * next_y - next_x * this_y);

            pos_1 = pos_2;
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge);
    }
    else
    {
        //3D case not supported
        EXCEPTION("Half-edge mesh in 3D not supported.");
    }
    // We take the absolute value just in case the nodes were really oriented clockwise
    mVolume = fabs(element_volume);

    return mVolume;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::ComputeSurfaceArea()
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    double surface_area = 0.0;
    if (SPACE_DIM == 2)
    {
        HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
        do
        {
            surface_area += next_edge->ComputeLength();
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge);
    }
    else
    {
        //3D case not supported
        EXCEPTION("Half-edge mesh in 3D not supported.");
    }
    mSurfaceArea = surface_area;
    return surface_area;
}



template class HEElement<1>;
template class HEElement<2>;
template class HEElement<3>;

