/*
 * HEElement.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HEElement.hpp"

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::CommonConstructor(const std::vector<HEVertex<SPACE_DIM>* > vertex_list)
{
    const unsigned int n_nodes = vertex_list.size();

    for (unsigned int i=0; i<n_nodes; ++i)
    {
        HalfEdge<SPACE_DIM>* out_edge = vertex_list[i]->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* edge;
        const unsigned int next_index = (i+1)%n_nodes;
        /**
         * If neighbouring halfedges are constructed, get the edge outgoing from this vertex
         */
        bool construct_edge = true;
        if (out_edge)
        {
            if (out_edge->GetPreviousHalfEdge()->GetTwinHalfEdge()->GetTargetVertex()==vertex_list[next_index])
            {
                edge = out_edge->GetPreviousHalfEdge()->GetTwinHalfEdge();
                construct_edge = false;
            }
        }

        if (construct_edge)
        {
            edge = new HalfEdge<SPACE_DIM>();

            edge->SetTargetVertex(vertex_list[next_index]);

            HalfEdge<SPACE_DIM>* twin_edge = new HalfEdge<SPACE_DIM>();
            twin_edge->SetTargetVertex(vertex_list[i]);
            twin_edge->SetTwinHalfEdge(edge);
            edge->SetTwinHalfEdge(twin_edge);
        }
        vertex_list[i]->SetOutgoingEdge(edge);
        edge->SetElement(this);
    }


    for (unsigned int i=0; i<n_nodes; ++i)
    {
        //Link internal edges
        HalfEdge<SPACE_DIM>* out_edge = vertex_list[i]->GetOutgoingEdge();
        const unsigned int next_index = (i+1)%n_nodes;
        const unsigned int previous_index = (i-1+n_nodes)%n_nodes;
        HalfEdge<SPACE_DIM>* next_edge = vertex_list[next_index]->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* previous_edge = vertex_list[previous_index]->GetOutgoingEdge();
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

    mHalfEdge = vertex_list[0]->GetOutgoingEdge();
    mNumNodes = vertex_list.size();
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index)
:AbstractElement<SPACE_DIM,SPACE_DIM>(index),
 mNumNodes(0)
{
    mHalfEdge = nullptr;
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<HEVertex<SPACE_DIM>* > vertex_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index)
{
    mHalfEdge = nullptr;
    CommonConstructor(vertex_list);
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<HalfEdge<SPACE_DIM>* > edge_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index)
{
    mHalfEdge = edge_list[0];
    const unsigned int n_edges = edge_list.size();
    for (unsigned int i=0; i<n_edges; ++i)
    {
        const unsigned int next_index = (i+1)%n_edges;
        const unsigned int prev_index = (i+n_edges-1)%n_edges;
        edge_list[i]->SetElement(this);
        edge_list[i]->SetPreviousHalfEdge(edge_list[prev_index]);
        edge_list[i]->SetNextHalfEdge(edge_list[next_index]);
    }
}

template <unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<Node<SPACE_DIM>* > node_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index)
 {
    std::vector<HEVertex<SPACE_DIM>* > vertex_list;
    for(auto node:node_list)
    {
        vertex_list.push_back(new HEVertex<SPACE_DIM>(*node));
    }
    mHalfEdge = nullptr;
    CommonConstructor(vertex_list);
 }

template <unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(const VertexElement<SPACE_DIM, SPACE_DIM> &element)
:AbstractElement<SPACE_DIM, SPACE_DIM>(element.GetIndex())
 {
    std::vector<HEVertex<SPACE_DIM>* > vertex_list;
    for (unsigned int i=0; i<element.GetNumNodes(); ++i)
    {
        vertex_list.push_back(new HEVertex<SPACE_DIM>(*element.GetNode(i)));
    }
    mHalfEdge = nullptr;
    CommonConstructor(vertex_list);
 }

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::~HEElement()
{}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::GetHalfEdge() const
{
    return mHalfEdge;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::SetHalfEdge(HalfEdge<SPACE_DIM>* edge)
{
    mHalfEdge = edge;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::UpdateNode(const unsigned& rIndex, HEVertex<SPACE_DIM>* pNode)
{

}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::MarkAsDeleted()
{
    this->mIsDeleted = true;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::RegisterWithNodes()
{

}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::AddNode(const unsigned int &rIndex, HEVertex<SPACE_DIM>* pNode)
{
    HalfEdge<SPACE_DIM>* prev_edge, next_edge;
    if (rIndex==0)
    {

    }
    HalfEdge<SPACE_DIM>* walker = mHalfEdge->GetNextHalfEdge();
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::UpdateNumNodes()
{
    mNumNodes = 0;
    HalfEdge<SPACE_DIM>* next_edge = mHalfEdge;
    do
    {
        mNumNodes++;
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != mHalfEdge);
}

template<unsigned int SPACE_DIM>
unsigned int HEElement<SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}
template class HEElement<1>;
template class HEElement<2>;
template class HEElement<3>;

