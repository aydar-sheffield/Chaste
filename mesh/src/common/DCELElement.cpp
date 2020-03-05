/*
 * DCELElement.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "DCELElement.hpp"

template<unsigned int SPACE_DIM>
DCELElement<SPACE_DIM>::DCELElement()
{
    mHalfEdge = nullptr;
}

template<unsigned int SPACE_DIM>
DCELElement<SPACE_DIM>::DCELElement(std::vector<DCELVertex<SPACE_DIM>* > vertex_list)
{
    const unsigned int n_nodes = vertex_list.size();

    for (unsigned int i=0; i<n_nodes; ++i)
    {
        DCELHalfEdge<SPACE_DIM>* out_edge = vertex_list[i]->GetOutgoingEdge();
        DCELHalfEdge<SPACE_DIM>* edge;
        const unsigned int next_index = (i+1)%n_nodes;
        bool construct_edge = true;
        if (!out_edge)
        {
            if (out_edge->GetPreviousHalfEdge()->GetTwinHalfEdge()->GetTargetVertex()==vertex_list[next_index])
            {
                edge = out_edge->GetPreviousHalfEdge()->GetTwinHalfEdge();
                construct_edge = false;
            }
        }
        if (construct_edge)
        {
            edge = new DCELHalfEdge<SPACE_DIM>();

            edge->SetTargetVertex(vertex_list[next_index]);

            DCELHalfEdge<SPACE_DIM>* twin_edge = new DCELHalfEdge<SPACE_DIM>();
            twin_edge->SetTargetVertex(vertex_list[i]);
            twin_edge->SetTwinHalfEdge(edge);
            edge->SetTwinHalfEdge(twin_edge);
        }
        vertex_list[i]->SetOutgoingEdge(edge);
        edge->SetElement(this);
    }

    //Linking edges
    for (unsigned int i=0; i<n_nodes; ++i)
    {
        DCELHalfEdge<SPACE_DIM>* out_edge = vertex_list[i]->GetOutgoingEdge();
        const unsigned int next_index = (i+1)%n_nodes;
        const unsigned int previous_index = (i-1+n_nodes)%n_nodes;
        DCELHalfEdge<SPACE_DIM>* next_edge = vertex_list[next_index]->GetOutgoingEdge();
        DCELHalfEdge<SPACE_DIM>* previous_edge = vertex_list[previous_index]->GetOutgoingEdge();
        out_edge->SetNextHalfEdge(next_edge);
        out_edge->SetPreviousHalfEdge(previous_edge);

/*        out_edge->GetTwinHalfEdge()->SetPreviousHalfEdge(next_edge->GetTwinHalfEdge());
        out_edge->GetTwinHalfEdge()->SetNextHalfEdge(previous_edge->GetTwinHalfEdge());*/
    }
    mHalfEdge = vertex_list[0]->GetOutgoingEdge();
}

template<unsigned int SPACE_DIM>
DCELElement<SPACE_DIM>::DCELElement(std::vector<DCELHalfEdge<SPACE_DIM>* > edge_list)
{
    mHalfEdge = edge_list[0];
    for (DCELHalfEdge<SPACE_DIM>* edge:edge_list)
        edge->SetElement(this);
}

template<unsigned int SPACE_DIM>
DCELElement<SPACE_DIM>::~DCELElement()
{}

template<unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>* DCELElement<SPACE_DIM>::GetHalfEdge() const
{
    return mHalfEdge;
}

template<unsigned int SPACE_DIM>
void DCELElement<SPACE_DIM>::SetHalfEdge(DCELHalfEdge<SPACE_DIM>* edge)
{
    mHalfEdge = edge;
}

template class DCELElement<1>;
template class DCELElement<2>;
template class DCELElement<3>;

