/*
 * DCELEdge.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "DCELHalfEdge.hpp"

#include "DCELElement.hpp"
#include "DCELVertex.hpp"

template<unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>::DCELHalfEdge()
{
    DCELHalfEdge<SPACE_DIM>* mTwin = nullptr;
    DCELHalfEdge<SPACE_DIM>* mNextEdge= nullptr;
    DCELHalfEdge<SPACE_DIM>* mPreviousEdge= nullptr;
    DCELVertex<SPACE_DIM>* mTargetVertex= nullptr;
    DCELElement<SPACE_DIM>* mElement = nullptr;
}

template<unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>::~DCELHalfEdge()
{}

template<unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>* DCELHalfEdge<SPACE_DIM>::GetTwinHalfEdge() const
{
    return mTwin;
}

template<unsigned int SPACE_DIM>
void DCELHalfEdge<SPACE_DIM>::SetTwinHalfEdge(DCELHalfEdge<SPACE_DIM>* edge)
{
    mTwin = edge;
}

template<unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>* DCELHalfEdge<SPACE_DIM>::GetNextHalfEdge() const
{
    return mNextEdge;
}

template<unsigned int SPACE_DIM>
void DCELHalfEdge<SPACE_DIM>::SetNextHalfEdge(DCELHalfEdge<SPACE_DIM>* edge)
{
    mNextEdge = edge;
}

template<unsigned int SPACE_DIM>
DCELHalfEdge<SPACE_DIM>* DCELHalfEdge<SPACE_DIM>::GetPreviousHalfEdge() const
{
    return mPreviousEdge;
}

template<unsigned int SPACE_DIM>
void DCELHalfEdge<SPACE_DIM>::SetPreviousHalfEdge(DCELHalfEdge<SPACE_DIM>* edge)
{
    mPreviousEdge = edge;
}

template<unsigned int SPACE_DIM>
DCELVertex<SPACE_DIM>* DCELHalfEdge<SPACE_DIM>::GetOriginVertex() const
{
    return mTwin->GetTargetVertex();
}

template<unsigned int SPACE_DIM>
void DCELHalfEdge<SPACE_DIM>::SetOriginVertex(DCELVertex<SPACE_DIM>* vertex)
{
    //nothing yet
}

template<unsigned int SPACE_DIM>
DCELVertex<SPACE_DIM>* DCELHalfEdge<SPACE_DIM>::GetTargetVertex() const
{
    return mTargetVertex;
}

template<unsigned int SPACE_DIM>
void DCELHalfEdge<SPACE_DIM>::SetTargetVertex(DCELVertex<SPACE_DIM>* vertex)
{
    mTargetVertex = vertex;
}

template<unsigned int SPACE_DIM>
DCELElement<SPACE_DIM>* DCELHalfEdge<SPACE_DIM>::GetElement() const
{
    return mElement;
}

template<unsigned int SPACE_DIM>
void DCELHalfEdge<SPACE_DIM>::SetElement(DCELElement<SPACE_DIM>* element)
{
    mElement = element;
}

