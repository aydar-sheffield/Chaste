/*
 * DCELHalfEdge.hpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef DCELHALFEDGE_HPP_
#define DCELHALFEDGE_HPP_

template<unsigned int SPACE_DIM>
class DCELVertex;
template<unsigned int SPACE_DIM>
class DCELElement;

template<unsigned int SPACE_DIM>
class DCELHalfEdge
{
public:
    DCELHalfEdge();
    ~DCELHalfEdge();
    DCELHalfEdge<SPACE_DIM>* GetTwinHalfEdge() const;
    void SetTwinHalfEdge(DCELHalfEdge<SPACE_DIM>* edge);

    DCELHalfEdge<SPACE_DIM>* GetNextHalfEdge() const;
    void SetNextHalfEdge(DCELHalfEdge<SPACE_DIM>* edge);

    DCELHalfEdge<SPACE_DIM>* GetPreviousHalfEdge() const;
    void SetPreviousHalfEdge(DCELHalfEdge<SPACE_DIM>* edge);

    DCELVertex<SPACE_DIM>* GetOriginVertex() const;
    void SetOriginVertex(DCELVertex<SPACE_DIM>* vertex);

    DCELVertex<SPACE_DIM>* GetTargetVertex() const;
    void SetTargetVertex(DCELVertex<SPACE_DIM>* vertex);

    DCELElement<SPACE_DIM>* GetElement() const;
    void SetElement(DCELElement<SPACE_DIM>* element);
private:
    DCELHalfEdge<SPACE_DIM>* mTwin;
    DCELHalfEdge<SPACE_DIM>* mNextEdge;
    DCELHalfEdge<SPACE_DIM>* mPreviousEdge;
    DCELVertex<SPACE_DIM>* mTargetVertex;
    DCELElement<SPACE_DIM>* mElement;
};

#endif /* DCELHALFEDGE_HPP_ */
