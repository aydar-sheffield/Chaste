/*
 * DCELHalfEdge.hpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef HALFEDGE_HPP_
#define HALFEDGE_HPP_

template<unsigned int SPACE_DIM>
class HEVertex;
template<unsigned int SPACE_DIM>
class HEElement;

template<unsigned int SPACE_DIM>
class HalfEdge
{
public:
    HalfEdge();
    ~HalfEdge();
    HalfEdge<SPACE_DIM>* GetTwinHalfEdge() const;
    void SetTwinHalfEdge(HalfEdge<SPACE_DIM>* edge);

    HalfEdge<SPACE_DIM>* GetNextHalfEdge() const;
    void SetNextHalfEdge(HalfEdge<SPACE_DIM>* edge);

    HalfEdge<SPACE_DIM>* GetPreviousHalfEdge() const;
    void SetPreviousHalfEdge(HalfEdge<SPACE_DIM>* edge);

    HEVertex<SPACE_DIM>* GetOriginVertex() const;

    HEVertex<SPACE_DIM>* GetTargetVertex() const;
    void SetTargetVertex(HEVertex<SPACE_DIM>* vertex);

    HEElement<SPACE_DIM>* GetElement() const;
    void SetElement(HEElement<SPACE_DIM>* element);

    unsigned int GetIndex() const;
    void SetIndex(const unsigned int new_index);
private:
    HalfEdge<SPACE_DIM>* mTwin;
    HalfEdge<SPACE_DIM>* mNextEdge;
    HalfEdge<SPACE_DIM>* mPreviousEdge;
    HEVertex<SPACE_DIM>* mTargetVertex;
    HEElement<SPACE_DIM>* mElement;

    unsigned int mIndex;
};

#endif /*HALFEDGE_HPP_ */
