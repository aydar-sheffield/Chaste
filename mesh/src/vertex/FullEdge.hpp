/*
 * FullEdge.hpp
 *
 *  Created on: 25 May 2020
 *      Author: aydar
 */

#ifndef FULLEDGE_HPP_
#define FULLEDGE_HPP_

#include "HalfEdge.hpp"
#include <set>
/**
 * Container class to represent a (full) edge.
 * Useful for edge traversal over the mesh, for example
 */
template<unsigned int SPACE_DIM>
class FullEdge
{
private:
    HalfEdge<SPACE_DIM>* mpFirst;
    HalfEdge<SPACE_DIM>* mpSecond;
    unsigned int mIndex;
public:
    FullEdge(const unsigned index = 0);
    FullEdge(HalfEdge<SPACE_DIM>* edge0, HalfEdge<SPACE_DIM>* edge1, const unsigned index = 0);
    FullEdge(HalfEdge<SPACE_DIM>* edge0, const unsigned index = 0);
    ~FullEdge();

    double GetLength() const;
    double ComputeLength();
    unsigned int GetIndex() const;
    void SetIndex(const unsigned int new_index);
    std::set<unsigned int> GetContainingElementIndices() const;
    HalfEdge<SPACE_DIM>* operator()(unsigned int index);
    HalfEdge<SPACE_DIM>* at(unsigned int index) const;
};

#endif /* FULLEDGE_HPP_ */
