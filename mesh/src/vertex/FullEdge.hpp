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
    double mLength;
public:
    FullEdge();
    FullEdge(HalfEdge<SPACE_DIM>* edge0, HalfEdge<SPACE_DIM>* edge1);
    FullEdge(HalfEdge<SPACE_DIM>* edge0);
    ~FullEdge();

    double GetLength() const;
    double ComputeLength();
    std::set<unsigned int> GetContainingElementIndices() const;
    HalfEdge<SPACE_DIM>* operator()(unsigned int index);
};

#endif /* FULLEDGE_HPP_ */
