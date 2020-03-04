/*
 * DCELElement.h
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef DCELELEMENT_HPP_
#define DCELELEMENT_HPP_
#include "DCELHalfEdge.hpp"
#include "DCELVertex.hpp"
template<unsigned int SPACE_DIM>
class DCELElement
{
public:
    DCELElement();
    //Vertex list given in CCW order
    DCELElement(std::vector<DCELVertex<SPACE_DIM>* > vertex_list);
    //Edge list given in CCW order
    DCELElement(std::vector<DCELHalfEdge<SPACE_DIM>* > edge_list);
    ~DCELElement();
    DCELHalfEdge<SPACE_DIM>* GetHalfEdge() const;
    void SetHalfEdge(DCELHalfEdge<SPACE_DIM>* edge);

private:
    DCELHalfEdge<SPACE_DIM>* mHalfEdge;
};

#endif /* _DCELELEMENT_HPP_ */
