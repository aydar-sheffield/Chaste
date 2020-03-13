/*
 * DCELVertex.hpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef HEVERTEX_HPP_
#define HEVERTEX_HPP_

#include "HalfEdge.hpp"
#include "UblasVectorInclude.hpp"
#include "Node.hpp"
template<unsigned int SPACE_DIM>
class HEVertex : public Node<SPACE_DIM>
{
public:
    /**
     * Constructor that takes in the node's location as a std::vector.
     *
     * @param index  the index of the node in the mesh
     * @param coords  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param edge the pointer to the outgoing edge
     */
    HEVertex(unsigned index, std::vector<double> coords,
             bool isBoundaryNode=false, HalfEdge<SPACE_DIM>* edge = nullptr);

    /**
     * Constructor that takes in the node's location as a c_vector.
     *
     * @param index  the index of the node in the mesh
     * @param location  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param edge the pointer to the outgoing edge
     */
    HEVertex(unsigned index, c_vector<double, SPACE_DIM> location,
             bool isBoundaryNode=false, HalfEdge<SPACE_DIM>* edge = nullptr);

    /**
     * Constructor that takes the coordinates of the node's location as separate input arguments.
     *
     * @param index  the index of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param v1 the x-coordinate of the node in the mesh (defaults to 0)
     * @param v2 the y-coordinate of the node in the mesh (defaults to 0)
     * @param v3 the z-coordinate of the node in the mesh (defaults to 0)
     * @param edge the pointer to the outgoing edge
     */
    HEVertex(unsigned index, bool isBoundaryNode=false,
             double v1=0, double v2=0, double v3=0,
             HalfEdge<SPACE_DIM>* edge = nullptr);
    /**
     * Constructor to convert Node object ot HEVertex object
     * @param node
     */
    HEVertex(const Node<SPACE_DIM> &node, HalfEdge<SPACE_DIM>* edge = nullptr);

    ~HEVertex();
    HalfEdge<SPACE_DIM>* GetOutgoingEdge() const override;
    void SetOutgoingEdge(HalfEdge<SPACE_DIM>* edge) override;

    HalfEdge<SPACE_DIM>* GetIncomingEdge() const;
    /** Get vertex the outgoing edge points to **/
    HEVertex<SPACE_DIM>* GetNextVertex() const;

    void UpdateElementIndices();
private:
    /** Outgoing edge **/
    HalfEdge<SPACE_DIM>* mEdge;
};


#endif /* DCELVERTEX_HPP_ */
