/*
 * DCELVertex.hpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef DCELVERTEX_HPP_
#define DCELVERTEX_HPP_

#include "DCELHalfEdge.hpp"
#include "UblasVectorInclude.hpp"
template<unsigned int SPACE_DIM>
class DCELVertex
{
public:
    DCELVertex();
    DCELVertex(const c_vector<double, SPACE_DIM> location);
    DCELVertex(const double x, const double y, const double z = 0.0);
    ~DCELVertex();
    DCELHalfEdge<SPACE_DIM>* GetOutgoingEdge() const;
    void SetOutgoingEdge(DCELHalfEdge<SPACE_DIM>* edge);
    DCELHalfEdge<SPACE_DIM>* GetIncomingEdge() const;
    DCELVertex<SPACE_DIM>* GetNextVertex() const;
    void SetNextVertex();
    DCELVertex<SPACE_DIM>* GetPreviousVertex() const;
    void SetPreviousVertex();

    /**
     * @return the node's location as a c_vector.
     *
     * The returned location may not be modified; if you want that functionality use
     * rGetModifiableLocation instead.
     */
    const c_vector<double, SPACE_DIM>& rGetLocation() const;

    /**
     * @return the node's location as a c_vector.
     *
     * If you modify the returned location,
     * Jacobian and JacobianDeterminant of elements need to be updated.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    c_vector<double, SPACE_DIM>& rGetModifiableLocation();
private:
    void CommonConstructor();
    /** Location **/
    c_vector<double, SPACE_DIM> mLocation;
    bool mIsDeleted;

    /** Outgoing edge **/
    DCELHalfEdge<SPACE_DIM>* mEdge;
    DCELVertex<SPACE_DIM>* mNextVertex;
    DCELVertex<SPACE_DIM>* mPreviousVertex;
};


#endif /* DCELVERTEX_HPP_ */
