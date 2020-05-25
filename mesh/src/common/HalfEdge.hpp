/*
 * DCELHalfEdge.hpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef HALFEDGE_HPP_
#define HALFEDGE_HPP_

template<unsigned int SPACE_DIM>
class HENode;
template<unsigned int SPACE_DIM>
class HEElement;

template<unsigned int SPACE_DIM>
class HalfEdge
{
private:
    HalfEdge<SPACE_DIM>* mpTwin;
    HalfEdge<SPACE_DIM>* mpNextEdge;
    HalfEdge<SPACE_DIM>* mpPreviousEdge;
    HENode<SPACE_DIM>* mpTargetNode;
    HEElement<SPACE_DIM>* mpElement;

    unsigned int mIndex;
    bool mIsDeleted;
    double mLength;
public:
    HalfEdge(HEElement<SPACE_DIM>* pElement = nullptr);
    ~HalfEdge();

    /**
     * @return twin edge
     */
    HalfEdge<SPACE_DIM>* GetTwinHalfEdge() const;

    /**
     * Set twin halfedge
     * @param edge
     * @param SetOtherEdgeToo if true, sets *this* as twin of half edge
     */
    void SetTwinHalfEdge(HalfEdge<SPACE_DIM>* pEdge, const bool SetOtherEdgeToo = false);

    /**
     * @return next halfedge
     */
    HalfEdge<SPACE_DIM>* GetNextHalfEdge() const;

    /**
     * Set next halfedge
     * @param edge
     * @param SetOtherEdgeToo if true, sets *this* as previous halfedge of *edge*
     */
    void SetNextHalfEdge(HalfEdge<SPACE_DIM>* pEdge, const bool SetOtherEdgeToo = false);

    /**
     * Get previous halfedge
     * @return
     */
    HalfEdge<SPACE_DIM>* GetPreviousHalfEdge() const;

    /**
     * Set previous halfedge
     * @param edge
     * @param SetOtherEdgeToo if true, sets *this* as next halfedge of *edge*
     */
    void SetPreviousHalfEdge(HalfEdge<SPACE_DIM>* pEdge, const bool SetOtherEdgeToo = false);

    /**
     * Gets origin Node of this halfedge.
     * @return
     */
    HENode<SPACE_DIM>* GetOriginNode() const;

    /**
     * Set origin Node of this halfedge.
     * This method modifies target vertices of halfedges adjacent to this edge
     * @param Node
     */
    void SetOriginNode(HENode<SPACE_DIM>* pNode);

    /**
     * Gets target Node of this halfedge
     * @return
     */
    HENode<SPACE_DIM>* GetTargetNode() const;

    /**
     * Sets target Node of this edge
     * @param Node
     * @param ModifyAdjacentEdges if true, sets *Node* as target Node of next halfedge's twin
     */
    void SetTargetNode(HENode<SPACE_DIM>* pNode, const bool ModifyAdjacentEdges = false);

    /**
     * Get the element to which this halfedge belongs to
     * @return
     */
    HEElement<SPACE_DIM>* GetElement() const;

    /**
     * Set the element to which this halfedge belongs to
     * @param element
     */
    void SetElement(HEElement<SPACE_DIM>* pElement);

    /**
     * Gets global index of this halfedge
     * @return
     */
    unsigned int GetIndex() const;

    /**
     * Sets global index of this halfedge
     * @param new_index
     */
    void SetIndex(const unsigned int new_index);

    bool IsDeleted() const;

    void SetDeletedStatus(const bool status, const bool SetTwin = false);

    bool IsFullyInitialized() const;

    bool IsOnBoundary() const;

    double ComputeLength() const;

    inline double GetLength();

    void UpdateLength(const double new_length = 0);
};

#endif /*HALFEDGE_HPP_ */
