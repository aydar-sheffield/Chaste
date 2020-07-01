/*
 * HEMutableVertexMesh.hpp
 *
 *  Created on: 30 Mar 2020
 *      Author: aydar
 */

#ifndef HEMUTABLEVERTEXMESH_HPP_
#define HEMUTABLEVERTEXMESH_HPP_

#include "HEVertexMesh.hpp"
#include "AbstractMutableVertexMesh.hpp"
#include "MutableVertexMesh.hpp"

template <unsigned int SPACE_DIM>
class HEMutableVertexMesh: public HEVertexMesh<SPACE_DIM>, public AbstractMutableVertexMesh<SPACE_DIM, SPACE_DIM>
{
    friend class TestHEMutableVertexMesh;
    friend class TestHEMutableVertexMeshReMesh;
protected:
    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    std::vector<HalfEdge<SPACE_DIM>* > mDeletedHalfEdges;

    /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;

    /**
     * Locations of T1 swaps (the mid point of the moving nodes), stored so they can be accessed and output by the cell population.
     * The locations are stored until they are cleared by ClearLocationsOfT1Swaps().
     */
    std::vector< c_vector<double, SPACE_DIM> > mLocationsOfT1Swaps;

    /**
     * The location of the last T2 swap (the centre of the removed triangle), stored so it can be accessed by the T2SwapCellKiller.
     */
    c_vector<double, SPACE_DIM> mLastT2SwapLocation;

    /**
     * Locations of T3 swaps (the location of the intersection with the edge), stored so they can be accessed and output by the cell population.
     * The locations are stored until they are cleared by ClearLocationsOfT3Swaps().
     */
    std::vector< c_vector<double, SPACE_DIM> > mLocationsOfT3Swaps;

    /**
     * Divide an element along the axis passing through the origin nodes of two halfedges.
     *
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     *
     * @param pElement the element to divide
     * @param pEdgeA halfedge of one of the nodes in this element
     * @param pEdgeB halfedge of another node in this element
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElement(HEElement<SPACE_DIM>* pElement,
                           HalfEdge<SPACE_DIM>* pEdgeA,
                           HalfEdge<SPACE_DIM>* pEdgeB,
                           bool placeOriginalElementBelow=false);

    /**
     * Helper method for ReMesh().
     *
     * Check if any neighbouring nodes in an element are closer than the mCellRearrangementThreshold
     * and are not contained in any triangular elements. If any such pair of nodes are found, then
     * call IdentifySwapType(), which in turn implements the appropriate local remeshing operation
     * (a T1 swap, void removal, or node merge).
     *
     * @return whether we need to check for, and implement, any further local remeshing operations
     *                   (true if any swaps are performed).
     */
    bool CheckForSwapsFromShortEdges();

    /**
     * Helper method for ReMesh().
     *
     * Check if any elements have become intersected and correct this by implementing the appropriate
     * local remeshing operation (a T3 swap or node merge).
     *
     * @return whether to recheck the mesh again
     */
    bool CheckForIntersections();

    /**
     * Helper method for ReMesh(), called by CheckForSwapsFromShortEdges() when
     * neighbouring nodes in an element have been found to be closer than the mCellRearrangementThreshold
     * and do not share any triangular elements.
     *
     * Identify the type of local remeshing operation required (T1 swap, void removal, or node merge).
     * @param full_edge edge to be swapped
     */
    void IdentifySwapType(FullEdge<SPACE_DIM>& full_edge);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Merge two given nodes in the mesh and update node/element ownership, by replacing
     * the node contained in the least number of elements with the other node. The merged
     * node is moved to the centre between the two old node positions.
     *
     * @param pNodeA one of the nodes to perform the merge with
     * @param pNodeB the other node to perform the merge with
     */
    void PerformNodeMerge(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB);

    /**
     * Collapses an edge by merging its two end nodes. By default, the new node is in the middle of the collapsed edge
     * @param pEdge collapsing edge
     */
    void CollapseEdge(HalfEdge<SPACE_DIM>* pEdge);

    /**
     * Helper method for merging edges during T3 swap, such that edge_A is the internal edge of the intersected element
     * and edge_p is the edge pointing to the node that is close to the intersected element.
     * The new edge in the intersected element will be common->pNode->edge_A target.
     * @param edge_A
     * @param edge_p
     * @param pNode
     */
    void MergeEdgesInT3Swap(HalfEdge<SPACE_DIM>* edge_A, HalfEdge<SPACE_DIM>* edge_p, HENode<SPACE_DIM>* pNode);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Perform a T1 swap on two given nodes contained in a given set of elements.
     * This involves replacing the two nodes with two new nodes placed on either
     * side of the previous shared edge, such that the edge formed by the two new nodes
     * is the perpendicular bisector of the previous shared edge, and 'just larger' (by a
     * factor mCellRearrangementRatio) than mThresholdDistance.
     *
     * @param pNodeA one of the nodes to perform the swap with
     * @param pNodeB the other node to perform the swap
     * @param rElementsContainingNodes set of common elements
     */
    void PerformT1Swap(HalfEdge<SPACE_DIM>* pEdge);

    /**
     * Helper method for ReMesh(), called by CheckForIntersections().
     *
     * Perform an element swap to resolve the situation where a given node
     * has been found to overlap a given element not containing it.
     *
     * @param pNode pointer to the node
     * @param elementIndex global index of the element in the mesh
     */
    void PerformIntersectionSwap(HENode<SPACE_DIM>* pNode, unsigned elementIndex);

    /**
     * Helper method for ReMesh(), called by CheckForT2Swaps().
     *
     * Perform a T2 swap on a given triangular element whose area is smaller than mT2Threshold
     * by replacing it with a new node located at the centroid of the former element and updating
     * node/element ownership.
     *
     * @param rElement the element to remove
     */
    void PerformT2Swap(HEElement<SPACE_DIM>& rElement);

    /**
     * Helper method for ReMesh(), called by CheckForIntersections().
     *
     * Perform a T3 swap on a given node that has been found to overlap a given element
     * by moving the node back onto the edge of that element, associating it with the
     * element and adding new nodes to maintain three elements per node.
     *
     * @param pNode pointer to the node
     * @param elementIndex global index of the element in the mesh
     */
    void PerformT3Swap(HENode<SPACE_DIM>* pNode, unsigned elementIndex);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Remove a triangular void bounded by three connected halfedges, in which one of the edges is
     * less than mCellRearrangementThreshold, through calls to PerformNodeMerge().
     * Note that all halfedges of the void must be accessible (i.e. linked)
     * @param pEdge arbitrary edge of the void
     */
    void PerformVoidRemoval(HalfEdge<SPACE_DIM>* pEdge);

    /**
     * Helper method for PerformT2Swap and PerformVoidRemoval.
     * Removes a triangle outlined by pEdge--pEdge->GetNextHalfEdge()--pEdge->GetPreviousHalfEdge()
     * Note that the triangle can either correspond to an element (i.e. we are performing T2 swap)
     * or to a void (i.e. we are performing void removal).
     * Note that pEdge->GetNextHalfEdge()->GetNextHalfEdge() == pEdge, i.e. correct halfedge must be suppled
     * @param pEdge arbitrary inner edge of the triangle outlined by its three (half)edges
     * @return the node at the intersection of the three edges, replacing the triangle
     */
    HENode<SPACE_DIM>* RemoveTriangle(HalfEdge<SPACE_DIM>* pEdge);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Handles the case where a swap involves a junction with more than three adjacent elements.
     * This is implemented in a separate method to allow child classes to override this behaviour
     * and implement junction remodelling with high-order nodes (see #2664).
     * @param full_edge the junction involved in the swap
     */
    virtual void HandleHighOrderJunctions(FullEdge<SPACE_DIM>& full_edge);

    /**
     * Helper method for ReMesh(), called by HandleHighOrderJunctions().
     *
     * Merge a node of a non-rosette cell with the central node
     * of a rosette, by replacing the node in the non-rosette
     * cell with that of the rosette centre, keeping the rosette
     * centre in the same position.
     *
     * @param pNodeA one of the nodes to perform the merge with
     * @param pNodeB the other node to perform the merge with
     */
    void PerformRosetteRankIncrease(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh(), called by CheckForRosettes().
     *
     * Split protorosette in random direction.
     * Create new node and redistribute nodes along line joining
     * centres of randomly selected cell and cell opposite it.
     *
     * @param pProtorosetteNode node at centre of protorosette
     */
    void PerformProtorosetteResolution(HENode<SPACE_DIM>* pProtorosetteNode);

    /**
     * Helper method for ReMesh(), called by CheckForRosettes().
     *
     * Split rosette by removing one cell at random.
     * Create new node positioned along line joining
     * rosette node and centre of randomly selected cell.
     * Rosette node will remain unmoved.
     *
     * @param pRosetteNode node at centre of rosette
     */
    void PerformRosetteRankDecrease(HENode<SPACE_DIM>* pRosetteNode);

    /**
     * Helper method for ReMesh().
     *
     * Check whether the mesh contains rosettes or protorosettes, and implement resolution events
     * if necessary.
     */
    void CheckForRosettes();

    /**
     * Helper method for ReMesh(), called by PerformT3Swap(). During T3 swaps nodes are merged onto edges. This
     * method checks if the edge is too short and moves its vertices apart if necessary in order to prevent T1 swaps
     * from happening right away. The method also checks that the location where the merged node is going to end up at
     * is not too close to one of the neighbouring vertices and moves it if necessary to prevent T1 swaps.
     *
     * @param pEdge the intersected edge
     * @param intersection  the intersection location, i.e. the location where we are planning to put the merged node
     *
     * @return intersection, the corrected location of where we are planning to put the merged node
     */
    c_vector<double, 2> WidenEdgeOrCorrectIntersectionLocationIfNecessary(HalfEdge<SPACE_DIM>* pEdge, c_vector<double,2> intersection);

    /**
     * Helper method to mark halfedge for deletion
     */
    void MarkHalfEdgeAsDeleted(HalfEdge<SPACE_DIM>* half_edge);
public:

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param cellRearrangementRatio ratio between the minimum threshold distance for element
     *                                rearrangement node separation after remeshing (defaults to 1.5)
     * @param protorosetteFormationProbability the probability of a protorosette formation event happening instead of
     *                                a T1 swap
     * @param protorosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a protorosette
     *                                will resolve (similar to the completion of a T1 swap)
     * @param rosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a rosette will
     *                                resolve (reduce the number of cells sharing a common vertex by 1)
     */
    HEMutableVertexMesh(std::vector<HENode<SPACE_DIM>*> nodes,
                      std::vector<HEElement< SPACE_DIM>*> vertexElements,
                      double cellRearrangementThreshold=0.01,
                      double t2Threshold=0.001,
                      double cellRearrangementRatio=1.5,
                      double protorosetteFormationProbability=0.0,
                      double protorosetteResolutionProbabilityPerTimestep=0.0,
                      double rosetteResolutionProbabilityPerTimestep=0.0);

    /**
     * Default constructor for use by serializer.
     */
    HEMutableVertexMesh();

    /**
     * Destructor.
     */
    virtual ~HEMutableVertexMesh();

    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* ConvertToMutableVertexMesh() const;

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of full edges
     */
    virtual unsigned GetNumFullEdges() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the locations of the T1 swaps
     */
    std::vector< c_vector<double, SPACE_DIM> > GetLocationsOfT1Swaps();

    /**
     * Helper method to clear the stored T1 swaps
     */
    void ClearLocationsOfT1Swaps();

    /**
     * @return the location of the last T2 swap
     */
    c_vector<double, SPACE_DIM> GetLastT2SwapLocation();

    /**
     * @return the locations of the T3 swaps
     */
    std::vector< c_vector<double, SPACE_DIM> > GetLocationsOfT3Swaps();

    /**
     * Helper method to clear the stored T3 swaps
     */
    void ClearLocationsOfT3Swaps();

    /**
     * Add a node to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     *
     * @param pNewNode pointer to the new node
     * @return the global index of the new node in the mesh.
     */
    unsigned AddNode(HENode<SPACE_DIM>* pNewNode);

    /**
     * Add an element to the mesh.
     *
     * @param pNewElement the new element
     *
     * @return the index of the new element in the mesh
     */
    unsigned AddElement(HEElement<SPACE_DIM>* pNewElement);

    /**
     * Bundle pEdge with its twin into FullEdge and add it to the mesh.
     * @param pEdge
     * @return new FullEdge index
     */
    unsigned int AddEdge(HalfEdge<SPACE_DIM>* pEdge);

    /**
     * Mark an element as deleted. Note that it DOES NOT deal with the associated
     * nodes and therefore should only be called immediately prior to a ReMesh()
     * being called.
     *
     * @param index  the global index of a specified vertex element
     */
    void DeleteElementPriorToReMesh(unsigned index);

    /**
     * Mark a given node as deleted. Note that this method DOES NOT deal with the
     * associated elements and therefore should only be called immediately prior
     * to a ReMesh() being called.
     *
     * @param index The index of the node to delete
     */
    void DeleteNodePriorToReMesh(unsigned index);

    /**
     * Divide an element along its short axis.
     *
     * @param pElement the element to divide
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongShortAxis(HEElement<SPACE_DIM>* pElement,
                                         bool placeOriginalElementBelow=false);

    /**
     * Divide an element along a specified axis.
     *
     * If the new nodes (intersections of axis with element) are within
     * mCellRearrangementThreshold of existing nodes then they are
     * moved 2*mCellRearrangementThreshold away.
     *
     * @param pElement the element to divide
     * @param axisOfDivision axis to divide the element by
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongGivenAxis(HEElement<SPACE_DIM>* pElement,
                                         c_vector<double, SPACE_DIM> axisOfDivision,
                                         bool placeOriginalElementBelow=false);

    /**
     * Helper method for ReMesh().
     *
     * Check for any triangular element whose area is smaller than mT2Threshold
     * and call PerformT2Swap() on any such element.
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     *
     * @return whether we need to check for, and implement, any further local remeshing operations
     */
    bool CheckForT2Swaps(VertexElementMap& rElementMap);

    /**
     * Delete mNodes and mElements.
     */
    void Clear();

    /**
     * Add a node on the edge between two nodes.
     *
     * @param pNodeA a pointer to one node
     * @param pNodeB a pointer to the other nodes
     */
    void DivideEdge(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh(). Removes the deleted nodes and elements from the mesh and updates the
     * rElementMap accordingly.
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    void RemoveDeletedNodesAndElements(VertexElementMap& rElementMap);

    /**
     * Helper method for ReMesh(). Removes the deleted nodes from the mesh and relabels the node indices.
     */
    void RemoveDeletedNodes();

    void RemoveDeletedEdges();

    /**
     * Update the state of the mesh by implementing any local remeshing operations (node merging,
     * or T1, T2 or T3 swaps) that are required, and store any changes in element indices using
     * the given VertexElementMap.
     *
     * This method calls several other methods, in particular CheckForT2Swaps(), CheckForSwapsFromShortEdges()
     * and CheckForIntersections().
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    virtual void ReMesh(VertexElementMap& rElementMap);

    /**
     * Alternative version of ReMesh which takes no parameters and does not require a VertexElementMap.
     * Note: inherited classes should overload ReMesh(VertexElementMap&).
     *
     * \todo This method seems to be redundant; remove it? (#2401)
     */
    void ReMesh();

    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);
};

#endif /* HEMUTABLEVERTEXMESH_HPP_ */
