/*
 * HEHEMutableVertexMesh.cpp
 *
 *  Created on: 30 Mar 2020
 *      Author: aydar
 */

#include "HEMutableVertexMesh.hpp"

#include "LogFile.hpp"
#include "UblasCustomFunctions.hpp"
#include "Warnings.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned SPACE_DIM>
HEMutableVertexMesh<SPACE_DIM>::HEMutableVertexMesh(std::vector<HENode<SPACE_DIM>*> HENodes,
                                                    std::vector<HEElement<SPACE_DIM>*> HEElements,
                                                    double cellRearrangementThreshold,
                                                    double t2Threshold,
                                                    double cellRearrangementRatio,
                                                    double protorosetteFormationProbability,
                                                    double protorosetteResolutionProbabilityPerTimestep,
                                                    double rosetteResolutionProbabilityPerTimestep)
                                                    :
    HEVertexMesh<SPACE_DIM>(HENodes, HEElements),
    AbstractMutableVertexMesh<SPACE_DIM, SPACE_DIM>(cellRearrangementThreshold,
                                                    t2Threshold,
                                                    cellRearrangementRatio,
                                                    protorosetteFormationProbability,
                                                    protorosetteResolutionProbabilityPerTimestep,
                                                    rosetteResolutionProbabilityPerTimestep)
{
    // If in 3D, then also populate mFaces
    if (SPACE_DIM == 3)
    {
        EXCEPTION("3D MutableVertexMesh is not supported");
    }

    this->mMeshChangesDuringSimulation = true;
}

template<unsigned SPACE_DIM>
HEMutableVertexMesh<SPACE_DIM>::HEMutableVertexMesh()
{
    // Note that the member variables initialised above will be overwritten as soon as archiving is complete
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

template<unsigned SPACE_DIM>
HEMutableVertexMesh<SPACE_DIM>::~HEMutableVertexMesh()
{

}

template<unsigned int SPACE_DIM>
MutableVertexMesh<SPACE_DIM, SPACE_DIM>* HEMutableVertexMesh<SPACE_DIM>::ConvertToMutableVertexMesh() const
{
    NodesAndElements<SPACE_DIM> nodes_elements = this->ConvertToVertexNodesAndElements();
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* mutable_vertex_mesh = new MutableVertexMesh<SPACE_DIM, SPACE_DIM>(nodes_elements.first, nodes_elements.second);
    return mutable_vertex_mesh;
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - this->mDeletedNodeIndices.size();
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumFullEdges() const
{
    return this->mFullEdges.size() - this->mDeletedHalfEdges.size()/2;
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - this->mDeletedElementIndices.size();
}

template<unsigned SPACE_DIM>
std::vector< c_vector<double, SPACE_DIM> > HEMutableVertexMesh<SPACE_DIM>::GetLocationsOfT1Swaps()
{
    return mLocationsOfT1Swaps;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::ClearLocationsOfT1Swaps()
{
    mLocationsOfT1Swaps.clear();
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> HEMutableVertexMesh<SPACE_DIM>::GetLastT2SwapLocation()
{
    return mLastT2SwapLocation;
}

template<unsigned SPACE_DIM>
std::vector< c_vector<double, SPACE_DIM> > HEMutableVertexMesh<SPACE_DIM>::GetLocationsOfT3Swaps()
{
    return mLocationsOfT3Swaps;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::ClearLocationsOfT3Swaps()
{
    mLocationsOfT3Swaps.clear();
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::Clear()
{
    this->mDeletedNodeIndices.clear();
    this->mDeletedElementIndices.clear();
    HEVertexMesh<SPACE_DIM>::Clear();
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::AddNode(HENode<SPACE_DIM>* pNewNode)
{
    /**
     * If no node is marked as deleted, add the node to the list of nodes.
     * Else, get the index of the last node marked, and set this index to pNewNode
     * Remove the marked node from mDeletedNodeIndices and delete it.
     */
    if (this->mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(this->mNodes.size());
        this->mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = this->mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        this->mDeletedNodeIndices.pop_back();
        //TODO: take edge deletion into account too
        delete this->mNodes[index];
        this->mNodes[index] = pNewNode;
    }
    return pNewNode->GetIndex();
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::AddElement(HEElement<SPACE_DIM>* pNewElement)
{
    unsigned new_element_index = pNewElement->GetIndex();

    if (new_element_index == this->mElements.size())
    {
        this->mElements.push_back(pNewElement);
    }
    else
    {
        this->mElements[new_element_index] = pNewElement;
    }
    return pNewElement->GetIndex();
}

template<unsigned SPACE_DIM>
unsigned int HEMutableVertexMesh<SPACE_DIM>::AddEdge(HalfEdge<SPACE_DIM>* pEdge)
{
    unsigned int new_edge_index = this->mFullEdges.size();
    this->mFullEdges.push_back(new FullEdge<SPACE_DIM>(pEdge));
    std::pair<HalfEdge<SPACE_DIM>*, FullEdge<SPACE_DIM>* > pair_0(pEdge, this->mFullEdges[new_edge_index]);
    std::pair<HalfEdge<SPACE_DIM>*, FullEdge<SPACE_DIM>* > pair_1(pEdge->GetTwinHalfEdge(), this->mFullEdges[new_edge_index]);
    this->mHalfToFullEdgeMap.insert(pair_0);
    this->mHalfToFullEdgeMap.insert(pair_1);
    return new_edge_index;
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::DivideElementAlongGivenAxis(HEElement<SPACE_DIM>* pElement,
                                                                     c_vector<double, SPACE_DIM> axisOfDivision,
                                                                     bool placeOriginalElementBelow)
{
    assert(SPACE_DIM == 2);                // LCOV_EXCL_LINE

    // Get the centroid of the element
    c_vector<double, SPACE_DIM> centroid = this->GetCentroidOfElement(pElement->GetIndex());

    // Create a vector perpendicular to the axis of division
    // NOTE: that axis of division is parallel to the newly formed edge
    c_vector<double, SPACE_DIM> perp_axis;
    perp_axis(0) = -axisOfDivision(1);
    perp_axis(1) = axisOfDivision(0);

    /*
     * Find which edges the axis of division crosses by finding any node
     * that lies on the opposite side of the axis of division to its next
     * neighbour.
     */
    std::vector<HalfEdge<SPACE_DIM>* > intersecting_edges;
    HalfEdge<SPACE_DIM>* next_edge = pElement->GetHalfEdge();
    //Left (relative to perp_axis) here means next in CW direction
    bool is_current_node_on_left
    = (inner_prod(this->GetVectorFromAtoB(next_edge->GetOriginNode()->rGetLocation(), centroid),perp_axis)>=0);

    do
    {
        bool is_next_node_on_left
        = (inner_prod(this->GetVectorFromAtoB(next_edge->GetTargetNode()->rGetLocation(), centroid), perp_axis) >= 0);
        if (is_current_node_on_left != is_next_node_on_left)
        {
            intersecting_edges.push_back(next_edge);
        }
        is_current_node_on_left = is_next_node_on_left;
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != pElement->GetHalfEdge());

    // If the axis of division does not cross two edges then we cannot proceed
    if (intersecting_edges.size() != 2)
    {
        EXCEPTION("Cannot proceed with element division: the given axis of division does not cross two edges of the element");
    }

    // Find the intersections between the axis of division and the element edges
    for (unsigned i=0; i<intersecting_edges.size(); i++)
    {
        /*
         * Get pointers to the nodes forming the edge into which one new node will be inserted.
         *
         * Note that when we use the first entry of intersecting_nodes to add a node,
         * we change the local index of the second entry of intersecting_nodes in
         * pElement, so must account for this by moving one entry further on.
         */
        Node<SPACE_DIM>* p_node_A = intersecting_edges[i]->GetOriginNode();
        Node<SPACE_DIM>* p_node_B = intersecting_edges[i]->GetTargetNode();

        c_vector<double, SPACE_DIM> position_a = p_node_A->rGetLocation();
        c_vector<double, SPACE_DIM> position_b = p_node_B->rGetLocation();
        c_vector<double, SPACE_DIM> a_to_b = this->GetVectorFromAtoB(position_a, position_b);

        c_vector<double, SPACE_DIM> intersection;

        if (norm_2(a_to_b) < 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold)
        {
            WARNING("Edge is too small for normal division; putting node in the middle of a and b. There may be T1 swaps straight away.");
            ///\todo or should we move a and b apart, it may interfere with neighbouring edges? (see #1399 and #2401)
            intersection = position_a + 0.5*a_to_b;
        }
        else
        {
            // Find the location of the intersection
            double determinant = a_to_b[0]*axisOfDivision[1] - a_to_b[1]*axisOfDivision[0];

            // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
            c_vector<double, SPACE_DIM> moved_centroid;
            moved_centroid = position_a + this->GetVectorFromAtoB(position_a, centroid);

            double alpha = (moved_centroid[0]*a_to_b[1] - position_a[0]*a_to_b[1]
                            -moved_centroid[1]*a_to_b[0] + position_a[1]*a_to_b[0])/determinant;

            intersection = moved_centroid + alpha*axisOfDivision;

            /*
             * If then new node is too close to one of the edge nodes, then reposition it
             * a distance this->mCellRearrangementRatio*this->mCellRearrangementThreshold further along the edge.
             */
            c_vector<double, SPACE_DIM> a_to_intersection = this->GetVectorFromAtoB(position_a, intersection);
            if (norm_2(a_to_intersection) < this->mCellRearrangementThreshold)
            {
                intersection = position_a + this->mCellRearrangementRatio*this->mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
            }

            c_vector<double, SPACE_DIM> b_to_intersection = this->GetVectorFromAtoB(position_b, intersection);
            if (norm_2(b_to_intersection) < this->mCellRearrangementThreshold)
            {
                assert(norm_2(a_to_intersection) > this->mCellRearrangementThreshold); // to prevent moving intersection back to original position

                intersection = position_b - this->mCellRearrangementRatio*this->mCellRearrangementThreshold*a_to_b/norm_2(a_to_b);
            }
        }

        // Add a new node to the mesh at the location of the intersection
        //The new node is on boundary if the edge is
        HENode<SPACE_DIM>* added_node = new HENode<SPACE_DIM>(0, intersecting_edges[i]->IsOnBoundary(), intersection[0], intersection[1]);
        AddNode(added_node);
        HalfEdge<SPACE_DIM>* new_edge = pElement->AddNode(intersecting_edges[i], added_node);
        AddEdge(new_edge);
    }

    // Now call DivideElement() to divide the element using the new nodes
    unsigned new_element_index = DivideElement(pElement,
                                               intersecting_edges[0],
                                               intersecting_edges[1],
                                               placeOriginalElementBelow);

    return new_element_index;
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::DivideElement(HEElement<SPACE_DIM>* pElement,
                                                       HalfEdge<SPACE_DIM>* pEdgeA,
                                                       HalfEdge<SPACE_DIM>* pEdgeB,
                                                       bool placeOriginalElementBelow)
{
    assert(SPACE_DIM == 2);                // LCOV_EXCL_LINE

    HalfEdge<SPACE_DIM>* new_edge = new HalfEdge<SPACE_DIM>();
    HalfEdge<SPACE_DIM>* new_edge_twin = new HalfEdge<SPACE_DIM>();
    AddEdge(new_edge);

    new_edge->SetTwinHalfEdge(new_edge_twin, true);

    HalfEdge<SPACE_DIM>* pEdgeA_prev = pEdgeA->GetPreviousHalfEdge();
    HalfEdge<SPACE_DIM>* pEdgeB_prev = pEdgeB->GetPreviousHalfEdge();

    new_edge->SetNextHalfEdge(pEdgeA,true);
    new_edge->SetPreviousHalfEdge(pEdgeB_prev,true);
    new_edge_twin->SetNextHalfEdge(pEdgeB,true);
    new_edge_twin->SetPreviousHalfEdge(pEdgeA_prev,true);

    new_edge->SetTargetNode(pEdgeA->GetOriginNode());
    new_edge_twin->SetTargetNode(pEdgeB->GetOriginNode());

    // Get the index of the new element
    unsigned new_element_index;
    if (this->mDeletedElementIndices.empty())
    {
        new_element_index = this->mElements.size();
    }
    else
    {
        new_element_index = this->mDeletedElementIndices.back();
        this->mDeletedElementIndices.pop_back();
        delete this->mElements[new_element_index];
    }

    // Add the new element to the mesh
    HEElement<SPACE_DIM>* new_element = new HEElement<SPACE_DIM>(new_element_index, new_edge);
    AddElement(new_element);
    pElement->SetHalfEdge(new_edge_twin);

    c_vector<double, SPACE_DIM> centroid_old, centroid_new;
    centroid_old = this->GetCentroidOfElement(pElement->GetIndex());
    centroid_new = this->GetCentroidOfElement(new_element_index);
    const bool is_original_below = centroid_old(1)<centroid_new(1);

    if (placeOriginalElementBelow)
    {
        if(!is_original_below)
        {
            pElement->SetHalfEdge(new_edge);
            new_element->SetHalfEdge(new_edge_twin);
        }
    }
    else
    {
        if(is_original_below)
        {
            pElement->SetHalfEdge(new_edge);
            new_element->SetHalfEdge(new_edge_twin);
        }
    }
    new_element->RegisterWithHalfEdges();
    pElement->RegisterWithHalfEdges();

    pElement->UpdateGeometry();
    new_element->UpdateGeometry();
    return new_element_index;
}

template<unsigned SPACE_DIM>
bool HEMutableVertexMesh<SPACE_DIM>::CheckForSwapsFromShortEdges()
{
    for (FullEdge<SPACE_DIM>* full_edge : this->mFullEdges)
    {
        //Assume edge lengths have been updated
        //Check to see if the edge is short enough for a swap...
        if (full_edge->GetLength() < this->mCellRearrangementThreshold)
        {
            // ...then check if a triangular element is adjacent to this edge
            bool edge_share_triangular_element= false;
            std::set<unsigned int> shared_elements = full_edge->GetContainingElementIndices();
            for (unsigned int element_index : shared_elements)
            {
                if (this->GetElement(element_index)->GetNumNodes()<=3)
                {
                    edge_share_triangular_element = true;
                    break;
                }
            }

            // ...and if none are, then perform the required type of swap and halt the search, returning true
            if (!edge_share_triangular_element)
            {
                IdentifySwapType(*full_edge);
            }
        }
    }
    return false;
}

template<unsigned SPACE_DIM>
bool HEMutableVertexMesh<SPACE_DIM>::CheckForIntersections()
{
    // If checking for internal intersections as well as on the boundary, then check that no nodes have overlapped any elements...
    if (this->mCheckForInternalIntersections)
    {
        ///\todo Change to only loop over neighbouring elements (see #2401)
        /// TODO: Parallelize
        for (typename AbstractMesh<SPACE_DIM,SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
                node_iter != this->GetNodeIteratorEnd();
                ++node_iter)
        {
            assert(!(node_iter->IsDeleted()));
            HENode<SPACE_DIM>* node = static_cast<HENode<SPACE_DIM>* >(&(*node_iter));
            std::set<unsigned int> containing_element_ind = node->GetContainingElementIndices();
            for (typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem_iter = this->GetElementIteratorBegin();
                    elem_iter != this->GetElementIteratorEnd();
                    ++elem_iter)
            {
                unsigned elem_index = elem_iter->GetIndex();

                // Check that the node is not part of this element
                if (containing_element_ind.count(elem_index) == 0)
                {
                    if (this->ElementIncludesPoint(node_iter->rGetLocation(), elem_index))
                    {
                        PerformIntersectionSwap(node, elem_index);
                        return true;
                    }
                }
            }
        }
    }
    else
    {
        // ...otherwise, just check that no boundary nodes have overlapped any boundary elements
        // First: find all boundary element and calculate their centroid only once
        std::vector<unsigned> boundary_element_indices;
        std::vector< c_vector<double, SPACE_DIM> > boundary_element_centroids;
        for (typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem_iter = this->GetElementIteratorBegin();
                elem_iter != this->GetElementIteratorEnd();
                ++elem_iter)
        {
            if (elem_iter->IsElementOnBoundary())
            {
                unsigned element_index = elem_iter->GetIndex();
                boundary_element_indices.push_back(element_index);
                // should be a map but I am too lazy to look up the syntax
                boundary_element_centroids.push_back(this->GetCentroidOfElement(element_index));
            }
        }
        // Second: Check intersections only for those nodes and elements within
        // this->mDistanceForT3SwapChecking within each other (node<-->element centroid)
        for (typename AbstractMesh<SPACE_DIM,SPACE_DIM>::NodeIterator node_iter = this->GetNodeIteratorBegin();
                node_iter != this->GetNodeIteratorEnd();
                ++node_iter)
        {
            if (node_iter->IsBoundaryNode())
            {
                assert(!(node_iter->IsDeleted()));
                HENode<SPACE_DIM>* node = static_cast<HENode<SPACE_DIM>* >(&(*node_iter));
                std::set<unsigned int> containing_element_ind = node->GetContainingElementIndices();

                // index in boundary_element_centroids and boundary_element_indices
                unsigned boundary_element_index = 0;
                for (std::vector<unsigned>::iterator elem_iter = boundary_element_indices.begin();
                        elem_iter != boundary_element_indices.end();
                        ++elem_iter)
                {
                    // Check that the node is not part of this element
                    if (containing_element_ind.count(*elem_iter) == 0)
                    {
                        c_vector<double, SPACE_DIM> node_location = node_iter->rGetLocation();
                        c_vector<double, SPACE_DIM> element_centroid = boundary_element_centroids[boundary_element_index];
                        double node_element_distance = norm_2(this->GetVectorFromAtoB(node_location, element_centroid));

                        if ( node_element_distance < this->mDistanceForT3SwapChecking )
                        {
                            if (this->ElementIncludesPoint(node_iter->rGetLocation(), *elem_iter))
                            {
                                this->PerformT3Swap(node, *elem_iter);

                                return true;
                            }
                        }
                    }
                    // increment the boundary element index
                    boundary_element_index +=1u;
                }
            }
        }
    }
    return false;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::IdentifySwapType(FullEdge<SPACE_DIM>& full_edge)
{
    HENode<SPACE_DIM>* pNodeA = full_edge(0)->GetTargetNode();
    HENode<SPACE_DIM>* pNodeB = full_edge(1)->GetTargetNode();
    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->GetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->GetContainingElementIndices();

    // Form the set union
    std::set<unsigned> all_indices, temp_union_set;
    std::set_union(nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                   nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                   std::inserter(temp_union_set, temp_union_set.begin()));
    all_indices.swap(temp_union_set); // temp_set will be deleted, all_indices now contains all the indices of elements
    // that touch the potentially swapping nodes

    if ((nodeA_elem_indices.size()>3) || (nodeB_elem_indices.size()>3))
    {
        /*
         * Looks like
         *
         *  \
         *   \ A   B
         * ---o---o---
         *   /
         *  /
         *
         */

        /*
         * This case is handled in a separate method to allow child classes to implement different
         * functionality for high-order-junction remodelling events (see #2664).
         */
        this->HandleHighOrderJunctions(full_edge);
    }
    else // each node is contained in at most three elements
    {
        switch (all_indices.size())
        {
        case 1:
        {
            /*
             * Each node is contained in a single element, so the nodes must lie on the boundary
             * of the mesh, as shown below. In this case, we merge the nodes and tidy up node
             * indices through calls to PerformNodeMerge() and RemoveDeletedNodes().
             *
             *    A   B
             * ---o---o---
             */
            assert(pNodeA->IsBoundaryNode());
            assert(pNodeB->IsBoundaryNode());

            CollapseEdge(full_edge(0));
            RemoveDeletedNodes();
            RemoveDeletedEdges();
            break;
        }
        case 2:
        {
            if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
            {
                if (pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode())
                {
                    /*
                     * The node configuration is as shown below, with voids on either side. In this case
                     * we perform a T1 swap, which separates the elements.
                     *
                     *   \   /
                     *    \ / Node A
                     * (1) |   (2)      (element number in brackets)
                     *    / \ Node B
                     *   /   \
                     */
                    PerformT1Swap(full_edge(0));
                }
                else if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
                {
                    /*
                     * The node configuration is as shown below, with a void on one side. We should not
                     * be able to reach this case at present, since we allow only for three-way junctions
                     * or boundaries, so we throw an exception.
                     *
                     *   \   /
                     *    \ / Node A
                     * (1) |   (2)      (element number in brackets)
                     *     x Node B
                     *     |
                     */
                    EXCEPTION("There is a non-boundary node contained only in two elements; something has gone wrong.");
                }
                else
                {
                    /*
                     * Each node is contained in two elements, so the nodes lie on an internal edge, as shown below.
                     * We should not be able to reach this case at present, since we allow only for three-way junctions
                     * or boundaries, so we throw an exception.
                     *
                     *    A   B
                     * ---o---o---
                     */
                    EXCEPTION("There are non-boundary nodes contained only in two elements; something has gone wrong.");
                }
            }// from [if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)]
            else
            {
                /*
                 * The node configuration looks like that shown below. In this case, we merge the nodes
                 * and tidy up node indices through calls to PerformNodeMerge() and  RemoveDeletedNodes().
                 *
                 * Outside
                 *         /
                 *   --o--o (2)
                 *     (1) \
                 *
                 * ///\todo this should be a T1 swap (see #1263 and #2401)
                 * Referring to the todo: this should probably stay a node-merge. If this is a T1 swap then
                 * the single boundary node will travel from element 1 to element 2, but still remain a single node.
                 * I.e. we would not reduce the total number of nodes in this situation.
                 */
                CollapseEdge(full_edge(0));
                RemoveDeletedNodes();
                RemoveDeletedEdges();
            }
            break;
        }
        case 3:
        {
            if (nodeA_elem_indices.size()==1 || nodeB_elem_indices.size()==1)
            {
                /*
                 * One node is contained in one element and the other node is contained in three elements.
                 * We should not be able to reach this case at present, since we allow each boundary node
                 * to be contained in at most two elements, so we throw an exception.
                 *
                 *    A   B
                 *
                 *  empty   /
                 *         / (3)
                 * ---o---o-----   (element number in brackets)
                 *  (1)    \ (2)
                 *          \
                 */
                assert(pNodeA->IsBoundaryNode());
                assert(pNodeB->IsBoundaryNode());

                EXCEPTION("There is a boundary node contained in three elements something has gone wrong.");
            }
            else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
            {
                // The short edge must be at the boundary. We need to check whether this edge is
                // adjacent to a triangular void before we swap. If it is a triangular void, we perform a T2-type swap.
                // If not, then we perform a normal T1 swap. I.e. in detail we need to check whether the
                // element in nodeA_elem_indices which is not in nodeB_elem_indices contains a shared node
                // with the element in nodeB_elem_indices which is not in nodeA_elem_indices.

                HalfEdge<SPACE_DIM>* boundary_edge = full_edge(0)->GetElement() == nullptr ? full_edge(0) : full_edge(1);
                HalfEdge<SPACE_DIM>* boundary_next = boundary_edge->GetNextHalfEdge();
                HalfEdge<SPACE_DIM>* boundary_previous = boundary_edge->GetPreviousHalfEdge();
                HalfEdge<SPACE_DIM>* in_element_edge = boundary_edge->GetTwinHalfEdge();
                assert(in_element_edge->GetElement() != nullptr);

                if (boundary_next->GetTargetNode()==boundary_previous->GetOriginNode())
                {
                    /*
                     * The node configuration looks like that shown below, and both nodes must be on the boundary.
                     * In this case we remove the void through a call to PerformVoidRemoval().
                     *
                     *    A  C  B                A      B
                     *      /\                 \        /
                     *     /v \                 \  (1) /
                     * (3)o----o (1)  or     (2) o----o (3)    (element number in brackets, v is a void)
                     *   /  (2) \                 \v /
                     *  /        \                 \/
                     *                             C
                     */
                    assert(pNodeA->IsBoundaryNode());
                    assert(pNodeB->IsBoundaryNode());

                    if (in_element_edge->GetElement()->GetNumNodes()==3)
                    {
                        /**
                         * Here, the triangular element would be along the short edge. Since we
                         * are already checking in CheckForSwapsFromShortEdges() whether the element
                         * is triangular, this exception is redundant for simulations. We leave it in for
                         * clarity.
                         * ///\todo: consider removing the checking for this exception (see #2401)
                         */
                        EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                    }

                    if (boundary_next->GetTwinHalfEdge()->GetElement()->GetNumNodes()==3
                            || boundary_previous->GetTwinHalfEdge()->GetElement()->GetNumNodes()==3)
                    {
                        /**
                         * If this is true then one of the elements adjacent to the triangular void
                         * is triangular. This element will then not share the short edge that is considered
                         * for a swap. Nevertheless, it would loose an edge during the swap. We are currently
                         * not able to deal with this situation.
                         * Related to #2533 and #2401.
                         */
                        EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                    }

                    PerformVoidRemoval(boundary_edge);
                }
                else
                {
                    /*
                     * The node configuration looks like that below, and both nodes must lie on the boundary.
                     * In this case we perform a T1 swap.
                     *
                     *     A  B                  A  B
                     *   \ empty/              \      /
                     *    \    /                \(1) /
                     * (3) o--o (1)  or      (2) o--o (3)    (element number in brackets)
                     *    / (2)\                /    \
                     *   /      \              /empty \
                     */
                    assert(pNodeA->IsBoundaryNode());
                    assert(pNodeB->IsBoundaryNode());

                    PerformT1Swap(boundary_edge);
                }
            } // from else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
            else
            {
                // In this case, one node must be contained in two elements and the other in three elements.
                assert (   (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==3)
                           || (nodeA_elem_indices.size()==3 && nodeB_elem_indices.size()==2) );

                // They can't both be boundary nodes
                assert(!(pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode()));

                if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
                {
                    /*
                     * The node configuration looks like that shown below. We perform a T1 swap in this case.
                     *
                     *     A  B                      A  B
                     *   \      /                  \      /
                     *    \ (1)/                    \(1) /
                     * (3) o--o (empty)  or  (empty) o--o (3)    (element number in brackets)
                     *    / (2)\                    /(2) \
                     *   /      \                  /      \
                     */
                    PerformT1Swap(full_edge(0));
                }
                else
                {
                    /*
                     * The node configuration looks like that shown below. We should not be able to reach this case
                     * at present, since we allow only for three-way junctions or boundaries, so we throw an exception.
                     *
                     *     A  B             A  B
                     *   \                       /
                     *    \  (1)           (1)  /
                     * (3) o--o---   or  ---o--o (3)    (element number in brackets)
                     *    /  (2)           (2)  \
                     *   /                       \
                     */
                    EXCEPTION("There are non-boundary nodes contained only in two elements; something has gone wrong.");
                }
            }
            break;
        }
        case 4:
        {
            /*
             * The node configuration looks like that shown below. We perform a T1 swap in this case.
             *
             *   \(1)/
             *    \ / Node A
             * (2) |   (4)      (element number in brackets)
             *    / \ Node B
             *   /(3)\
             */

            /*
             * This case is handled in a separate method to allow child classes to implement different
             * functionality for junction remodelling events (see #2664).
             */
            if (this->mProtorosetteFormationProbability > RandomNumberGenerator::Instance()->ranf())
            {
                CollapseEdge(full_edge(0));
                RemoveDeletedNodes();
                RemoveDeletedEdges();
            }
            else
            {
                this->PerformT1Swap(full_edge(0));
            }
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;

        }
    }
}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::CollapseEdge(HalfEdge<SPACE_DIM>* pEdge)
{
    //Here we delete node A
    HENode<SPACE_DIM>* pNodeA = pEdge->GetTargetNode();
    HENode<SPACE_DIM>* pNodeB = pEdge->GetOriginNode();
    // Move node A to the mid-point
    pNodeA->rGetModifiableLocation() += 0.5 * this->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation());

    HalfEdge<SPACE_DIM>* previous_edge = pEdge->GetPreviousHalfEdge();
    HalfEdge<SPACE_DIM>* next_edge = pEdge->GetNextHalfEdge();
    HalfEdge<SPACE_DIM>* twin = pEdge->GetTwinHalfEdge();
    HalfEdge<SPACE_DIM>* twin_previous = twin->GetPreviousHalfEdge();
    HalfEdge<SPACE_DIM>* twin_next = twin->GetNextHalfEdge();

    //If the outgoing halfedges of nodes A and/or B are pEdge or its twin, respectively,
    //then set outgoing halfedges to the next outgoing edge
    HalfEdge<SPACE_DIM>* node_A_out_edge = pNodeA->GetOutgoingEdge();
    if (node_A_out_edge == twin)
    {
        pNodeA->SetOutgoingEdge(twin_next);
    }

    //Since node B will be deleted, all its incoming edges must be pointing to node A
    HalfEdge<SPACE_DIM>* node_B_in_edge = pNodeB->GetOutgoingEdge()->GetTwinHalfEdge();
    node_B_in_edge->SetTargetNode(pNodeA,true);
    assert(previous_edge->GetTargetNode()==pNodeA);
    assert(twin->GetTargetNode()==pNodeA);

    //Change adjacency relations. Traversal of element edges will not include collapsed edge
    previous_edge->SetNextHalfEdge(next_edge,true);
    twin_previous->SetNextHalfEdge(twin_next,true);

    this->mDeletedNodeIndices.push_back(pNodeB->GetIndex());
    pNodeB->MarkAsDeleted();
    this->mDeletedHalfEdges.push_back(pEdge);
    this->mDeletedHalfEdges.push_back(pEdge->GetTwinHalfEdge());
    pEdge->SetDeletedStatus(true, true);

    HEElement<SPACE_DIM>* pEdgeElement = pEdge->GetElement();
    if (pEdgeElement)
    {
        if (pEdgeElement->GetHalfEdge() == pEdge)
            pEdge->GetElement()->SetHalfEdge(next_edge);
    }
    HEElement<SPACE_DIM>* pEdgeTwinElement = twin->GetElement();
    if (pEdgeTwinElement)
    {
        if (pEdgeTwinElement->GetHalfEdge()==twin)
            twin->GetElement()->SetHalfEdge(twin_next);
    }
    //TODO:Is updating adjacent elements necessary? Volume, or node numbers, for example
    next_edge = pNodeA->GetOutgoingEdge();
    do
    {
        HEElement<SPACE_DIM>* element = next_edge->GetElement();
        if (element)
        {
            element->UpdateGeometry();
        }
        next_edge = next_edge->GetTwinHalfEdge()->GetNextHalfEdge();
    }while(next_edge != twin_next);

}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::MergeEdgesInT3Swap(HalfEdge<SPACE_DIM>* edge_A, HalfEdge<SPACE_DIM>* edge_p, HENode<SPACE_DIM>* pNode)
{


    HEElement<SPACE_DIM>* intersected_element = edge_A->GetElement();
    const bool is_edge_p_inside = edge_p->GetElement();
    HEElement<SPACE_DIM>* intersecting_element = is_edge_p_inside ? edge_p->GetElement() : edge_p->GetTwinHalfEdge()->GetElement();

    //Node A is the origin of edge_A and node_B is the target of edge_A
    //Here we merge edges and modify adjacency relations
    if (is_edge_p_inside)
    {
        assert(edge_A->GetTargetNode()==edge_p->GetOriginNode());
        //Edges are merged such that the intersected element retains A to pNode edge
        //and pNode to B edge is inserted into it
        //NodeB is the common node
        HalfEdge<SPACE_DIM>* pNode_to_B = edge_p->GetTwinHalfEdge();
        assert(!pNode_to_B->GetElement());
        assert(pNode_to_B->GetTargetNode() == edge_A->GetTargetNode());

        pNode_to_B->SetElement(intersected_element);
        assert(!edge_A->GetTwinHalfEdge()->GetElement());
        assert(edge_A->GetTwinHalfEdge()->GetPreviousHalfEdge()==pNode_to_B);
        pNode_to_B->GetPreviousHalfEdge()->SetNextHalfEdge(edge_A->GetTwinHalfEdge(), true);
        pNode_to_B->SetNextHalfEdge(edge_A->GetNextHalfEdge(), true);
        pNode_to_B->SetPreviousHalfEdge(edge_A, true);

        edge_A->SetTargetNode(pNode);
    }
    else
    {
        assert(edge_A->GetOriginNode()==edge_p->GetOriginNode());
        //Edges are merged such that the intersected element retains pNode to B edge (i.e. A to B edge)
        //and pNode to A edge is inserted into it
        //nodeA is the common node
        assert(edge_p->GetTwinHalfEdge()->GetElement()==intersecting_element);

        assert(!edge_A->GetTwinHalfEdge()->GetElement());
        assert(edge_A->GetTwinHalfEdge()->GetNextHalfEdge()==edge_p);
        edge_p->GetNextHalfEdge()->SetPreviousHalfEdge(edge_A->GetTwinHalfEdge(), true);

        edge_p->SetElement(intersected_element);
        edge_p->SetPreviousHalfEdge(edge_A->GetPreviousHalfEdge(), true);
        edge_p->SetNextHalfEdge(edge_A, true);

        edge_A->GetTwinHalfEdge()->SetTargetNode(pNode);
    }
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::HandleHighOrderJunctions(FullEdge<SPACE_DIM>& full_edge)
{

}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::RemoveDeletedNodesAndElements(VertexElementMap& rElementMap)
{

}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::RemoveDeletedNodes()
{
    // Remove any nodes that have been marked for deletion and store all other nodes in a temporary structure
    std::vector<Node<SPACE_DIM>*> live_nodes;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            delete this->mNodes[i];
        }
        else
        {
            live_nodes.push_back(this->mNodes[i]);
        }
    }

    // Sanity check
    assert(this->mDeletedNodeIndices.size() == this->mNodes.size() - live_nodes.size());

    // Repopulate the nodes vector and reset the list of deleted node indices
    this->mNodes = live_nodes;
    this->mDeletedNodeIndices.clear();

    // Finally, reset the node indices to run from zero
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(i);
    }
}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::RemoveDeletedEdges()
{
    std::vector<FullEdge<SPACE_DIM>* > live_full_edges;
    for (unsigned int i=0; i<this->mFullEdges.size(); ++i)
    {
        HalfEdge<SPACE_DIM>* edge = this->mFullEdges[i]->at(0);
        if (edge->IsDeleted())
        {
            //Make sure the element does not contain pointer to the deleted edge
            if (edge->GetElement())
            {
                assert(edge->GetElement()->GetHalfEdge()!=edge);
            }
            HalfEdge<SPACE_DIM>* twin = edge->GetTwinHalfEdge();
            //Remove edges from the map
            this->mHalfToFullEdgeMap.erase(edge);
            this->mHalfToFullEdgeMap.erase(twin);
            delete edge;
            delete twin;
            delete this->mFullEdges[i];
        }
        else
        {
            live_full_edges.push_back(this->mFullEdges[i]);
        }
    }

    // Sanity check
    assert(this->mDeletedHalfEdges.size()/2 == this->mFullEdges.size() - live_full_edges.size());

    // Repopulate the nodes vector and reset the list of deleted node indices
    this->mFullEdges = live_full_edges;
    this->mDeletedHalfEdges.clear();

    // TODO: do we need to reset full edge indices?
}



template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM==2);     // LCOV_EXCL_LINE

    if (SPACE_DIM == 2)
    {
        // Make sure the map is big enough
        rElementMap.Resize(this->GetNumAllElements());

        /*
         * To begin the remeshing process, we do not need to call Clear() and remove all current data,
         * since cell birth, rearrangement and death result only in local remeshing of a vertex-based
         * mesh. Instead, we just remove any deleted elements and nodes.
         */
        RemoveDeletedNodesAndElements(rElementMap);
        bool recheck_mesh = true;
        while (recheck_mesh == true)
        {
            // We check for any short edges and perform swaps if necessary and possible.
            recheck_mesh = CheckForSwapsFromShortEdges();
        }

        // Check for element intersections
        recheck_mesh = true;
        while (recheck_mesh == true)
        {
            // Check mesh for intersections, and perform T3 swaps where required
            recheck_mesh = CheckForIntersections();
        }

        RemoveDeletedNodes();

        /*
         * This is handled in a separate method to allow child classes to implement additional ReMeshing functionality
         * (see #2664).
         */
        this->CheckForRosettes();
    }
    else // 3D
    {
        // LCOV_EXCL_START
        EXCEPTION("Remeshing has not been implemented in 3D (see #827 and #860)\n");
        // LCOV_EXCL_STOP
        ///\todo Implement ReMesh() in 3D (see #1422)
    }
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::ReMesh()
{
    VertexElementMap map(GetNumElements());
    ReMesh(map);
}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::PerformT1Swap(HalfEdge<SPACE_DIM>* pEdge)
{
    // First compute and store the location of the T1 swap, which is at the midpoint of nodes A and B
    double distance_between_nodes_CD = this->mCellRearrangementRatio*this->mCellRearrangementThreshold;
    HENode<SPACE_DIM>* pNodeA = pEdge->GetTargetNode();
    HENode<SPACE_DIM>* pNodeB = pEdge->GetOriginNode();
    c_vector<double, SPACE_DIM> nodeA_location = pNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> nodeB_location = pNodeB->rGetLocation();
    c_vector<double, SPACE_DIM> vector_AB = this->GetVectorFromAtoB(nodeA_location, nodeB_location);
    mLocationsOfT1Swaps.push_back(nodeA_location + 0.5*vector_AB);

    double distance_AB = norm_2(vector_AB);
    if (distance_AB < 1e-10) ///\todo remove magic number? (see #1884 and #2401)
    {
        EXCEPTION("Nodes are too close together, this shouldn't happen");
    }

    /*
     * Compute the locations of two new nodes C, D, placed on either side of the
     * edge E_old formed by nodes A and B, such that the edge E_new formed by the
     * new nodes is the perpendicular bisector of E_old, with |E_new| 'just larger'
     * (this->mCellRearrangementRatio) than mThresholdDistance.
     *
     * We implement the following changes to the mesh:
     *
     * The element whose index was in nodeA_elem_indices but not nodeB_elem_indices,
     * and the element whose index was in nodeB_elem_indices but not nodeA_elem_indices,
     * should now both contain nodes A and B.
     *
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which
     * node C lies inside, should now only contain node A.
     *
     * The element whose index was in nodeA_elem_indices and nodeB_elem_indices, and which
     * node D lies inside, should now only contain node B.
     *
     * Iterate over all elements involved and identify which element they are
     * in the diagram then update the nodes as necessary.
     *
     *   \(1)/
     *    \ / Node A
     * (2) |   (4)     elements in brackets
     *    / \ Node B
     *   /(3)\
     */

    // Move nodes A and B to C and D respectively
    c_vector<double, SPACE_DIM> vector_CD;
    vector_CD(0) = -vector_AB(1) * distance_between_nodes_CD / distance_AB;
    vector_CD(1) =  vector_AB(0) * distance_between_nodes_CD / distance_AB;

    c_vector<double, SPACE_DIM> nodeC_location = nodeA_location + 0.5*vector_AB - 0.5*vector_CD;
    c_vector<double, SPACE_DIM> nodeD_location = nodeC_location + vector_CD;

    pNodeA->rGetModifiableLocation() = nodeC_location;
    pNodeB->rGetModifiableLocation() = nodeD_location;

    nodeA_location = pNodeA->rGetLocation();
    nodeB_location = pNodeB->rGetLocation();

    HalfEdge<SPACE_DIM>* next_edge = pEdge->GetNextHalfEdge();
    HalfEdge<SPACE_DIM>* previous_edge = pEdge->GetPreviousHalfEdge();
    HalfEdge<SPACE_DIM>* twin = pEdge->GetTwinHalfEdge();
    HalfEdge<SPACE_DIM>* twin_next_edge = twin->GetNextHalfEdge();
    HalfEdge<SPACE_DIM>* twin_previous_edge = twin->GetPreviousHalfEdge();

    //Change adjacency relations. Traversal of element edges will not include collapsed edge
    previous_edge->SetNextHalfEdge(next_edge,true);
    twin_previous_edge->SetNextHalfEdge(twin_next_edge,true);
    HEElement<SPACE_DIM>* edge_A_element = pEdge->GetElement();
    HEElement<SPACE_DIM>* edge_B_element = twin->GetElement();
    if (edge_A_element)
    {
        if (edge_A_element->GetHalfEdge()==pEdge)
            edge_A_element->SetHalfEdge(previous_edge);
    }
    if (edge_B_element)
    {
        if (edge_B_element->GetHalfEdge()==twin)
            edge_B_element->SetHalfEdge(twin_previous_edge);
    }

    /*
     * If elements 1 and 3 are voids, we separate two cells by creating a void between
     *
     * |\   /|     |\      /|
     * | \ / |     | \    / |
     * |  |  |  => | /    \ |
     * | / \ |     |/      \|
     * |/   \|
     */
    if (!next_edge->GetTwinHalfEdge()->GetElement()&&!twin_next_edge->GetTwinHalfEdge()->GetElement())
    {
        assert(pEdge->GetTargetNode()->IsBoundaryNode());
        assert(twin->GetTargetNode()->IsBoundaryNode());
        //Modify adjacency relations
        next_edge->GetTwinHalfEdge()->SetNextHalfEdge(previous_edge->GetTwinHalfEdge(), true);
        twin_next_edge->GetTwinHalfEdge()->SetNextHalfEdge(twin_previous_edge->GetTwinHalfEdge(), true);

        //Modify outgoing edges
        pNodeA->SetOutgoingEdge(next_edge);
        pNodeB->SetOutgoingEdge(twin_next_edge);

        //Modify the target nodes
        //Edge in element 4 pointing to node A must now point to node B
        twin_previous_edge->SetTargetNode(pNodeB);

        //Edge in element 2 pointing to node B must now point to node A
        previous_edge->SetTargetNode(pNodeA);

        edge_A_element->UpdateGeometry();
        edge_B_element->UpdateGeometry();
    }
    else
    {
        //Create halfedges pointing to node C and D, respectively
        HalfEdge<SPACE_DIM>* pEdgeC = new HalfEdge<SPACE_DIM>();
        HalfEdge<SPACE_DIM>* pEdgeD = new HalfEdge<SPACE_DIM>();
        pEdgeC->SetTwinHalfEdge(pEdgeD,true);
        pEdgeC->SetTargetNode(pNodeA);
        pEdgeD->SetTargetNode(pNodeB);

        //Insert halfedge pointing to node D into element 1 and modify adjacency relations
        HalfEdge<SPACE_DIM>* out_A_edge = pNodeA->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* next_out_A_edge = out_A_edge;
        do
        {
            HEElement<SPACE_DIM>* out_A_element = next_out_A_edge->GetElement();
            //Find element 1
            if (out_A_element != edge_A_element && out_A_element != edge_B_element)
            {
                pEdgeD->SetElement(out_A_element);
                pEdgeD->SetPreviousHalfEdge(next_out_A_edge->GetPreviousHalfEdge(),true);
                pEdgeD->SetNextHalfEdge(next_out_A_edge, true);
                //Edge previously pointing to node A now must point to B
                pEdgeD->GetNextHalfEdge()->GetTwinHalfEdge()->SetTargetNode(pNodeB);
                pNodeA->SetOutgoingEdge(pEdgeD);
                break;
            }
            next_out_A_edge = next_out_A_edge->GetTwinHalfEdge()->GetNextHalfEdge();
        }while(next_out_A_edge != out_A_edge);

        //Insert halfedge pointing to node C into element 3 and modify adjacency relations
        HalfEdge<SPACE_DIM>* out_B_edge = pNodeB->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* next_out_B_edge = out_B_edge;
        do
        {
            HEElement<SPACE_DIM>* out_B_element = next_out_B_edge->GetElement();
            //Find element 3
            if (out_B_element != edge_A_element && out_B_element != edge_B_element)
            {
                pEdgeC->SetElement(out_B_element);
                pEdgeC->SetPreviousHalfEdge(next_out_B_edge->GetPreviousHalfEdge(),true);
                pEdgeC->SetNextHalfEdge(next_out_B_edge, true);
                //Edge previously pointing to node B now must point to A
                pEdgeC->GetNextHalfEdge()->GetTwinHalfEdge()->SetTargetNode(pNodeA);
                pNodeB->SetOutgoingEdge(pEdgeC);
                break;
            }
            next_out_B_edge = next_out_B_edge->GetTwinHalfEdge()->GetNextHalfEdge();
        }while(next_out_B_edge != out_B_edge);

        AddEdge(pEdgeC);
        this->MarkHalfEdgeAsDeleted(pEdge);
        // Sort out boundary nodes
        //Count the number of elements containing nodes A and B to determine if they should be set as boundary nodes
        unsigned int node_A_num_elements= edge_A_element ? pNodeA->GetContainingElementIndices().size() : 0;
        unsigned int node_B_num_elements= edge_B_element ? pNodeB->GetContainingElementIndices().size() : 0;
        if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
        {
            if (node_A_num_elements == 3)
            {
                pNodeA->SetAsBoundaryNode(false);
            }
            else
            {
                pNodeA->SetAsBoundaryNode(true);
            }
            if (node_B_num_elements == 3)
            {
                pNodeB->SetAsBoundaryNode(false);
            }
            else
            {
                pNodeB->SetAsBoundaryNode(true);
            }
        }

        if (edge_A_element)
        {
            edge_A_element->UpdateGeometry();
        }
        if (edge_B_element)
        {
            edge_B_element->UpdateGeometry();
        }
        if (pEdgeC->GetElement())
        {
            pEdgeC->GetElement()->UpdateGeometry();
        }
        if (pEdgeD->GetElement())
        {
            pEdgeD->GetElement()->UpdateGeometry();
        }
    }
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::PerformIntersectionSwap(HENode<SPACE_DIM>* pNode, unsigned elementIndex)
{
    assert(SPACE_DIM == 2);                    // LCOV_EXCL_LINE

    HEElement<SPACE_DIM>* p_element = this->GetElement(elementIndex);

    HalfEdge<SPACE_DIM>* edge = p_element->GetHalfEdge();
    HalfEdge<SPACE_DIM>* next_edge = edge;
    std::set<unsigned> elements_containing_intersecting_node;

    do
    {
        HalfEdge<SPACE_DIM>* twin = next_edge->GetTwinHalfEdge();
        HEElement<SPACE_DIM>* neigh_element = twin->GetElement();
        if (neigh_element)
        {
            HalfEdge<SPACE_DIM>* neigh_next_edge = twin;
            do
            {
                if (neigh_next_edge->GetTargetNode()->GetIndex() == pNode->GetIndex())
                {
                    elements_containing_intersecting_node.insert(neigh_element->GetIndex());
                }
                neigh_next_edge = neigh_next_edge->GetNextHalfEdge();
            }while(neigh_next_edge != twin);
        }
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != edge);

    /*
     * If there are not two elements containing the intersecting node then the node is coming from the other side of the element
     * and there is no way to fix it unless you want to make two new elements.
     */
    assert(elements_containing_intersecting_node.size() == 2);

    std::set<unsigned> all_elements_containing_intersecting_node = pNode->GetContainingElementIndices();

    std::set<unsigned> intersecting_element;

    std::set_difference(all_elements_containing_intersecting_node.begin(), all_elements_containing_intersecting_node.end(),
                        elements_containing_intersecting_node.begin(), elements_containing_intersecting_node.end(),
                        std::inserter(intersecting_element, intersecting_element.begin()));

    /*
     * Identify nodes and elements to perform switch on
     * Intersecting node is node A
     * Other node is node B
     *
     * Element 1 only contains node A
     * Element 2 has nodes B and A (in that order)
     * Element 3 only contains node B
     * Element 4 has nodes A and B (in that order)
     */
    unsigned node_A_index = pNode->GetIndex();
    unsigned node_B_index;
    unsigned element_1_index = *(intersecting_element.begin());
    unsigned element_2_index;
    unsigned element_3_index = elementIndex;
    unsigned element_4_index;

    std::set<unsigned>::iterator iter = elements_containing_intersecting_node.begin();
    unsigned element_a_index = *(iter);
    iter++;
    unsigned element_b_index = *(iter);

    HEElement<SPACE_DIM>* p_element_a = this->GetElement(element_a_index);
    HEElement<SPACE_DIM>* p_element_b = this->GetElement(element_b_index);

    std::set<unsigned> element_a_nodes;
    HalfEdge<SPACE_DIM>* start_edge = p_element_a->GetHalfEdge();
    next_edge = start_edge;
    do
    {
        element_a_nodes.insert(next_edge->GetTargetNode()->GetIndex());
        next_edge = next_edge ->GetNextHalfEdge();
    }while(next_edge != start_edge);

    std::set<unsigned> element_b_nodes;
    start_edge = p_element_b->GetHalfEdge();
    next_edge = start_edge;
    do
    {
        element_b_nodes.insert(next_edge->GetTargetNode()->GetIndex());
        next_edge = next_edge ->GetNextHalfEdge();
    }while(next_edge != start_edge);

    std::set<unsigned> switching_nodes;
    std::set_intersection(element_a_nodes.begin(), element_a_nodes.end(),
                          element_b_nodes.begin(), element_b_nodes.end(),
                          std::inserter(switching_nodes, switching_nodes.begin()));

    assert(switching_nodes.size() == 2);

    // Check intersecting node is this set
    assert(switching_nodes.find(node_A_index) != switching_nodes.end());
    switching_nodes.erase(node_A_index);

    assert(switching_nodes.size() == 1);

    node_B_index = *(switching_nodes.begin());

    // Now identify elements 2 and 4

    unsigned node_A_local_index_in_a = p_element_a->GetNodeLocalIndex(node_A_index);
    unsigned node_B_local_index_in_a = p_element_a->GetNodeLocalIndex(node_B_index);

    if ((node_B_local_index_in_a+1)%p_element_a->GetNumNodes() == node_A_local_index_in_a)
    {
        assert((p_element_b->GetNodeLocalIndex(node_A_index)+1)%p_element_b->GetNumNodes()
               == p_element_b->GetNodeLocalIndex(node_B_index));

        // Element 2 is element a, element 4 is element b
        element_2_index = element_a_index;
        element_4_index = element_b_index;
    }
    else
    {
        assert((p_element_b->GetNodeLocalIndex(node_B_index)+1)%p_element_b->GetNumNodes()
               == p_element_b->GetNodeLocalIndex(node_A_index));

        // Element 2 is element b, element 4 is element a
        element_2_index = element_b_index;
        element_4_index = element_a_index;
    }

    HEElement<SPACE_DIM>* element_1 = this->mElements[element_1_index];
    HEElement<SPACE_DIM>* element_2 = this->mElements[element_2_index];
    HEElement<SPACE_DIM>* element_3 = this->mElements[element_3_index];
    HEElement<SPACE_DIM>* element_4 = this->mElements[element_4_index];
    unsigned int element_1_nodes = element_1->GetNumNodes();
    unsigned int element_2_nodes = element_2->GetNumNodes();
    unsigned int element_3_nodes = element_3->GetNumNodes();
    unsigned int element_4_nodes = element_4->GetNumNodes();

    unsigned intersected_edge = this->GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), elementIndex);

    unsigned node_A_local_index_in_1 = element_1->GetNodeLocalIndex(node_A_index);
    unsigned node_B_local_index_in_2 = element_2->GetNodeLocalIndex(node_B_index);
    unsigned node_B_local_index_in_3 = element_3->GetNodeLocalIndex(node_B_index);
    unsigned node_A_local_index_in_4 = element_4->GetNodeLocalIndex(node_A_index);

    HENode<SPACE_DIM>* node_B = static_cast<HENode<SPACE_DIM>* >(this->mNodes[node_B_index]);
    HENode<SPACE_DIM>* node_A = static_cast<HENode<SPACE_DIM>* >(this->mNodes[node_A_index]);
    HalfEdge<SPACE_DIM>* edge_A_to_B = element_2->GetHalfEdge(node_A_local_index_in_4);
    if (edge_A_to_B->GetTargetNode()==node_A)
        edge_A_to_B = edge_A_to_B->GetTwinHalfEdge();
    assert(edge_A_to_B->GetTargetNode()==node_B);
    HalfEdge<SPACE_DIM>* edge_B_to_A = edge_A_to_B->GetTwinHalfEdge();

    unsigned node_before_A_in_1 = (node_A_local_index_in_1 + element_1_nodes - 1)%element_1_nodes;
    unsigned node_before_B_in_2 = (node_B_local_index_in_2+element_2_nodes-1)%element_2_nodes;
    unsigned node_before_B_in_3 = (node_B_local_index_in_3 + element_3_nodes- 1)%element_3_nodes;
    unsigned node_before_A_in_4 = (node_A_local_index_in_4+element_4_nodes-1)%element_4_nodes;

    HalfEdge<SPACE_DIM>* edge_in_1_to_A = element_1->GetHalfEdge(node_before_A_in_1);
    HalfEdge<SPACE_DIM>* twin_edge_to_A = edge_in_1_to_A->GetTwinHalfEdge();
    assert(twin_edge_to_A==edge_B_to_A->GetNextHalfEdge());

    HalfEdge<SPACE_DIM>* edge_in_3_to_B = element_3->GetHalfEdge(node_before_B_in_3);
    HalfEdge<SPACE_DIM>* twin_edge_to_B = edge_in_3_to_B->GetTwinHalfEdge();
    assert(twin_edge_to_B==edge_A_to_B->GetNextHalfEdge());

    //Element 2
    HalfEdge<SPACE_DIM>* edge_in_2_to_B = element_2->GetHalfEdge(node_before_B_in_2);
    edge_in_2_to_B->SetNextHalfEdge(twin_edge_to_A, true);

    //Element 4
    HalfEdge<SPACE_DIM>* edge_in_4_to_A = element_4->GetHalfEdge(node_before_A_in_4);
    edge_in_4_to_A->SetNextHalfEdge(twin_edge_to_B, true);
    if (intersected_edge==node_B_local_index_in_3)
    {
        /*
         * Add node B to element 1 after node A
         * Add node A to element 3 after node B
         *
         * Remove node B from element 2
         * Remove node A from element 4
         */

        //Element 1

        node_A->SetOutgoingEdge(edge_A_to_B);
        edge_A_to_B->SetNextHalfEdge(edge_in_1_to_A->GetNextHalfEdge(),true);
        edge_A_to_B->SetPreviousHalfEdge(edge_in_1_to_A, true);
        edge_A_to_B->SetElement(element_1);

        //Element 3
        node_B->SetOutgoingEdge(edge_B_to_A);
        edge_B_to_A->SetNextHalfEdge(edge_in_3_to_B->GetNextHalfEdge(),true);
        edge_B_to_A->SetPreviousHalfEdge(edge_in_3_to_B, true);
        edge_B_to_A->SetElement(element_3);

        //Element 2
        edge_in_2_to_B->SetTargetNode(node_A);

        //Element 4
        edge_in_4_to_A->SetTargetNode(node_B);
    }
    else
    {
        assert((intersected_edge+1)%p_element->GetNumNodes()==node_B_local_index_in_3);

        // Add node B to element 1 before node A and add node A to element 3 before node B

        //Element 1
        edge_in_1_to_A->SetTargetNode(node_B);
        node_B->SetOutgoingEdge(edge_B_to_A);
        edge_B_to_A->SetNextHalfEdge(edge_in_1_to_A->GetNextHalfEdge(),true);
        edge_B_to_A->SetPreviousHalfEdge(edge_in_1_to_A, true);
        edge_B_to_A->SetElement(element_1);

        //Element 3
        edge_in_3_to_B->SetTargetNode(node_A);
        node_A->SetOutgoingEdge(edge_A_to_B);
        edge_A_to_B->SetNextHalfEdge(edge_in_3_to_B->GetNextHalfEdge(),true);
        edge_A_to_B->SetPreviousHalfEdge(edge_in_3_to_B, true);
        edge_A_to_B->SetElement(element_3);
    }


    element_1->UpdateGeometry();
    element_2->UpdateGeometry();
    element_3->UpdateGeometry();
    element_4->UpdateGeometry();
}


template<unsigned SPACE_DIM>
bool HEMutableVertexMesh<SPACE_DIM>::CheckForT2Swaps(VertexElementMap& rElementMap)
{
    // Loop over elements to check for T2 swaps
    for (typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem_iter = this->GetElementIteratorBegin();
            elem_iter != this->GetElementIteratorEnd();
            ++elem_iter)
    {
        // If this element is triangular...
        if (elem_iter->GetNumNodes() == 3)
        {
            // ...and smaller than the threshold area...
            if (elem_iter->GetVolume() < this->GetT2Threshold())
            {
                // ...then perform a T2 swap and break out of the loop
                PerformT2Swap(*elem_iter);
                ///\todo: cover this line in a test
                rElementMap.SetDeleted(elem_iter->GetIndex());
                return true;
            }
        }

    }
    return false;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::PerformT2Swap(HEElement<SPACE_DIM>& rElement)
{
    assert(rElement.GetNumNodes() == 3);
    // Note that we define this vector before setting it, as otherwise the profiling build will break (see #2367)
    c_vector<double, SPACE_DIM> new_node_location;
    new_node_location = this->GetCentroidOfElement(rElement.GetIndex());
    mLastT2SwapLocation = new_node_location;

    bool is_node_on_boundary = false;
    //Loop over the nodes
    HalfEdge<SPACE_DIM>* edge = rElement.GetHalfEdge();
    HalfEdge<SPACE_DIM>* next_edge = edge;
    do
    {
        //If any of the nodes are on the boundary, the new node will also be on the boundary
        if (next_edge->GetTargetNode()->IsBoundaryNode())
        {
            is_node_on_boundary = true;
            break;
        }
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != edge);

    HENode<SPACE_DIM>* p_new_node = RemoveTriangle(edge);
    p_new_node->SetAsBoundaryNode(is_node_on_boundary);
    p_new_node->rGetModifiableLocation() = new_node_location;

    this->mDeletedElementIndices.push_back(rElement.GetIndex());
    rElement.MarkAsDeleted();
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::PerformT3Swap(HENode<SPACE_DIM>* pNode, unsigned elementIndex)
{
    assert(SPACE_DIM == 2);                 // LCOV_EXCL_LINE - code will be removed at compile time
    assert(pNode->IsBoundaryNode());

    // Store the index of the elements containing the intersecting node
    std::set<unsigned> elements_containing_intersecting_node = pNode->GetContainingElementIndices();

    std::set<unsigned int> affected_element_indices;
    affected_element_indices.insert(elementIndex);
    for (unsigned int index: elements_containing_intersecting_node)
        affected_element_indices.insert(index);

    // Get the local index of the node in the intersected element after which the new node is to be added
    unsigned node_A_local_index = this->GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), elementIndex);

    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
    c_vector<double, SPACE_DIM> node_location;
    node_location = pNode->rGetModifiableLocation();

    // Get element
    HEElement<SPACE_DIM>* p_element = this->GetElement(elementIndex);
    unsigned num_nodes = p_element->GetNumNodes();

    HENode<SPACE_DIM>* node_A = p_element->GetNode(node_A_local_index);
    HalfEdge<SPACE_DIM>* edge_from_A = p_element->GetHalfEdge(node_A_local_index);
    HENode<SPACE_DIM>* node_B = edge_from_A->GetTargetNode();
    // Get the nodes at either end of the edge to be divided
    unsigned vertexA_index = node_A->GetIndex();
    unsigned vertexB_index = node_B->GetIndex();

    // Check these nodes are also boundary nodes if this fails then the elements have become concave and you need a smaller timestep
    if (!this->mNodes[vertexA_index]->IsBoundaryNode() || !this->mNodes[vertexB_index]->IsBoundaryNode())
    {
        EXCEPTION("A boundary node has intersected a non-boundary edge; this is because the boundary element has become concave. You need to rerun the simulation with a smaller time step to prevent this.");
    }

    // Get the nodes at either end of the edge to be divided and calculate intersection
    c_vector<double, SPACE_DIM> vertexA = node_A->rGetLocation();
    c_vector<double, SPACE_DIM> vertexB = node_B->rGetLocation();
    c_vector<double, SPACE_DIM> vector_a_to_point = this->GetVectorFromAtoB(vertexA, node_location);

    c_vector<double, SPACE_DIM> vector_a_to_b = this->GetVectorFromAtoB(vertexA, vertexB);

    c_vector<double, SPACE_DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);
    c_vector<double, SPACE_DIM> intersection = vertexA + edge_ab_unit_vector*inner_prod(vector_a_to_point, edge_ab_unit_vector);

    // Store the location of the T3 swap, the location of the intersection with the edge
    ///\todo the intersection location is sometimes overwritten when WidenEdgeOrCorrectIntersectionLocationIfNecessary
    // is called (see #2401) - we should correct this in these cases!

    mLocationsOfT3Swaps.push_back(intersection);

    if (elements_containing_intersecting_node.size() == 1)
    {
        // Get the index of the element containing the intersecting node
        unsigned intersecting_element_index = *elements_containing_intersecting_node.begin();

        // Get element
        HEElement<SPACE_DIM>* p_intersecting_element = this->GetElement(intersecting_element_index);

        unsigned local_index = p_intersecting_element->GetNodeLocalIndex(pNode->GetIndex());
        HalfEdge<SPACE_DIM>* edge_from_pNode = p_intersecting_element->GetHalfEdge(local_index);
        unsigned next_node = edge_from_pNode->GetTargetNode()->GetIndex();
        unsigned previous_node = edge_from_pNode->GetPreviousHalfEdge()->GetOriginNode()->GetIndex();

        // Check to see if the nodes adjacent to the intersecting node are contained in the intersected element between vertices A and B
        if (next_node == vertexA_index || previous_node == vertexA_index || next_node == vertexB_index || previous_node == vertexB_index)
        {
            unsigned common_vertex_index;

            if (next_node == vertexA_index || previous_node == vertexA_index)
            {
                common_vertex_index = vertexA_index;
            }
            else
            {
                common_vertex_index = vertexB_index;
            }

            HENode<SPACE_DIM>* common_vertex = static_cast<HENode<SPACE_DIM>* >(this->mNodes[common_vertex_index]);
            std::set<unsigned> elements_containing_common_vertex = common_vertex->GetContainingElementIndices();
            assert(elements_containing_common_vertex.size()>1);

            //Register affected elements before deleting the common vertex
            for (unsigned int index : elements_containing_common_vertex)
                affected_element_indices.insert(index);

            std::set<unsigned>::const_iterator it = elements_containing_common_vertex.begin();
            HEElement<SPACE_DIM>* p_element_common_1 = this->GetElement(*it);
            it++;
            HEElement<SPACE_DIM>* p_element_common_2 = this->GetElement(*it);

            // Find the number and indices of common vertices between element_1 and element_2
            unsigned num_common_vertices = 0;
            std::vector<unsigned> common_vertex_indices;
            for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
            {
                for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                {
                    if (p_element_common_1->GetNodeGlobalIndex(i)==p_element_common_2->GetNodeGlobalIndex(j))
                    {
                        num_common_vertices++;
                        common_vertex_indices.push_back(p_element_common_1->GetNodeGlobalIndex(i));
                    }
                }
            }

            if (num_common_vertices == 1 || elements_containing_common_vertex.size() > 2)
            {
                /*
                 * This is the situation here.
                 *
                 *  From          To
                 *   _             _
                 *    |  <---       |
                 *    |  /\         |\
                 *    | /  \        | \
                 *   _|/____\      _|__\
                 *
                 * The edge goes from vertexA--vertexB to vertexA--pNode--vertexB
                 */

                // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(edge_from_A, intersection);

                // Move original node
                pNode->rGetModifiableLocation() = intersection;

                HalfEdge<SPACE_DIM>* out_edge = common_vertex->GetOutgoingEdge();
                HalfEdge<SPACE_DIM>* next_out_edge = out_edge;
                HalfEdge<SPACE_DIM>* edge_to_pNode = nullptr;
                do
                {
                    if (next_out_edge->GetTargetNode()==pNode)
                    {
                        edge_to_pNode = next_out_edge;
                        break;
                    }
                    next_out_edge = next_out_edge->GetTwinHalfEdge()->GetNextHalfEdge();
                }while(next_out_edge != out_edge);
                assert(edge_to_pNode);

                MergeEdgesInT3Swap(edge_from_A, edge_to_pNode, pNode);
                //If common vertex is A, then merge as above. I.e. edge_to_pNode->pNode and pNode->edge_from_A->vertex_B
                //Else, vertex_A->edge_from_A->pNode->edge_from_pNode->Vertex_B

                // Check the nodes are updated correctly
                assert(pNode->GetContainingElementIndices().size() == 2);
            }
            else if (num_common_vertices == 2)
            {
                // The two elements must have an edge in common.  Find whether the common edge is the same as the
                // edge that is merged onto.

                if ((common_vertex_indices[0]==vertexA_index && common_vertex_indices[1]==vertexB_index) ||
                        (common_vertex_indices[1]==vertexA_index && common_vertex_indices[0]==vertexB_index))
                {
                    /*
                     * Due to a previous T3 swap the situation looks like this.
                     *
                     *              pNode
                     *     \         |\    /
                     *      \        | \  /
                     *       \_______|__\/
                     *       /A      |     B
                     *      /         \
                     *
                     * A T3 Swap would merge pNode onto an edge of its own element.
                     * We prevent this by just removing pNode. By doing this we also avoid the
                     * intersecting element to be concave.
                     */

                    // Delete pNode in the intersecting element
                    std::set<HalfEdge<SPACE_DIM>* > deleted_edge = p_intersecting_element->DeleteNode(pNode);
                    pNode->MarkAsDeleted();
                    this->mDeletedNodeIndices.push_back(pNode->GetIndex());
                    for (auto edge:deleted_edge)
                    this->MarkHalfEdgeAsDeleted(edge);
                }
                else
                {
                    /*
                     * This is the situation here.
                     *
                     * C is common_vertex D is the other one.
                     *
                     *  From          To
                     *   _ D          _
                     *    | <---       |
                     *    | /\         |\
                     *   C|/  \        | \
                     *   _|____\      _|__\
                     *
                     *  The edge goes from vertexC--vertexB to vertexC--pNode--vertexD
                     *  then vertex B is removed as it is no longer needed.
                     */

                    //Reduce the situation to the case with one common vertex by deleting the common edge

                    HENode<SPACE_DIM>* deleted_node = edge_from_A->GetPreviousHalfEdge()->GetOriginNode();
                    if (edge_from_A->GetNextHalfEdge()->GetTwinHalfEdge()->GetElement() == p_intersecting_element)
                    {
                        deleted_node = edge_from_A->GetNextHalfEdge()->GetTargetNode();
                    }

                    std::set<unsigned int> deleted_node_element_indices = deleted_node->GetContainingElementIndices();
                    for (unsigned int index : deleted_node_element_indices)
                        affected_element_indices.insert(index);

                    common_vertex->rGetModifiableLocation() = deleted_node->rGetLocation();
                    std::set<HalfEdge<SPACE_DIM>* > deleted_edges = p_element->DeleteNode(deleted_node);
                    for (unsigned int index: affected_element_indices)
                    {
                        this->GetElement(index)->UpdateGeometry();
                    }

                    this->mDeletedNodeIndices.push_back(deleted_node->GetIndex());
                    deleted_node->MarkAsDeleted();
                    for (auto edges : deleted_edges)
                        this->MarkHalfEdgeAsDeleted(edges);

                    //Make sure we are in a situation with one common vertex
                    assert(pNode->GetNumContainingElements()==1);
                    num_common_vertices = 0;
                    for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                    {
                        for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                        {
                            if (p_element_common_1->GetNodeGlobalIndex(i)==p_element_common_2->GetNodeGlobalIndex(j))
                            {
                                num_common_vertices++;
                                common_vertex_indices.push_back(p_element_common_1->GetNodeGlobalIndex(i));
                            }
                        }
                    }
                    elements_containing_common_vertex=common_vertex->GetContainingElementIndices();
                    assert(num_common_vertices == 1 || elements_containing_common_vertex.size() > 2);
                    assert(num_common_vertices != 2);
                    local_index = p_intersecting_element->GetNodeLocalIndex(pNode->GetIndex());
                     edge_from_pNode = p_intersecting_element->GetHalfEdge(local_index);
                    next_node = edge_from_pNode->GetTargetNode()->GetIndex();
                    previous_node = edge_from_pNode->GetPreviousHalfEdge()->GetOriginNode()->GetIndex();
                    const bool prev_case
                    = next_node == vertexA_index || previous_node == vertexA_index || next_node == vertexB_index || previous_node == vertexB_index;
                    assert(prev_case);

                    this->PerformT3Swap(pNode, elementIndex);
                }
            }
            else if (num_common_vertices == 4)
            {
                /*
                 * The two elements share edges CA and BD due to previous swaps but not the edge AB
                 *
                 *  From          To
                 *  D___         D___
                 *    |            |
                 *   B|\           |
                 *    | \          |
                 *    | /          |
                 *   A|/           |
                 *  C_|__        C_|__
                 *
                 *  We just remove the intersecting node as well as vertices A and B.
                 */

                //pNode to A edge is deleted, leaving B to A edge intact in the intersecting element
                std::set<HalfEdge<SPACE_DIM>* > pNode_deleted_edges = p_intersecting_element->DeleteNode(pNode);
                assert(pNode_deleted_edges.size()==1);
                assert((*pNode_deleted_edges.begin())->GetTargetNode()==node_A);

                //We merge the duplicate edges, such that the elements share AB edge, and delete the nodes A and B
                HalfEdge<SPACE_DIM>* inter_BA_edge = (*pNode_deleted_edges.begin())->GetPreviousHalfEdge();
                HalfEdge<SPACE_DIM>* twin_edge_from_A = edge_from_A->GetTwinHalfEdge();

                twin_edge_from_A->SetPreviousHalfEdge(inter_BA_edge->GetPreviousHalfEdge(), true);
                twin_edge_from_A->SetNextHalfEdge(inter_BA_edge->GetNextHalfEdge(), true);
                twin_edge_from_A->SetElement(p_intersecting_element);
                twin_edge_from_A->GetOriginNode()->SetOutgoingEdge(twin_edge_from_A);
                node_A->SetOutgoingEdge(edge_from_A);
                if (p_intersecting_element->GetHalfEdge()==inter_BA_edge)
                    p_intersecting_element->SetHalfEdge(twin_edge_from_A);

                //Delete nodes A and B
                std::set<HalfEdge<SPACE_DIM>* > node_A_deleted_edges = p_intersecting_element->DeleteNode(node_A);

                std::set<HalfEdge<SPACE_DIM>* > node_B_deleted_edges = p_intersecting_element->DeleteNode(node_B);

                // Mark all three nodes as deleted
                pNode->MarkAsDeleted();
                this->mDeletedNodeIndices.push_back(pNode->GetIndex());
                this->mNodes[vertexA_index]->MarkAsDeleted();
                this->mDeletedNodeIndices.push_back(vertexA_index);
                this->mNodes[vertexB_index]->MarkAsDeleted();
                this->mDeletedNodeIndices.push_back(vertexB_index);

                //Mark deleted edges
                this->MarkHalfEdgeAsDeleted(inter_BA_edge);
                this->MarkHalfEdgeAsDeleted(*pNode_deleted_edges.begin());
                for(auto edge: node_A_deleted_edges)
                    this->MarkHalfEdgeAsDeleted(edge);
                for(auto edge: node_B_deleted_edges)
                    this->MarkHalfEdgeAsDeleted(edge);
            }
            else
            {
                // This can't happen as nodes can't be on the internal edge of 2 elements.
                NEVER_REACHED;
            }
        }
        else
        {
            /*
             *  From          To
             *   ____        _______
             *                 / \
             *    /\   ^      /   \
             *   /  \  |
             *
             *  The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
             */

            // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
            intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(edge_from_A, intersection);
            edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

            // Move original node
            pNode->rGetModifiableLocation() = intersection + 0.5*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;

            // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
            c_vector<double, SPACE_DIM> new_node_location;
            new_node_location = intersection - 0.5*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;

            // Add new node which will always be a boundary node
            unsigned new_node_global_index = AddNode(new HENode<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]));
            HENode<SPACE_DIM>* new_node = static_cast<HENode<SPACE_DIM>* >(this->mNodes[new_node_global_index]);

            HalfEdge<SPACE_DIM>* edge_to_new_node =  p_element->AddNode(edge_from_A, new_node);
            this->AddEdge(edge_to_new_node);
            //Get the edge to pNode before adding the node to the intersected element
            HalfEdge<SPACE_DIM>* edge_to_pNode = p_intersecting_element->GetHalfEdge(pNode);

            HalfEdge<SPACE_DIM>* edge_from_new_to_pNode = p_element->AddNode(edge_to_new_node->GetNextHalfEdge(), pNode);
            this->AddEdge(edge_from_new_to_pNode);

            //Modify adjacency relations
            //edge_from_A is now pNode to B edge
            HalfEdge<SPACE_DIM>* next_edge_to_pNode = edge_to_pNode->GetNextHalfEdge();
            edge_from_A->GetTwinHalfEdge()->SetNextHalfEdge(edge_to_pNode->GetTwinHalfEdge(), true);
            //Intersecting element adjacencies
            edge_from_new_to_pNode->GetTwinHalfEdge()->SetElement(p_intersecting_element);
            edge_from_new_to_pNode->GetTwinHalfEdge()->SetPreviousHalfEdge(edge_to_pNode,true);
            edge_from_new_to_pNode->GetTwinHalfEdge()->SetNextHalfEdge(next_edge_to_pNode, true);
            edge_to_new_node->GetTwinHalfEdge()->SetPreviousHalfEdge(next_edge_to_pNode->GetTwinHalfEdge(), true);
            next_edge_to_pNode->GetTwinHalfEdge()->SetTargetNode(new_node);
            // The nodes must have been updated correctly
            assert(pNode->GetContainingElementIndices().size() == 2);
            assert(new_node->GetContainingElementIndices().size() == 2);
        }
    }
    else if (elements_containing_intersecting_node.size() == 2)
    {
        // Find the nodes contained in elements containing the intersecting node
        std::set<unsigned>::const_iterator it = elements_containing_intersecting_node.begin();

        HEElement<SPACE_DIM>* p_element_1 = this->GetElement(*it);
        unsigned num_nodes_elem_1 = p_element_1->GetNumNodes();
        it++;

        HEElement<SPACE_DIM>* p_element_2 = this->GetElement(*it);
        unsigned num_nodes_elem_2 = p_element_2->GetNumNodes();

        unsigned node_global_index = pNode->GetIndex();

        unsigned local_index_1 = p_element_1->GetNodeLocalIndex(node_global_index);
        unsigned next_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + 1)%num_nodes_elem_1);
        unsigned previous_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + num_nodes_elem_1 - 1)%num_nodes_elem_1);

        unsigned local_index_2 = p_element_2->GetNodeLocalIndex(node_global_index);
        unsigned next_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + 1)%num_nodes_elem_2);
        unsigned previous_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + num_nodes_elem_2 - 1)%num_nodes_elem_2);

        // Check to see if the nodes adjacent to the intersecting node are contained in the intersected element between vertices A and B
        if ((next_node_1 == vertexA_index || previous_node_1 == vertexA_index || next_node_2 == vertexA_index || previous_node_2 == vertexA_index) &&
                (next_node_1 == vertexB_index || previous_node_1 == vertexB_index || next_node_2 == vertexB_index || previous_node_2 == vertexB_index))
        {
            /*
             * Here we have
             *        __
             *      /|             /
             *  __ / |     --> ___/
             *     \ |            \
             *      \|__           \
             *
             * Where the node on the left has overlapped the edge A B
             *
             * Move p_node to the intersection on A B and merge AB and p_node
             */

            // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
            intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(edge_from_A, intersection);
            edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

            // Check they are all boundary nodes
            assert(pNode->IsBoundaryNode());
            assert(this->mNodes[vertexA_index]->IsBoundaryNode());
            assert(this->mNodes[vertexB_index]->IsBoundaryNode());

            //Register affected elements
            std::set<unsigned int> node_A_elements = node_A->GetContainingElementIndices();
            std::set<unsigned int> node_B_elements = node_B->GetContainingElementIndices();
            for (unsigned int index:node_A_elements)
                affected_element_indices.insert(index);
            for (unsigned int index:node_B_elements)
                affected_element_indices.insert(index);

            // Move p_node to the intersection with the edge AB
            pNode->rGetModifiableLocation() = intersection;
            pNode->SetAsBoundaryNode(false);

            //Remove AB edge, modify adjacency relations, and delete nodes A,B
            HalfEdge<SPACE_DIM>* prev_edge_from_A = edge_from_A->GetPreviousHalfEdge();
            HalfEdge<SPACE_DIM>* next_edge_from_A = edge_from_A->GetNextHalfEdge();
            HalfEdge<SPACE_DIM>* twin_edge_from_A = edge_from_A->GetTwinHalfEdge();
            assert(twin_edge_from_A->GetNextHalfEdge()->GetTargetNode()==pNode);
            assert(twin_edge_from_A->GetPreviousHalfEdge()->GetOriginNode()==pNode);

            //Remove AB edge and modify adjacency relations
            prev_edge_from_A->SetNextHalfEdge(twin_edge_from_A->GetNextHalfEdge(),true);
            next_edge_from_A->SetPreviousHalfEdge(twin_edge_from_A->GetPreviousHalfEdge(),true);
            node_A->SetOutgoingEdge(twin_edge_from_A->GetNextHalfEdge());
            node_B->SetOutgoingEdge(next_edge_from_A);
            if (p_element->GetHalfEdge()==edge_from_A)
                p_element->SetHalfEdge(prev_edge_from_A);

            //AP, PB edges belong to the intersected element
            prev_edge_from_A->GetNextHalfEdge()->SetElement(p_element);
            next_edge_from_A->GetPreviousHalfEdge()->SetElement(p_element);

            //Delete nodes A,B
            std::set<HalfEdge<SPACE_DIM>* > deleted_edges_A = p_element->DeleteNode(node_A);
            std::set<HalfEdge<SPACE_DIM>* > deleted_edges_B = p_element->DeleteNode(node_B);
            for (auto edges: deleted_edges_A)
                this->MarkHalfEdgeAsDeleted(edges);
            for (auto edges: deleted_edges_B)
                this->MarkHalfEdgeAsDeleted(edges);

            // Remove vertex A from the mesh
            this->mNodes[vertexA_index]->MarkAsDeleted();
            this->mDeletedNodeIndices.push_back(vertexA_index);

            // Remove vertex B from the mesh
            this->mNodes[vertexB_index]->MarkAsDeleted();
            this->mDeletedNodeIndices.push_back(vertexB_index);
        }
        else
        {
            if (next_node_1 == vertexA_index || previous_node_1 == vertexA_index || next_node_2 == vertexA_index || previous_node_2 == vertexA_index)
            {
                // Get elements containing vertexA_index (the common vertex)
                std::set<unsigned> elements_containing_vertex_A = node_A->GetContainingElementIndices();
                assert(elements_containing_vertex_A.size() > 1);

                std::set<unsigned>::const_iterator iter = elements_containing_vertex_A.begin();
                HEElement<SPACE_DIM>* p_element_common_1 = this->GetElement(*iter);
                iter++;
                HEElement<SPACE_DIM>* p_element_common_2 = this->GetElement(*iter);

                // Calculate the number of common vertices between element_1 and element_2
                unsigned num_common_vertices = 0;
                for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                {
                    for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                    {
                        if (p_element_common_1->GetNodeGlobalIndex(i) == p_element_common_2->GetNodeGlobalIndex(j))
                        {
                            num_common_vertices++;
                        }
                    }
                }

                if (num_common_vertices == 1 || node_A->GetNumContainingElements() > 2)
                {
                    /*
                     *  From          To
                     *   _ B              _ B
                     *    |  <---          |
                     *    |   /|\          |\
                     *    |  / | \         | \
                     *    | /  |  \        |\ \
                     *   _|/___|___\      _|_\_\
                     *     A                A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--pNode--new_node--vertexB
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(edge_from_A, intersection);

                    edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                    // Move original node and change to non-boundary node
                    pNode->rGetModifiableLocation() = intersection - 0.5*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;
                    pNode->SetAsBoundaryNode(false);

                    //Merge A--pNode edge first
                    HalfEdge<SPACE_DIM>* out_edge = node_A->GetOutgoingEdge();
                    HalfEdge<SPACE_DIM>* next_out_edge = out_edge;
                    HalfEdge<SPACE_DIM>* edge_to_pNode = nullptr;
                    do
                    {
                        if (next_out_edge->GetTargetNode()==pNode)
                        {
                            edge_to_pNode = next_out_edge;
                            break;
                        }
                        next_out_edge = next_out_edge->GetTwinHalfEdge()->GetNextHalfEdge();
                    }while(next_out_edge != out_edge);
                    assert(edge_to_pNode);
                    assert(!edge_to_pNode->GetNextHalfEdge()->GetElement());
                    MergeEdgesInT3Swap(edge_from_A, edge_to_pNode, pNode);

                    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                    c_vector<double, SPACE_DIM> new_node_location;
                    new_node_location = intersection + 0.5*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;

                    // Add new node, which will always be a boundary node
                    HENode<SPACE_DIM>* new_node = new HENode<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]);
                    AddNode(new_node);

                    //We need to store the edges below, because inserting the node modifies adjacency relations
                    //The name of the new variable should be read from right to left
                    HalfEdge<SPACE_DIM>* temp_next_twin_edge_from_A = edge_from_A->GetTwinHalfEdge()->GetNextHalfEdge();
                    HalfEdge<SPACE_DIM>* temp_twin_next_twin_edge_from_A = temp_next_twin_edge_from_A->GetTwinHalfEdge();

                    //Add a node to the intersected element first
                    HalfEdge<SPACE_DIM>* pNode_to_new_node_edge = p_element->AddNode(edge_from_A,new_node);
                    HalfEdge<SPACE_DIM>* twin_pNode_to_new_node_edge = pNode_to_new_node_edge->GetTwinHalfEdge();

                    assert(temp_twin_next_twin_edge_from_A->GetTargetNode()==pNode);
                    temp_twin_next_twin_edge_from_A->SetTargetNode(new_node);

                    //Adding the node does not yield correct adjacency relations in this case, so we fix it
                    assert(edge_from_A->GetTwinHalfEdge()->GetTargetNode()==new_node);
                    edge_from_A->GetTwinHalfEdge()->SetNextHalfEdge(temp_next_twin_edge_from_A, true);

                    //Insert twin_pNode_to_new_node_edge into the element neighbouring the intersecting element
                    twin_pNode_to_new_node_edge->SetNextHalfEdge(temp_twin_next_twin_edge_from_A->GetNextHalfEdge(),true);
                    twin_pNode_to_new_node_edge->SetPreviousHalfEdge(temp_twin_next_twin_edge_from_A, true);
                    twin_pNode_to_new_node_edge->SetElement(temp_twin_next_twin_edge_from_A->GetElement());

                    this->AddEdge(pNode_to_new_node_edge);
                    // Check the nodes are updated correctly
                    assert(pNode->GetContainingElementIndices().size() == 3);
                    assert(new_node->GetContainingElementIndices().size() == 2);
                }
                else if (num_common_vertices == 2)
                {
                    /*
                     *  From          To
                     *   _ B              _ B
                     *    |<---          |
                     *    | /|\          |\
                     *    |/ | \         | \
                     *    |  |  \        |\ \
                     *   _|__|___\      _|_\_\
                     *     A              A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--pNode--new_node--vertexB
                     * then vertexA is removed
                     */

                    //Reduce the situation to the case with one common vertex by deleting the common edge
                    HENode<SPACE_DIM>* deleted_node = edge_from_A->GetPreviousHalfEdge()->GetOriginNode();

                    std::set<unsigned int> deleted_node_element_indices = deleted_node->GetContainingElementIndices();
                    for (unsigned int index:deleted_node_element_indices)
                        affected_element_indices.insert(index);
                    assert(deleted_node_element_indices.size() >= 2);

                    node_A->rGetModifiableLocation() = deleted_node->rGetLocation();
                    std::set<HalfEdge<SPACE_DIM>* > deleted_edges = p_element->DeleteNode(deleted_node);
                    for (unsigned int index: affected_element_indices)
                    {
                        this->GetElement(index)->UpdateGeometry();
                    }

                    this->mDeletedNodeIndices.push_back(deleted_node->GetIndex());
                    deleted_node->MarkAsDeleted();
                    for (auto edges : deleted_edges)
                        this->MarkHalfEdgeAsDeleted(edges);

                    //Make sure we are in a situation with one common vertex
                    assert(pNode->GetNumContainingElements()==2);
                    num_common_vertices = 0;
                    for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                    {
                        for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                        {
                            if (p_element_common_1->GetNodeGlobalIndex(i) == p_element_common_2->GetNodeGlobalIndex(j))
                            {
                                num_common_vertices++;
                            }
                        }
                    }
                    assert(node_A->GetNumContainingElements()>2 || num_common_vertices==1);
                    assert(num_common_vertices != 2);

                    node_global_index = pNode->GetIndex();
                    num_nodes_elem_1 = p_element_1->GetNumNodes();
                    num_nodes_elem_2 = p_element_2->GetNumNodes();
                    local_index_1 = p_element_1->GetNodeLocalIndex(node_global_index);
                    next_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + 1)%num_nodes_elem_1);
                    previous_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + num_nodes_elem_1 - 1)%num_nodes_elem_1);

                    local_index_2 = p_element_2->GetNodeLocalIndex(node_global_index);
                    next_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + 1)%num_nodes_elem_2);
                    previous_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + num_nodes_elem_2 - 1)%num_nodes_elem_2);
                    const bool prev_case
                    = next_node_1 == vertexA_index || previous_node_1 == vertexA_index || next_node_2 == vertexA_index || previous_node_2 == vertexA_index;
                    assert(prev_case);

                    this->PerformT3Swap(pNode, elementIndex);
                }
                else
                {
                    // This can't happen as nodes can't be on the internal edge of two elements
                    NEVER_REACHED;
                }
            }
            else if (next_node_1 == vertexB_index || previous_node_1 == vertexB_index || next_node_2 == vertexB_index || previous_node_2 == vertexB_index)
            {
                // Get elements containing vertexB_index (the common vertex)

                assert(node_B->GetContainingElementIndices().size()>1);

                std::set<unsigned> elements_containing_vertex_B = node_B->GetContainingElementIndices();
                std::set<unsigned>::const_iterator iter = elements_containing_vertex_B.begin();
                HEElement<SPACE_DIM>* p_element_common_1 = this->GetElement(*iter);
                iter++;
                HEElement<SPACE_DIM>* p_element_common_2 = this->GetElement(*iter);

                // Calculate the number of common vertices between element_1 and element_2
                unsigned num_common_vertices = 0;
                for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                {
                    for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                    {
                        if (p_element_common_1->GetNodeGlobalIndex(i) == p_element_common_2->GetNodeGlobalIndex(j))
                        {
                            num_common_vertices++;
                        }
                    }
                }

                if (num_common_vertices == 1 || node_B->GetContainingElementIndices().size() > 2)
                {
                    /*
                     *  From          To
                     *   _B_________      _B____
                     *    |\   |   /       | / /
                     *    | \  |  /        |/ /
                     *    |  \ | /         | /
                     *    |   \|/          |/
                     *   _|   <---        _|
                     *    A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
                     */

                    // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                    intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(edge_from_A, intersection);
                    edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                    // Move original node and change to non-boundary node
                    pNode->rGetModifiableLocation() = intersection + 0.5*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;
                    pNode->SetAsBoundaryNode(false);

                    //Merge A--pNode edge first
                    HalfEdge<SPACE_DIM>* out_edge = node_B->GetOutgoingEdge();
                    HalfEdge<SPACE_DIM>* next_out_edge = out_edge;
                    HalfEdge<SPACE_DIM>* edge_to_pNode = nullptr;
                    do
                    {
                        if (next_out_edge->GetTargetNode()==pNode)
                        {
                            edge_to_pNode = next_out_edge;
                            break;
                        }
                        next_out_edge = next_out_edge->GetTwinHalfEdge()->GetNextHalfEdge();
                    }while(next_out_edge != out_edge);
                    assert(edge_to_pNode);
                    assert(edge_to_pNode->GetNextHalfEdge()->GetElement());
                    assert(!edge_to_pNode->GetTwinHalfEdge()->GetElement());
                    assert(edge_to_pNode->GetElement());

                    MergeEdgesInT3Swap(edge_from_A, edge_to_pNode, pNode);

                    // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                    c_vector<double, SPACE_DIM> new_node_location;
                    new_node_location = intersection - 0.5*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;

                    // Add new node, which will always be a boundary node
                    HENode<SPACE_DIM>* new_node = new HENode<SPACE_DIM>(0, true, new_node_location[0], new_node_location[1]);
                    AddNode(new_node);

                    //We need to store the edges below, because inserting the node modifies adjacency relations
                    //The name of the new variable should be read from right to left
                    HalfEdge<SPACE_DIM>* temp_prev_twin_edge_from_A = edge_from_A->GetTwinHalfEdge()->GetPreviousHalfEdge();
                    HalfEdge<SPACE_DIM>* temp_twin_prev_twin_edge_from_A = temp_prev_twin_edge_from_A->GetTwinHalfEdge();

                    //Add a node to the intersected element first
                    HalfEdge<SPACE_DIM>* A_to_new_node_edge = p_element->AddNode(edge_from_A,new_node);
                    assert(new_node->GetOutgoingEdge()->GetTargetNode()==node_A);
                    //HalfEdge<SPACE_DIM>* twin_A_to_new_node_edge = A_to_new_node_edge->GetTwinHalfEdge();

                    std::cout<<"NEW"<<std::endl;
                    std::cout<<"AB: "<<edge_from_A->GetNodePair().first<<" "<<edge_from_A->GetNodePair().second<<std::endl;
                    std::cout<<"to P: "<<A_to_new_node_edge->GetNodePair().first<<" "<<A_to_new_node_edge->GetNodePair().second<<std::endl;

                    //Check if the edge has been created correctly
                    assert(A_to_new_node_edge->GetNextHalfEdge()==edge_from_A);

                    //Check if the half edge points to the wrong node and change accordingly
                    assert(temp_prev_twin_edge_from_A->GetTargetNode()==pNode);
                    assert(temp_twin_prev_twin_edge_from_A->GetPreviousHalfEdge()->GetTargetNode()==pNode);
                    temp_prev_twin_edge_from_A->SetTargetNode(new_node);

                    //Adding the node does not yield correct adjacency relations in this case, so we fix it
                    A_to_new_node_edge->GetTwinHalfEdge()->SetPreviousHalfEdge(temp_prev_twin_edge_from_A, true);

                    //Insert edge_from_A into the element neighbouring the intersecting element
                    edge_from_A->GetTwinHalfEdge()->SetPreviousHalfEdge(temp_twin_prev_twin_edge_from_A->GetPreviousHalfEdge(), true);
                    edge_from_A->GetTwinHalfEdge()->SetNextHalfEdge(temp_twin_prev_twin_edge_from_A,true);
                    edge_from_A->GetTwinHalfEdge()->SetElement(temp_twin_prev_twin_edge_from_A->GetElement());

                    std::cout<<"from P: "<<pNode->GetOutgoingEdge()->GetNodePair().first<<" "<<pNode->GetOutgoingEdge()->GetNodePair().second<<std::endl;
                    HalfEdge<SPACE_DIM>* e1 =  pNode->GetOutgoingEdge()->GetTwinHalfEdge()->GetNextHalfEdge();
                    HalfEdge<SPACE_DIM>* e2 =  e1->GetTwinHalfEdge()->GetNextHalfEdge();
                    std::cout<<"from P 1 : "<<e1->GetNodePair().first<<" "<<e1->GetNodePair().second<<std::endl;
                    std::cout<<"from P 2 : "<<e2->GetNodePair().first<<" "<<e2->GetNodePair().second<<std::endl;
                    this->AddEdge(A_to_new_node_edge);
                    // Check the nodes are updated correctly
                    assert(pNode->GetContainingElementIndices().size() == 3);
                    assert(new_node->GetContainingElementIndices().size() == 2);
                }
                else if (num_common_vertices == 2)
                {
                    /*
                     *  From          To
                     *   _B_______      _B____
                     *    |  |   /       | / /
                     *    |  |  /        |/ /
                     *    |\ | /         | /
                     *    | \|/          |/
                     *   _| <---        _|
                     *    A
                     *
                     * The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
                     * then vertexB is removed
                     */

                    //Reduce the situation to the case with one common vertex by deleting the common edge
                    HENode<SPACE_DIM>* deleted_node = edge_from_A->GetNextHalfEdge()->GetTargetNode();

                    std::set<unsigned int> deleted_node_element_indices = deleted_node->GetContainingElementIndices();
                    for (unsigned int index:deleted_node_element_indices)
                        affected_element_indices.insert(index);
                    assert(deleted_node_element_indices.size() >= 2);

                    node_B->rGetModifiableLocation() = deleted_node->rGetLocation();
                    std::set<HalfEdge<SPACE_DIM>* > deleted_edges = p_element->DeleteNode(deleted_node);
                    for (unsigned int index: affected_element_indices)
                    {
                        this->GetElement(index)->UpdateGeometry();
                    }

                    this->mDeletedNodeIndices.push_back(deleted_node->GetIndex());
                    deleted_node->MarkAsDeleted();
                    for (auto edges : deleted_edges)
                        this->MarkHalfEdgeAsDeleted(edges);

                    //Make sure we are in a situation with one common vertex
                    assert(pNode->GetNumContainingElements()==2);
                    num_common_vertices = 0;
                    for (unsigned i=0; i<p_element_common_1->GetNumNodes(); i++)
                    {
                        for (unsigned j=0; j<p_element_common_2->GetNumNodes(); j++)
                        {
                            if (p_element_common_1->GetNodeGlobalIndex(i) == p_element_common_2->GetNodeGlobalIndex(j))
                            {
                                num_common_vertices++;
                            }
                        }
                    }
                    assert(node_B->GetNumContainingElements()>2 || num_common_vertices==1);
                    assert(num_common_vertices != 2);

                    node_global_index = pNode->GetIndex();
                    num_nodes_elem_1 = p_element_1->GetNumNodes();
                    num_nodes_elem_2 = p_element_2->GetNumNodes();
                    local_index_1 = p_element_1->GetNodeLocalIndex(node_global_index);
                    next_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + 1)%num_nodes_elem_1);
                    previous_node_1 = p_element_1->GetNodeGlobalIndex((local_index_1 + num_nodes_elem_1 - 1)%num_nodes_elem_1);

                    local_index_2 = p_element_2->GetNodeLocalIndex(node_global_index);
                    next_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + 1)%num_nodes_elem_2);
                    previous_node_2 = p_element_2->GetNodeGlobalIndex((local_index_2 + num_nodes_elem_2 - 1)%num_nodes_elem_2);
                    const bool prev_case
                    = next_node_1 == vertexB_index || previous_node_1 == vertexB_index || next_node_2 == vertexB_index || previous_node_2 == vertexB_index;
                    assert(prev_case);

                    this->PerformT3Swap(pNode, elementIndex);
                }
                else
                {
                    // This can't happen as nodes can't be on the internal edge of two elements
                    NEVER_REACHED;
                }
            }
            else
            {
                /*
                 *  From          To
                 *   _____         _______
                 *                  / | \
                 *    /|\   ^      /  |  \
                 *   / | \  |
                 *
                 * The edge goes from vertexA--vertexB to vertexA--new_node_1--pNode--new_node_2--vertexB
                 */
                // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                intersection = this->WidenEdgeOrCorrectIntersectionLocationIfNecessary(edge_from_A, intersection);
                edge_ab_unit_vector = this->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                // Move original node and change to non-boundary node
                pNode->rGetModifiableLocation() = intersection;
                pNode->SetAsBoundaryNode(false);

                c_vector<double, SPACE_DIM> new_node_1_location;
                new_node_1_location = intersection - this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;
                c_vector<double, SPACE_DIM> new_node_2_location;
                new_node_2_location = intersection + this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;

                // Add new nodes which will always be boundary nodes
                HENode<SPACE_DIM>* new_node_1 = new HENode<SPACE_DIM>(0, true, new_node_1_location[0], new_node_1_location[1]);
                HENode<SPACE_DIM>* new_node_2 = new HENode<SPACE_DIM>(0, true, new_node_2_location[0], new_node_2_location[1]);
                AddNode(new_node_1);
                AddNode(new_node_2);

                //Find an edge incoming pointing to pNode that is on the boundary, i.e. does
                //NOT have an element associated to it. There should be only one such edge
                HalfEdge<SPACE_DIM>* in_edge = pNode->GetOutgoingEdge()->GetTwinHalfEdge();
                HalfEdge<SPACE_DIM>* next_in_edge = in_edge;
                //HalfEdge<SPACE_DIM>* edge_to_pNode = nullptr;
                do
                {
                    if (!next_in_edge->GetElement())
                        break;
                    next_in_edge = next_in_edge->GetNextHalfEdge()->GetTwinHalfEdge();
                }while(next_in_edge != in_edge && next_in_edge->GetElement());
                if (next_in_edge == in_edge && next_in_edge->GetElement())
                {
                    EXCEPTION("All edges pointing to the intersecting node during T3 Swap are internal");
                }

                //We need to store the edges below, because inserting the new nodes modifies adjacency relations
                //The name of the new variable should be read from right to left

                //Edges of the first intersecting element, depicted on the left of the diagram above
                //next_in_edge is the edge pointing to pNode
                HalfEdge<SPACE_DIM>* temp_twin_in_edge = next_in_edge->GetTwinHalfEdge();
                HalfEdge<SPACE_DIM>* temp_prev_twin_in_edge = temp_twin_in_edge->GetPreviousHalfEdge();

                //Edges of the second intersecting element, depicted on the right of the diagram above
                HalfEdge<SPACE_DIM>* temp_next_in_edge = next_in_edge->GetNextHalfEdge();
                HalfEdge<SPACE_DIM>* temp_twin_next_in_edge = temp_next_in_edge->GetTwinHalfEdge();
                HalfEdge<SPACE_DIM>* temp_next_twin_next_in_edge = temp_twin_next_in_edge->GetNextHalfEdge();
                assert(!temp_next_in_edge->GetElement());
                assert(temp_prev_twin_in_edge->GetTwinHalfEdge()==temp_next_twin_next_in_edge);

                //Add the nodes to the intersected element
                HalfEdge<SPACE_DIM>* A_to_new_node_1_edge = p_element->AddNode(edge_from_A, new_node_1);
                this->AddEdge(A_to_new_node_1_edge);
                assert(A_to_new_node_1_edge->GetTargetNode()==new_node_1);

                HalfEdge<SPACE_DIM>* new_node_1_to_pNode = p_element->AddNode(edge_from_A, pNode);
                this->AddEdge(new_node_1_to_pNode);
                assert(new_node_1_to_pNode->GetTargetNode()==pNode);

                HalfEdge<SPACE_DIM>* pNode_to_new_node_2 = p_element->AddNode(edge_from_A, new_node_2);
                this->AddEdge(pNode_to_new_node_2);
                assert(pNode_to_new_node_2->GetTargetNode()==new_node_2);

                A_to_new_node_1_edge->GetTwinHalfEdge()->SetPreviousHalfEdge(next_in_edge,true);
                next_in_edge->SetTargetNode(new_node_1);

                //Adjacency relations of the "left" element.
                HalfEdge<SPACE_DIM>* twin_new_node_1_to_pNode = new_node_1_to_pNode->GetTwinHalfEdge();
                twin_new_node_1_to_pNode->SetElement(temp_twin_in_edge->GetElement());
                twin_new_node_1_to_pNode->SetPreviousHalfEdge(temp_prev_twin_in_edge,true);
                twin_new_node_1_to_pNode->SetNextHalfEdge(temp_twin_in_edge, true);

                //Adjacency relations of the "right" element.
                HalfEdge<SPACE_DIM>* twin_pNode_to_new_node_2 = pNode_to_new_node_2->GetTwinHalfEdge();
                twin_pNode_to_new_node_2->SetElement(temp_twin_next_in_edge->GetElement());
                twin_pNode_to_new_node_2->SetPreviousHalfEdge(temp_twin_next_in_edge,true);
                twin_pNode_to_new_node_2->SetNextHalfEdge(temp_next_twin_next_in_edge,true);

                edge_from_A->GetTwinHalfEdge()->SetNextHalfEdge(temp_next_in_edge,true);

                // Check the nodes are updated correctly
                assert(pNode->GetContainingElementIndices().size() == 3);
                assert(new_node_1->GetContainingElementIndices().size() == 2);
                assert(new_node_2->GetContainingElementIndices().size() == 2);
            }
        }
    }
    else
    {
        EXCEPTION("Trying to merge a node, contained in more than 2 elements, into another element, this is not possible with the vertex mesh.");
    }

    for (unsigned int index: affected_element_indices)
    {
        HEElement<SPACE_DIM>* update_element = this->GetElement(index);
        update_element->UpdateGeometry();
    }

}

template <unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::PerformVoidRemoval(HalfEdge<SPACE_DIM>* pEdge)
{
    // Note that we define this vector before setting it, as otherwise the profiling build will break (see #2367)
    c_vector<double, SPACE_DIM> node_0_location = pEdge->GetTargetNode()->rGetLocation();
    c_vector<double, SPACE_DIM> node_1_location = pEdge->GetNextHalfEdge()->GetTargetNode()->rGetLocation();
    c_vector<double, SPACE_DIM> node_2_location = pEdge->GetOriginNode()->rGetLocation();

    // Calculate void centroid
    c_vector<double, SPACE_DIM> nodes_midpoint = node_0_location
                        + this->GetVectorFromAtoB(node_0_location, node_1_location) / 3.0
                        + this->GetVectorFromAtoB(node_0_location, node_2_location) / 3.0;

    HENode<SPACE_DIM>* p_new_node = RemoveTriangle(pEdge);
    p_new_node->rGetModifiableLocation() = nodes_midpoint;

    // Remove the deleted nodes/edges and re-index
    RemoveDeletedNodes();
    RemoveDeletedEdges();
}

template <unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HEMutableVertexMesh<SPACE_DIM>::RemoveTriangle(HalfEdge<SPACE_DIM>* pEdge)
{
    // Note that we define this vector before setting it, as otherwise the profiling build will break (see #2367)
    c_vector<double, SPACE_DIM> new_node_location;
    new_node_location = zero_vector<double>(SPACE_DIM);

    // Create a new node at the element's centroid. By default, we set the boundary status to false;
    HENode<SPACE_DIM>* p_new_node = new HENode<SPACE_DIM>(GetNumNodes(), new_node_location, false);
    AddNode(p_new_node);

    //Set target vertices of the incoming edges to the new node
    //without changing neighbour relations of edges (next or previous) so that the case with rosettes is resolved
    //correctly (i.e. if the vertex is contained in more than 3 elements)

    //Loop over the nodes
    HalfEdge<SPACE_DIM>* next_edge = pEdge;
    HEElement<SPACE_DIM>* element = pEdge->GetElement();
    do
    {
        //Loop over incoming edges
        HalfEdge<SPACE_DIM>* next_in_edge = next_edge;
        do
        {
            //Only change target vertices of the edges that are not contained in this element
            bool edge_in_this_element
            = next_in_edge->GetElement()==element||next_in_edge->GetTwinHalfEdge()->GetElement()==element;
            if (!edge_in_this_element)
                next_in_edge->SetTargetNode(p_new_node);
            next_in_edge = next_in_edge->GetNextHalfEdge()->GetTwinHalfEdge();
        }while(next_in_edge != next_edge);
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != pEdge);

    //Change adjacency relations of the neighbouring elements' edges
    next_edge = pEdge;
    std::vector<unsigned int> neighbouring_element_indices;
    do
    {
        HalfEdge<SPACE_DIM>* twin = next_edge->GetTwinHalfEdge();
        if (twin->GetElement())
        {
            if (twin->GetElement()->GetNumNodes()<4)
            {
                EXCEPTION("One of the neighbours of a small triangular element is also a triangle - dealing with this has not been implemented yet");
            }
            neighbouring_element_indices.push_back(twin->GetElement()->GetIndex());
            if (twin->GetElement()->GetHalfEdge()==twin)
            {
                twin->GetElement()->SetHalfEdge(twin->GetNextHalfEdge());
            }
        }
        HalfEdge<SPACE_DIM>* twin_previous = twin->GetPreviousHalfEdge();
        HalfEdge<SPACE_DIM>* twin_next = twin->GetNextHalfEdge();

        //Edge is now inaccessible (effectively deleted) from the neighbouring element
        twin_previous->SetNextHalfEdge(twin_next,true);
        p_new_node->SetOutgoingEdge(twin_next);

        HENode<SPACE_DIM>* target_node = next_edge->GetTargetNode();

        //Mark element vertex and edges as deleted
        this->mDeletedNodeIndices.push_back(target_node->GetIndex());
        target_node->MarkAsDeleted();
        this->MarkHalfEdgeAsDeleted(next_edge);

        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != pEdge);

    //Update neighbouring elements
    for (unsigned int index : neighbouring_element_indices)
    {
        this->GetElement(index)->UpdateGeometry();
    }
    return p_new_node;
}

template <unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::CheckForRosettes()
{

}

template<unsigned SPACE_DIM>
c_vector<double, 2>
HEMutableVertexMesh<SPACE_DIM>::WidenEdgeOrCorrectIntersectionLocationIfNecessary(HalfEdge<SPACE_DIM>* pEdge, c_vector<double,2> intersection)
{
    /**
     * If the edge is shorter than 4.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold move vertexA and vertexB
     * 4.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold apart.
     * \todo investigate if moving A and B causes other issues with nearby nodes (see #2401)
     *
     * Note: this distance is so that there is always enough room for new nodes (if necessary)
     * \todo currently this assumes a worst case scenario of 3 nodes between A and B could be less movement for other cases
     *       (see #1399 and #2401)
     */
    HENode<SPACE_DIM>* nodeA = pEdge->GetOriginNode();
    HENode<SPACE_DIM>* nodeB = pEdge->GetTargetNode();
    c_vector<double, SPACE_DIM> vertexA = pEdge->GetOriginNode()->rGetLocation();
    c_vector<double, SPACE_DIM> vertexB = pEdge->GetTargetNode()->rGetLocation();
    c_vector<double, SPACE_DIM> vector_a_to_b = pEdge->GetVector();

    if (pEdge->GetLength() < 4.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold)
    {
        WARNING("Trying to merge a node onto an edge which is too small.");

        c_vector<double, SPACE_DIM> centre_a_and_b = vertexA + 0.5*vector_a_to_b;

        vertexA = centre_a_and_b  - 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*vector_a_to_b/norm_2(vector_a_to_b);
        ChastePoint<SPACE_DIM> vertex_A_point(vertexA);
        nodeA->SetPoint(vertex_A_point);

        vertexB = centre_a_and_b  + 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*vector_a_to_b/norm_2(vector_a_to_b);
        ChastePoint<SPACE_DIM> vertex_B_point(vertexB);
        nodeB->SetPoint(vertex_B_point);

        intersection = centre_a_and_b;
        pEdge->ComputeLength();
    }

    // Reset distances
    vector_a_to_b = pEdge->GetVector();
    c_vector<double,2> edge_ab_unit_vector = vector_a_to_b/pEdge->GetLength();

    // Reset the intersection away from vertices A and B to allow enough room for new nodes
    /**
     * If the intersection is within this->mCellRearrangementRatio^2*this->mCellRearrangementThreshold of vertexA or vertexB move it
     * this->mCellRearrangementRatio^2*this->mCellRearrangementThreshold away.
     *
     * Note: this distance so that there is always enough room for new nodes (if necessary).
     * \todo currently this assumes a worst case scenario of 3 nodes between A and B; could be less movement for other cases
     *       (see #2401)
     */
    if (norm_2(intersection - vertexA) < 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold)
    {
        intersection = vertexA + 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;
    }
    if (norm_2(intersection - vertexB) < 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold)
    {
        intersection = vertexB - 2.0*this->mCellRearrangementRatio*this->mCellRearrangementThreshold*edge_ab_unit_vector;
    }
    return intersection;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::MarkHalfEdgeAsDeleted(HalfEdge<SPACE_DIM>* half_edge)
{
    this->mDeletedHalfEdges.push_back(half_edge);
    this->mDeletedHalfEdges.push_back(half_edge->GetTwinHalfEdge());
    half_edge->SetDeletedStatus(true, true);
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}

// Explicit instantiation
template class HEMutableVertexMesh<1>;
template class HEMutableVertexMesh<2>;
template class HEMutableVertexMesh<3>;
