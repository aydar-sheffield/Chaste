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
                                                    cellRearrangementRatio,
                                                    t2Threshold,
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

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - this->mDeletedNodeIndices.size();
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumFullEdges() const
{
    return this->mFullEdges.size() - mDeletedHalfEdges.size()/2;
}

template<unsigned SPACE_DIM>
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - this->mDeletedElementIndices.size();
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> HEMutableVertexMesh<SPACE_DIM>::GetLastT2SwapLocation()
{
    return mLastT2SwapLocation;
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
    return new_element_index;
}

template<unsigned SPACE_DIM>
bool HEMutableVertexMesh<SPACE_DIM>::CheckForSwapsFromShortEdges()
{
    for (FullEdge<SPACE_DIM>* full_edge : this->mFullEdges)
    {
        //Assume edge lengths have been updated
        if (full_edge->GetLength() < this->mCellRearrangementThreshold)
        {
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

            if (!edge_share_triangular_element)
            {
                IdentifySwapType(*full_edge);
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

        }
    }
}

template<unsigned int SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::PerformNodeMerge(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB)
{

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

    previous_edge->SetNextHalfEdge(next_edge,true);
    twin_previous->SetNextHalfEdge(twin_next,true);

    this->mDeletedNodeIndices.push_back(pNodeB->GetIndex());
    pNodeB->MarkAsDeleted();
    this->mDeletedHalfEdges.push_back(pEdge);
    this->mDeletedHalfEdges.push_back(pEdge->GetTwinHalfEdge());
    next_edge->SetDeletedStatus(true, true);

    if (pEdge->GetElement())
    {
        pEdge->GetElement()->SetHalfEdge(next_edge);
    }
    if (twin->GetElement())
    {
        twin->GetElement()->SetHalfEdge(twin_next);
    }
    //TODO:Is updating adjacent elements necessary? Volume, or node numbers, for example
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::HandleHighOrderJunctions(FullEdge<SPACE_DIM>& full_edge)
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
        HalfEdge<SPACE_DIM>* edge = this->mFullEdges[i]->operator()(0);
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

    // Create a new node at the element's centroid; this will be a boundary node if any existing nodes were on the boundary
    bool is_node_on_boundary = false;
    HENode<SPACE_DIM>* p_new_node = new HENode<SPACE_DIM>(GetNumNodes(), new_node_location, is_node_on_boundary);
    unsigned new_node_global_index = AddNode(p_new_node);

    //Set target vertices of the incoming edges to the new node
    //without changing neighbour relations of edges (next or previous) so that the case with rosettes is resolved
    //correctly (i.e. if the vertex is contained in more than 3 elements)

    //Loop over the nodes
    HalfEdge<SPACE_DIM>* edge = rElement.GetHalfEdge();
    HalfEdge<SPACE_DIM>* next_edge = edge;
    do
    {
        //Loop over incoming edges
        HalfEdge<SPACE_DIM>* next_in_edge = next_edge;
        do
        {
            //Only change target vertices of the edges that are not contained in this element
            bool edge_in_this_element
            = next_in_edge->GetElement()==&rElement||next_in_edge->GetTwinHalfEdge()->GetElement()==&rElement;
            if (!edge_in_this_element)
                next_in_edge->SetTargetNode(p_new_node);
            next_in_edge = next_in_edge->GetNextHalfEdge()->GetTwinHalfEdge();
        }while(next_in_edge != next_edge);
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != edge);

    //Change adjacency relations of the neighbouring elements' edges
    next_edge = edge;
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
        //If any of the nodes are on the boundary, the new node will also be on the boundary
        if (target_node->IsBoundaryNode())
        {
            is_node_on_boundary = true;
        }
        //Mark element vertex and edges as deleted
        this->mDeletedNodeIndices.push_back(target_node->GetIndex());
        target_node->MarkAsDeleted();
        this->mDeletedHalfEdges.push_back(next_edge);
        this->mDeletedHalfEdges.push_back(next_edge->GetTwinHalfEdge());
        next_edge->SetDeletedStatus(true, true);

        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != edge);

    p_new_node->SetAsBoundaryNode(is_node_on_boundary);

    this->mDeletedElementIndices.push_back(rElement.GetIndex());
    rElement.MarkAsDeleted();

    //Update neighbouring elements
    for (unsigned int index : neighbouring_element_indices)
    {
        this->GetElement(index)->RegisterWithHalfEdges();
        this->GetElement(index)->ComputeSurfaceArea();
        this->GetElement(index)->ComputeVolume();
    }
}


// Explicit instantiation
template class HEMutableVertexMesh<1>;
template class HEMutableVertexMesh<2>;
template class HEMutableVertexMesh<3>;
