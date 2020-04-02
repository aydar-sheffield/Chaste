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
unsigned HEMutableVertexMesh<SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - this->mDeletedElementIndices.size();
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
    return false;
}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::IdentifySwapType(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB)
{

}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::HandleHighOrderJunctions(HENode<SPACE_DIM>* pNodeA, HENode<SPACE_DIM>* pNodeB)
{

}

template<unsigned SPACE_DIM>
void HEMutableVertexMesh<SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{

}

// Explicit instantiation
template class HEMutableVertexMesh<1>;
template class HEMutableVertexMesh<2>;
template class HEMutableVertexMesh<3>;
