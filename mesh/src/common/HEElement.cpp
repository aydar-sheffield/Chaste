/*
 * HEElement.cpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#include "HEElement.hpp"
#include <map>
template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::CommonConstructor(const std::vector<HENode<SPACE_DIM>* > node_list)
{
    const unsigned int n_nodes = node_list.size();

    std::map<unsigned int, std::pair<HalfEdge<SPACE_DIM>*, HalfEdge<SPACE_DIM>*> > common_node_to_in_out_edges;
    for (unsigned int i=0; i<n_nodes; ++i)
    {
        HalfEdge<SPACE_DIM>* out_edge = node_list[i]->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* edge;
        const unsigned int next_index = (i+1)%n_nodes;
        const unsigned int previous_index = (i-1+n_nodes)%n_nodes;
        /**
         * If neighbouring halfedges are constructed, get the edge outgoing from this node
         */
        bool construct_edge = true;
        //If the node is a common to another edge
        if (out_edge)
        {
            //Iterate over outgoing edges to check if the halfedge has already been constructed
            typename HENode<SPACE_DIM>::OutgoingEdgeIterator iter = node_list[i]->GetOutgoingEdgeIteratorBegin();
            typename HENode<SPACE_DIM>::OutgoingEdgeIterator end_iter = node_list[i]->GetOutgoingEdgeIteratorEnd();

            for (; iter!=end_iter; ++iter)
            {
                bool is_edge_constructed = iter->GetTargetNode() == node_list[next_index];
                if (is_edge_constructed)
                {
                    edge = *iter;
                    construct_edge =false;
                    break;
                }

            }

            //If a new outgoing edge needs to be constructed, check to see if an existing outgoing edge
            //is shared with this element. If it is not, this element must only have one common vertex with at least one element.
            //Therefore, an external outgoing edge of this element must be the next edge of another element's external edge pointing to this node
            //I.e. we have:
            // \ V /
            //  \ /
            //   A  this
            //  / \
            //_/ V \
            //where V means void and A is the common node
            //Note that there could only be one external incoming or outgoing edge
            iter = node_list[i]->GetOutgoingEdgeIteratorBegin();
            bool is_out_edge_common = false;
            for (; iter!=end_iter; ++iter)
            {
                if (iter->GetTargetNode() == node_list[previous_index])
                {
                    is_out_edge_common = true;
                    break;
                }
            }
            if (construct_edge&&!is_out_edge_common)
            {
                iter = node_list[i]->GetOutgoingEdgeIteratorBegin();
                std::pair<HalfEdge<SPACE_DIM>*, HalfEdge<SPACE_DIM>* > in_out_pair(nullptr, nullptr);
                for (; iter!=end_iter; ++iter)
                {
                    if (!iter->GetTwinHalfEdge()->GetElement())
                        in_out_pair.first = iter->GetTwinHalfEdge();
                    if (!iter->GetElement())
                        in_out_pair.second = *iter;
                }
                common_node_to_in_out_edges[i] = in_out_pair;
            }

        }

        if (construct_edge)
        {
            edge = new HalfEdge<SPACE_DIM>();

            edge->SetTargetNode(node_list[next_index]);

            HalfEdge<SPACE_DIM>* twin_edge = new HalfEdge<SPACE_DIM>();
            twin_edge->SetTargetNode(node_list[i]);
            twin_edge->SetTwinHalfEdge(edge,true);
        }
        node_list[i]->SetOutgoingEdge(edge);
        edge->SetElement(this);
    }


    for (unsigned int i=0; i<n_nodes; ++i)
    {
        //Link internal edges
        HalfEdge<SPACE_DIM>* out_edge = node_list[i]->GetOutgoingEdge();
        const unsigned int next_index = (i+1)%n_nodes;
        const unsigned int previous_index = (i-1+n_nodes)%n_nodes;
        HalfEdge<SPACE_DIM>* next_edge = node_list[next_index]->GetOutgoingEdge();
        HalfEdge<SPACE_DIM>* previous_edge = node_list[previous_index]->GetOutgoingEdge();
        out_edge->SetNextHalfEdge(next_edge);
        out_edge->SetPreviousHalfEdge(previous_edge);

        //Linking external (CW or boundary) edges
        HalfEdge<SPACE_DIM>* out_edge_twin = out_edge->GetTwinHalfEdge();
        HalfEdge<SPACE_DIM>* next_edge_twin = next_edge -> GetTwinHalfEdge();
        HalfEdge<SPACE_DIM>* previous_edge_twin = previous_edge -> GetTwinHalfEdge();

        //If on boundary
        if (!out_edge_twin->GetElement())
        {
            while(previous_edge_twin->GetElement())
                previous_edge_twin = previous_edge_twin->GetPreviousHalfEdge()->GetTwinHalfEdge();
            out_edge_twin->SetNextHalfEdge(previous_edge_twin, true);

            while(next_edge_twin->GetElement())
                next_edge_twin = next_edge_twin->GetNextHalfEdge()->GetTwinHalfEdge();
            out_edge_twin->SetPreviousHalfEdge(next_edge_twin, true);
        }
    }
    for (unsigned int i=0; i<n_nodes; ++i)
    {
        //If there is a common node (but not edge) between elements
        if (common_node_to_in_out_edges.count(i)>0)
        {
            HalfEdge<SPACE_DIM>* out_edge_twin = node_list[i]->GetOutgoingEdge()->GetTwinHalfEdge();
            HalfEdge<SPACE_DIM>* previous_edge_twin = node_list[i]->GetOutgoingEdge()->GetPreviousHalfEdge() -> GetTwinHalfEdge();

            if (common_node_to_in_out_edges[i].second)
                out_edge_twin->SetNextHalfEdge(common_node_to_in_out_edges[i].second,true);
            if (common_node_to_in_out_edges[i].first)
                previous_edge_twin->SetPreviousHalfEdge(common_node_to_in_out_edges[i].first,true);
        }
    }

    mpHalfEdge = node_list[0]->GetOutgoingEdge();
    mNumNodes = n_nodes;
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index)
:AbstractElement<SPACE_DIM,SPACE_DIM>(index),
 mpHalfEdge(nullptr),
 mNumNodes(0),
 mVolume(0),
 mSurfaceArea(0)
{
    mpHalfEdge = nullptr;
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<HENode<SPACE_DIM>* > node_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index),
 mpHalfEdge(nullptr),
 mNumNodes(0),
 mVolume(0),
 mSurfaceArea(0)
{
    CommonConstructor(node_list);
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, HalfEdge<SPACE_DIM>* edge)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index),
mpHalfEdge(edge)
{
    RegisterWithHalfEdges();
}

template <unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(unsigned int index, const std::vector<Node<SPACE_DIM>* > node_list)
:AbstractElement<SPACE_DIM, SPACE_DIM>(index),
 mpHalfEdge(nullptr),
 mNumNodes(0),
 mVolume(0),
 mSurfaceArea(0)
 {
    std::vector<HENode<SPACE_DIM>* > henode_list;
    for(auto node:node_list)
    {
        henode_list.push_back(new HENode<SPACE_DIM>(*node));
    }
    mpHalfEdge = nullptr;
    CommonConstructor(henode_list);
 }

template <unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::HEElement(const VertexElement<SPACE_DIM, SPACE_DIM> &rElement)
:AbstractElement<SPACE_DIM, SPACE_DIM>(rElement.GetIndex()),
 mpHalfEdge(nullptr),
 mNumNodes(0),
 mVolume(0),
 mSurfaceArea(0)
 {
    std::vector<HENode<SPACE_DIM>* > node_list;
    for (unsigned int i=0; i<rElement.GetNumNodes(); ++i)
    {
        node_list.push_back(new HENode<SPACE_DIM>(*rElement.GetNode(i)));
    }
    mpHalfEdge = nullptr;
    CommonConstructor(node_list);
 }

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>::~HEElement()
{}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::GetHalfEdge(const unsigned int local_index) const
{
    assert(local_index<mNumNodes);
    //Outgoing internal edge from node local_index
    HalfEdge<SPACE_DIM>* out_edge = mpHalfEdge;
    for (unsigned int i=0; i<local_index; ++i)
    {
        out_edge = out_edge->GetNextHalfEdge();
    }
    return out_edge;
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::GetHalfEdge(HENode<SPACE_DIM>* pTarget) const
{
    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
    if (next_edge->GetTargetNode()!=pTarget)
    {
        do
        {
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge&&pTarget!=next_edge->GetTargetNode());
        if (next_edge == mpHalfEdge)
        {
            EXCEPTION("Halfedge search by target node not found.");
        }
    }
    return next_edge;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::SetHalfEdge(HalfEdge<SPACE_DIM>* pEdge)
{
    mpHalfEdge = pEdge;
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HEElement<SPACE_DIM>::GetNode(unsigned int local_index) const
{
    assert(local_index<mNumNodes);
    HalfEdge<SPACE_DIM>* edge = GetHalfEdge(local_index);
    return edge->GetOriginNode();
}

template<unsigned int SPACE_DIM>
unsigned HEElement<SPACE_DIM>::GetNodeGlobalIndex(unsigned int local_index) const
{
    return GetNode(local_index)->GetIndex();
}

template<unsigned int SPACE_DIM>
unsigned HEElement<SPACE_DIM>::GetNodeLocalIndex(unsigned globalIndex) const
{
    unsigned local_index = UINT_MAX;
    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge->GetPreviousHalfEdge();
    unsigned int counter = 0;
    do
    {
        if (next_edge->GetTargetNode()->GetIndex() == globalIndex)
        {
            local_index = counter;
        }
        counter++;
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge != mpHalfEdge->GetPreviousHalfEdge());

    return local_index;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::UpdateNode(const unsigned& rIndex, HENode<SPACE_DIM>* pNode)
{
    assert(rIndex<mNumNodes);
    //Outgoing internal edge from node rIndex
    HalfEdge<SPACE_DIM>* out_edge = GetHalfEdge(rIndex);
    out_edge->SetOriginNode(pNode);
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::ReplaceNode(HENode<SPACE_DIM>* pOldNode, HENode<SPACE_DIM>* pNewNode)
{
    assert(pOldNode != pNewNode);

    HalfEdge<SPACE_DIM>* next_edge = GetHalfEdge(pOldNode)->GetNextHalfEdge();

    next_edge->SetOriginNode(pNewNode);
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::MarkAsDeleted()
{
    this->mIsDeleted = true;

    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;

    do
    {
        next_edge->SetElement(nullptr);
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge!=mpHalfEdge);

}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::RegisterWithHalfEdges()
{
    assert(mpHalfEdge);
    HalfEdge<SPACE_DIM>* edge = mpHalfEdge;
    mNumNodes = 0;
    do
    {
        mNumNodes++;
        edge->SetElement(this);
        edge = edge->GetNextHalfEdge();
    }while(edge!=mpHalfEdge);
}

template<unsigned int SPACE_DIM>
HalfEdge<SPACE_DIM>* HEElement<SPACE_DIM>::AddNode(HalfEdge<SPACE_DIM>* pEdge, HENode<SPACE_DIM>* pNode)
{
    HalfEdge<SPACE_DIM>* prev_edge = pEdge->GetPreviousHalfEdge();

    //A new halfedge between pEdge->OriginNode() and pNode is created with pNode as its target index.
    HalfEdge<SPACE_DIM>* new_edge = new HalfEdge<SPACE_DIM>(this);
    //Set the element to which the new twin edge belongs to
    HalfEdge<SPACE_DIM>* new_edge_twin = new HalfEdge<SPACE_DIM>(pEdge->GetTwinHalfEdge()->GetElement());

    //Set edge adjacency relations
    new_edge->SetTwinHalfEdge(new_edge_twin, true);
    new_edge->SetNextHalfEdge(pEdge, true);
    new_edge->SetPreviousHalfEdge(prev_edge, true);

    //External halfedges
    new_edge_twin->SetNextHalfEdge(pEdge->GetTwinHalfEdge()->GetNextHalfEdge(), true);
    new_edge_twin->SetPreviousHalfEdge(pEdge->GetTwinHalfEdge(), true);

    //Set target Node relations
    new_edge_twin->SetTargetNode(pEdge->GetOriginNode());
    if (pEdge->GetOriginNode()->GetOutgoingEdge()==pEdge)
    {
        pEdge->GetOriginNode()->SetOutgoingEdge(new_edge);
    }
    pEdge->SetOriginNode(pNode);
    pNode->SetOutgoingEdge(new_edge_twin);

    assert(new_edge->IsFullyInitialized());
    assert(new_edge_twin->IsFullyInitialized());

    if (pEdge==mpHalfEdge)
        mpHalfEdge = new_edge;
    mNumNodes++;
    return new_edge;
}

template<unsigned int SPACE_DIM>
std::set<HalfEdge<SPACE_DIM>* > HEElement<SPACE_DIM>::DeleteNode(HENode<SPACE_DIM>* pNode)
{
    //Deletion of a node in triangular element not supported
    assert(mNumNodes>3);

    std::set<HalfEdge<SPACE_DIM>* > deleted_edges;

    //Find incoming edge
    HalfEdge<SPACE_DIM>* in_edge = GetHalfEdge(pNode);
    HalfEdge<SPACE_DIM>* in_edge_twin = in_edge->GetTwinHalfEdge();

    HalfEdge<SPACE_DIM>* next_in_edge = in_edge;
    unsigned int num_in_edges= 0;
    do
    {
        num_in_edges++;
        next_in_edge = next_in_edge->GetNextHalfEdge()->GetTwinHalfEdge();
    }while(next_in_edge != in_edge);
    assert(num_in_edges>=2);

    //If the deleted node only connects two edges...
    if (num_in_edges==2)
    {
        in_edge->SetTargetNode(in_edge->GetNextHalfEdge()->GetTargetNode());
        HalfEdge<SPACE_DIM>* deleted_edge = in_edge->GetNextHalfEdge();
        deleted_edges.insert(deleted_edge);

        //If the deleted edge or its twin are associated with their respective elements
        //shift the edges to the next ones
        if (deleted_edge == mpHalfEdge)
        {
            mpHalfEdge = in_edge;
        }
        if (deleted_edge->GetTwinHalfEdge()->GetElement())
        {
            HEElement<SPACE_DIM>* neighbouring_element = deleted_edge->GetTwinHalfEdge()->GetElement();
            if (neighbouring_element->GetHalfEdge() == deleted_edge->GetTwinHalfEdge())
            {
                neighbouring_element->SetHalfEdge(deleted_edge->GetTwinHalfEdge()->GetNextHalfEdge());
            }
        }
        //If the twin of the deleted edge is an outgoing edge of its origin node, chage the outgoing edge of the node
        if (deleted_edge->GetTwinHalfEdge()==deleted_edge->GetTwinHalfEdge()->GetOriginNode()->GetOutgoingEdge())
            deleted_edge->GetTwinHalfEdge()->GetOriginNode()->SetOutgoingEdge(deleted_edge->GetTwinHalfEdge()->GetNextHalfEdge());

        in_edge->SetNextHalfEdge(deleted_edge->GetNextHalfEdge(), true);
        in_edge_twin->SetPreviousHalfEdge(in_edge_twin->GetPreviousHalfEdge()->GetPreviousHalfEdge(),true);
    }
    else
    {
        next_in_edge = in_edge;
        do
        {
            HalfEdge<SPACE_DIM>* temp_next_in_edge = next_in_edge->GetNextHalfEdge()->GetTwinHalfEdge();
            if (next_in_edge->GetElement())
            {
                //Next affected element is null if the next_in_edge->GetNExtHalfEdge() is on boundary
                HEElement<SPACE_DIM>* next_affected_element = temp_next_in_edge->GetElement();

                next_in_edge->SetTargetNode(next_in_edge->GetNextHalfEdge()->GetTargetNode());

                if (next_in_edge->GetNextHalfEdge() == next_in_edge->GetElement()->GetHalfEdge())
                {
                    next_in_edge->GetElement()->SetHalfEdge(next_in_edge);
                }

                next_in_edge->SetNextHalfEdge(next_in_edge->GetNextHalfEdge()->GetNextHalfEdge(), true);

                HalfEdge<SPACE_DIM>* twin_next_in_edge = next_in_edge->GetTwinHalfEdge();
                twin_next_in_edge->SetElement(nullptr);
                if (next_affected_element)
                {
                    twin_next_in_edge->SetPreviousHalfEdge(temp_next_in_edge->GetTwinHalfEdge(),true);
                }
                else
                {
                    twin_next_in_edge->SetPreviousHalfEdge(temp_next_in_edge->GetPreviousHalfEdge(), true);
                    //Marking internal edges
                    HalfEdge<SPACE_DIM>* deleted_edge = temp_next_in_edge->GetTwinHalfEdge();
                    deleted_edges.insert(deleted_edge);
                    if (deleted_edge->GetTwinHalfEdge()==deleted_edge->GetTwinHalfEdge()->GetOriginNode()->GetOutgoingEdge())
                        deleted_edge->GetTwinHalfEdge()->GetOriginNode()->SetOutgoingEdge(deleted_edge->GetTwinHalfEdge()->GetNextHalfEdge());
                }

                if (next_in_edge->GetElement() != this)
                    next_in_edge->GetElement()->GetNumNodes(true);
            }
            next_in_edge = temp_next_in_edge;
        }while(next_in_edge != in_edge);
    }
    mNumNodes--;
    return deleted_edges;
}

template<unsigned int SPACE_DIM>
unsigned int HEElement<SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned int SPACE_DIM>
unsigned int HEElement<SPACE_DIM>::GetNumNodes(const bool compute)
{
    if (compute)
    {
        HalfEdge<SPACE_DIM>* edge = mpHalfEdge;
        mNumNodes = 0;
        do
        {
            mNumNodes++;
            edge = edge->GetNextHalfEdge();
        }while(edge!=mpHalfEdge);
        return mNumNodes;
    }
    else
    {
        return mNumNodes;
    }
}

template<unsigned int SPACE_DIM>
bool HEElement<SPACE_DIM>::IsElementOnBoundary() const
{
    bool is_element_on_boundary = false;
    HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
    do
    {
        is_element_on_boundary = next_edge->GetTargetNode()->IsBoundaryNode();
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge!=mpHalfEdge &&!is_element_on_boundary);

    return is_element_on_boundary;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::GetVolume() const
{
    return mVolume;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::GetSurfaceArea() const
{
    return mSurfaceArea;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::ComputeVolume()
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    double element_volume = 0.0;
    if (SPACE_DIM == 2)
    {
        // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location;
        first_node_location = mpHalfEdge->GetOriginNode()->rGetLocation();
        c_vector<double, SPACE_DIM> pos_1;
        pos_1 = zero_vector<double>(SPACE_DIM);

        // Loop over vertices
        HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
        do
        {
            c_vector<double, SPACE_DIM> next_node_location = next_edge->GetTargetNode()->rGetLocation();
            c_vector<double, SPACE_DIM> pos_2 = next_node_location-first_node_location;

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            element_volume += 0.5 * (this_x * next_y - next_x * this_y);

            pos_1 = pos_2;
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge);
    }
    else
    {
        //3D case not supported
        EXCEPTION("Half-edge mesh in 3D not supported.");
    }
    // We take the absolute value just in case the nodes were really oriented clockwise
    mVolume = fabs(element_volume);

    return mVolume;
}

template<unsigned int SPACE_DIM>
double HEElement<SPACE_DIM>::ComputeSurfaceArea()
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    double surface_area = 0.0;
    if (SPACE_DIM == 2)
    {
        HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
        do
        {
            surface_area += next_edge->ComputeLength();
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge);
    }
    else
    {
        //3D case not supported
        EXCEPTION("Half-edge mesh in 3D not supported.");
    }
    mSurfaceArea = surface_area;
    return surface_area;
}

template<unsigned int SPACE_DIM>
void HEElement<SPACE_DIM>::UpdateGeometry()
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    double element_volume = 0.0;
    double surface_area = 0.0;
    mNumNodes = 0;
    if (SPACE_DIM == 2)
    {
        // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location;
        first_node_location = mpHalfEdge->GetOriginNode()->rGetLocation();
        c_vector<double, SPACE_DIM> pos_1;
        pos_1 = zero_vector<double>(SPACE_DIM);

        // Loop over vertices
        HalfEdge<SPACE_DIM>* next_edge = mpHalfEdge;
        do
        {
            c_vector<double, SPACE_DIM> next_node_location = next_edge->GetTargetNode()->rGetLocation();
            c_vector<double, SPACE_DIM> pos_2 = next_node_location-first_node_location;

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            element_volume += 0.5 * (this_x * next_y - next_x * this_y);

            pos_1 = pos_2;
            surface_area += next_edge->ComputeLength(); //also updates edge lengths
            mNumNodes++;

            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=mpHalfEdge);
    }
    else
    {
        //3D case not supported
        EXCEPTION("Half-edge mesh in 3D not supported.");
    }
    // We take the absolute value just in case the nodes were really oriented clockwise
    mVolume = fabs(element_volume);
    mSurfaceArea = surface_area;
}

template class HEElement<1>;
template class HEElement<2>;
template class HEElement<3>;

