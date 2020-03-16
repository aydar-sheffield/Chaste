/*
 * HEVertexMesh.cpp
 *
 *  Created on: 11 Mar 2020
 *      Author: aydar
 */

#include "HEVertexMesh.hpp"

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEVertexMesh(std::vector<HEVertex<SPACE_DIM>* > vertices,
                                      std::vector<HEElement<SPACE_DIM>* > elements)
{
    Clear();
    // Populate mNodes and mElements
    for (unsigned index = 0; index < vertices.size(); index++)
    {
        Node<SPACE_DIM>* p_temp_node = vertices[index];
        this->mNodes.push_back(p_temp_node);
    }

    for (unsigned index = 0; index < elements.size(); index++)
    {
        HEElement<SPACE_DIM>* p_temp_element = elements[index];
        mElements.push_back(p_temp_element);
    }


}

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                      std::vector<VertexElement<SPACE_DIM, SPACE_DIM>*> vertexElements)
{
    Clear();

    std::set<HEVertex<SPACE_DIM>* > temp_vertex_list;
    for (unsigned index = 0; index < vertexElements.size(); index++)
    {
        HEElement<SPACE_DIM>* p_temp_element = new HEElement<SPACE_DIM>(*vertexElements[index]);
        mElements.push_back(p_temp_element);

        HalfEdge<SPACE_DIM>* edge = p_temp_element->GetHalfEdge();
        HalfEdge<SPACE_DIM>* next_edge = edge->GetNextHalfEdge();
        while (edge!=next_edge)
        {
            temp_vertex_list.insert(next_edge->GetTargetVertex());
            next_edge = next_edge->GetNextHalfEdge();
        }
        temp_vertex_list.insert(edge->GetTargetVertex());
    }

    for (auto vertex:temp_vertex_list)
    {
        this->mNodes.push_back(vertex);
    }

}

template<unsigned int SPACE_DIM>
void HEVertexMesh<SPACE_DIM>::Clear()
{
    // Delete elements
    for (unsigned i = 0; i < mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete edges
    for (unsigned int i=0; i< mEdges.size(); ++i)
    {
        delete mEdges[i];
    }
    mEdges.clear();

    // Delete nodes
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();
}

template <unsigned int SPACE_DIM>
unsigned int HEVertexMesh<SPACE_DIM>::SolveNodeMapping(unsigned int index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template<unsigned int SPACE_DIM>
unsigned int HEVertexMesh<SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template<unsigned int SPACE_DIM>
unsigned int HEVertexMesh<SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>* HEVertexMesh<SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index<mElements.size());
    return mElements[index];
}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> HEVertexMesh<SPACE_DIM>::GetCentroidOfElement(unsigned index)
{
    HEElement<SPACE_DIM>* p_element = GetElement(index);
    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

    switch (SPACE_DIM)
    {
        case 1:
        {
            centroid = 0.5 * (p_element->GetNodeLocation(0) + p_element->GetNodeLocation(1));
        }
        break;
        case 2:
        {
            double centroid_x = 0;
            double centroid_y = 0;

            // Note that we cannot use GetVolumeOfElement() below as it returns the absolute, rather than signed, area
            double element_signed_area = 0.0;

            HalfEdge<SPACE_DIM>* edge = p_element->GetHalfEdge();
            // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
            c_vector<double, SPACE_DIM> first_node_location = edge->GetTargetVertex()->rGetLocation();
            c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

            // Loop over vertices
            HalfEdge<SPACE_DIM>* next_edge = edge->GetNextHalfEdge();
            do
            {
                c_vector<double, SPACE_DIM> next_node_location = next_edge->GetTargetVertex()->rGetLocation();
                c_vector<double, SPACE_DIM> pos_2 = this->GetVectorFromAtoB(first_node_location, next_node_location);

                double this_x = pos_1[0];
                double this_y = pos_1[1];
                double next_x = pos_2[0];
                double next_y = pos_2[1];

                double signed_area_term = this_x * next_y - this_y * next_x;

                centroid_x += (this_x + next_x) * signed_area_term;
                centroid_y += (this_y + next_y) * signed_area_term;
                element_signed_area += 0.5 * signed_area_term;

                pos_1 = pos_2;
                next_edge = next_edge->GetNextHalfEdge();
            }while(next_edge!=edge->GetNextHalfEdge());

            assert(element_signed_area != 0.0);

            // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
            centroid = first_node_location;
            centroid(0) += centroid_x / (6.0 * element_signed_area);
            centroid(1) += centroid_y / (6.0 * element_signed_area);
        }
        break;
        case 3:
        {
            EXCEPTION("Centroid computation for 3D mesh is not supported yet");
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}

template class Edge<1>;
template class Edge<2>;
template class Edge<3>;
template class HEVertexMesh<1>;
template class HEVertexMesh<2>;
template class HEVertexMesh<3>;
