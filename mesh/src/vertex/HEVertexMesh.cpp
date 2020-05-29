/*
 * HEVertexMesh.cpp
 *
 *  Created on: 11 Mar 2020
 *      Author: aydar
 */

#include "HEVertexMesh.hpp"
#include <map>
#include <set>
#include "RandomNumberGenerator.hpp"

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEVertexMesh(std::vector<HENode<SPACE_DIM>* > vertices,
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
    ConstructFullEdges();
}

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                      std::vector<VertexElement<SPACE_DIM, SPACE_DIM>*> vertexElements)
{
    Clear();

    ConvertFromVertexMesh(nodes, vertexElements);
    ConstructFullEdges();
}

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEVertexMesh(VertexMesh<SPACE_DIM, SPACE_DIM>* vertex_mesh)
{
    Clear();

    const unsigned int n_vertex_elements = vertex_mesh->GetNumElements();
    std::vector<VertexElement<SPACE_DIM, SPACE_DIM>* > vertex_elements(n_vertex_elements);
    for (unsigned int i=0; i<n_vertex_elements; ++i)
    {
        vertex_elements[i] = vertex_mesh->GetElement(i);
    }

    const unsigned int n_nodes = vertex_mesh->GetNumNodes();
    std::vector<Node<SPACE_DIM>* > nodes(n_nodes);
    for (unsigned int i=0; i<n_nodes; ++i)
    {
        nodes[i] = vertex_mesh->GetNode(i);
    }

    ConvertFromVertexMesh(nodes, vertex_elements);
    ConstructFullEdges();
}

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEVertexMesh()
{
}

template<unsigned int SPACE_DIM>
HEVertexMesh<SPACE_DIM>::~HEVertexMesh()
{
}

template<unsigned int SPACE_DIM>
void HEVertexMesh<SPACE_DIM>::ConvertFromVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                    std::vector<VertexElement<SPACE_DIM, SPACE_DIM>*> vertexElements)
{
    Clear();

    std::set<HENode<SPACE_DIM>* > temp_node_list;
    for (unsigned index = 0; index < vertexElements.size(); index++)
    {
        //Create HEElement from VertexElement...
        HEElement<SPACE_DIM>* p_temp_element = new HEElement<SPACE_DIM>(*vertexElements[index]);
        mElements.push_back(p_temp_element);

        //Populate HEVertex list
        HalfEdge<SPACE_DIM>* edge = p_temp_element->GetHalfEdge();
        HalfEdge<SPACE_DIM>* next_edge = edge->GetNextHalfEdge();
        while (edge!=next_edge)
        {
            temp_node_list.insert(next_edge->GetTargetNode());
            next_edge = next_edge->GetNextHalfEdge();
        }
        temp_node_list.insert(edge->GetTargetNode());
    }

    for (auto node:temp_node_list)
    {
        this->mNodes.push_back(node);
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
    for (unsigned int i=0; i< mFullEdges.size(); ++i)
    {
        delete mFullEdges[i];
    }
    mFullEdges.clear();

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
unsigned int HEVertexMesh<SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned int SPACE_DIM>
HEElement<SPACE_DIM>* HEVertexMesh<SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index<mElements.size());
    return mElements[index];
}

template<unsigned int SPACE_DIM>
HENode<SPACE_DIM>* HEVertexMesh<SPACE_DIM>::GetNode(unsigned index) const
{
    assert(index<this->mNodes.size());
    return static_cast<HENode<SPACE_DIM>* >(this->mNodes[index]);
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

            // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
            HalfEdge<SPACE_DIM>* edge = p_element->GetHalfEdge();
            c_vector<double, SPACE_DIM> first_node_location;
            first_node_location = edge->GetOriginNode()->rGetLocation();
            c_vector<double, SPACE_DIM> pos_1;
            pos_1 = zero_vector<double>(SPACE_DIM);

            // Loop over vertices
            HalfEdge<SPACE_DIM>* next_edge = edge;
            do
            {
                c_vector<double, SPACE_DIM> next_node_location = next_edge->GetTargetNode()->rGetLocation();
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
            }while(next_edge!=edge);

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

template <unsigned SPACE_DIM>
double HEVertexMesh<SPACE_DIM>::GetVolumeOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    HEElement<SPACE_DIM>* p_element = GetElement(index);

    return p_element->GetVolume();
}

template<unsigned int SPACE_DIM>
double HEVertexMesh<SPACE_DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    HEElement<SPACE_DIM>* p_element = GetElement(index);

    return p_element->GetSurfaceArea();
}

template<unsigned int SPACE_DIM>
double HEVertexMesh<SPACE_DIM>::ComputeVolumeOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3);
    return GetElement(index)->ComputeVolume();
}

template<unsigned int SPACE_DIM>
double HEVertexMesh<SPACE_DIM>::ComputeSurfaceAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3);
    return GetElement(index)->ComputeSurfaceArea();
}

template <unsigned SPACE_DIM>
c_vector<double, 3> HEVertexMesh<SPACE_DIM>::CalculateMomentsOfElement(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    // Define helper variables
    HEElement<SPACE_DIM>* p_element = GetElement(index);
    c_vector<double, 3> moments = zero_vector<double>(3);

    // Since we compute I_xx, I_yy and I_xy about the centroid, we must shift each vertex accordingly
    c_vector<double, SPACE_DIM> centroid = GetCentroidOfElement(index);

    // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
    HalfEdge<SPACE_DIM>* edge = p_element->GetHalfEdge();
    c_vector<double, SPACE_DIM> first_node_location;
    first_node_location = edge->GetOriginNode()->rGetLocation();

    c_vector<double, SPACE_DIM> pos_1;
    pos_1 = this->GetVectorFromAtoB(centroid, first_node_location);

    // Loop over vertices
    HalfEdge<SPACE_DIM>* next_edge = edge;
    do
    {
        c_vector<double, SPACE_DIM> next_node_location = next_edge->GetTargetNode()->rGetLocation();
        c_vector<double, SPACE_DIM> pos_2 = this->GetVectorFromAtoB(centroid, next_node_location);

        double signed_area_term = pos_1(0) * pos_2(1) - pos_2(0) * pos_1(1);
        // Ixx
        moments(0) += (pos_1(1) * pos_1(1) + pos_1(1) * pos_2(1) + pos_2(1) * pos_2(1)) * signed_area_term;

        // Iyy
        moments(1) += (pos_1(0) * pos_1(0) + pos_1(0) * pos_2(0) + pos_2(0) * pos_2(0)) * signed_area_term;

        // Ixy
        moments(2) += (pos_1(0) * pos_2(1) + 2 * pos_1(0) * pos_1(1) + 2 * pos_2(0) * pos_2(1) + pos_2(0) * pos_1(1)) * signed_area_term;

        pos_1 = pos_2;
        next_edge = next_edge->GetNextHalfEdge();
    }while(next_edge!=edge);

    moments(0) /= 12;
    moments(1) /= 12;
    moments(2) /= 24;

    /*
     * If the nodes owned by the element were supplied in a clockwise rather
     * than anticlockwise manner, or if this arose as a result of enforcing
     * periodicity, then our computed quantities will be the wrong sign, so
     * we need to fix this.
     */
    if (moments(0) < 0.0)
    {
        moments(0) = -moments(0);
        moments(1) = -moments(1);
        moments(2) = -moments(2);
    }
    return moments;
}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> HEVertexMesh<SPACE_DIM>::GetShortAxisOfElement(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);

    // Calculate the moments of the element about its centroid (recall that I_xx and I_yy must be non-negative)
    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    // Normalise the moments vector to remove problem of a very small discriminant (see #2874)
    moments /= norm_2(moments);

    // If the principal moments are equal...
    double discriminant = (moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2);
    if (fabs(discriminant) < DBL_EPSILON)
    {
        // ...then every axis through the centroid is a principal axis, so return a random unit vector
        short_axis(0) = RandomNumberGenerator::Instance()->ranf();
        short_axis(1) = sqrt(1.0 - short_axis(0) * short_axis(0));
    }
    else
    {
        // If the product of inertia is zero, then the coordinate axes are the principal axes
        if (fabs(moments(2)) < DBL_EPSILON)
        {
            if (moments(0) < moments(1))
            {
                short_axis(0) = 0.0;
                short_axis(1) = 1.0;
            }
            else
            {
                short_axis(0) = 1.0;
                short_axis(1) = 0.0;
            }
        }
        else
        {
            // Otherwise we find the eigenvector of the inertia matrix corresponding to the largest eigenvalue
            double lambda = 0.5 * (moments(0) + moments(1) + sqrt(discriminant));

            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - lambda) / moments(2);

            // Normalise the short axis before returning it
            short_axis /= norm_2(short_axis);
        }
    }

    return short_axis;
}

template <unsigned SPACE_DIM>
double HEVertexMesh<SPACE_DIM>::GetElongationShapeFactorOfElement(unsigned index)
{
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    c_vector<double, 3> moments = CalculateMomentsOfElement(index);

    double discriminant = sqrt((moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2));

    // Note that as the matrix of second moments of area is symmetric, both its eigenvalues are real
    double largest_eigenvalue = (moments(0) + moments(1) + discriminant) * 0.5;
    double smallest_eigenvalue = (moments(0) + moments(1) - discriminant) * 0.5;

    double elongation_shape_factor = sqrt(largest_eigenvalue / smallest_eigenvalue);
    return elongation_shape_factor;
}

template<unsigned int SPACE_DIM>
void HEVertexMesh<SPACE_DIM>::ConstructFullEdges()
{
    //The set below keeps track of halfedges that we already visited
    std::set<HalfEdge<SPACE_DIM>*> traversal_set;
    for (typename HEVertexMesh<SPACE_DIM>::HEElementIterator elem_iter = this->GetElementIteratorBegin();
             elem_iter != this->GetElementIteratorEnd();
             ++elem_iter)
    {
        HalfEdge<SPACE_DIM>* next_edge = elem_iter->GetHalfEdge();
        do
        {
            //If the halfedge has not been bundled into a full edge...
            if (traversal_set.count(next_edge)==0)
            {
                //... construct full edge of this halfedge and its twin
                FullEdge<SPACE_DIM>* edge = new FullEdge<SPACE_DIM>(next_edge);
                mFullEdges.push_back(edge);
                traversal_set.insert(next_edge);
                traversal_set.insert(next_edge->GetTwinHalfEdge());
                mHalfToFullEdgeMap.insert(std::pair<HalfEdge<SPACE_DIM>*, FullEdge<SPACE_DIM>* >(next_edge, edge));
                mHalfToFullEdgeMap.insert(std::pair<HalfEdge<SPACE_DIM>*, FullEdge<SPACE_DIM>* >(next_edge->GetTwinHalfEdge(), edge));
            }
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge!=elem_iter->GetHalfEdge());
    }
}

template<unsigned int SPACE_DIM>
FullEdge<SPACE_DIM>* HEVertexMesh<SPACE_DIM>::GetFullEdge(const unsigned int index) const
{
    assert(index<mFullEdges.size());
    return mFullEdges[index];
}

template<unsigned int SPACE_DIM>
unsigned int HEVertexMesh<SPACE_DIM>::GetNumFullEdges() const
{
    return mFullEdges.size();
}

template<unsigned int SPACE_DIM>
unsigned int HEVertexMesh<SPACE_DIM>::GetNumAllFullEdges() const
{
    return mFullEdges.size();
}

template class HEVertexMesh<1>;
template class HEVertexMesh<2>;
template class HEVertexMesh<3>;
