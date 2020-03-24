#ifndef TESTDCELMESH_HPP_
#define TESTDCELMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "HEElement.hpp"
#include "HalfEdge.hpp"
#include "HENode.hpp"
#include "VertexMesh.hpp"
#include "HEVertexMesh.hpp"
#include "FakePetscSetup.hpp"
class TestHEVertexMesh : public CxxTest::TestSuite
{
public:
    void TestBasicMesh()
    {
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));

        vertices.push_back(new HENode<2>(4,false,2.0,0.0));
        vertices.push_back(new HENode<2>(5,false,2.0,1.0));

        std::vector<HENode<2>* > element_1_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        std::vector<HENode<2>* > element_2_vertices = {vertices[1], vertices[4], vertices[5], vertices[2]};
        HEElement<2> element_1(0,element_1_vertices);
        HEElement<2> element_2(1,element_2_vertices);
        std::vector<HEElement<2>* > elements = {&element_1, &element_2};
        HEVertexMesh<2> mesh(vertices, elements);

        HalfEdge<2>* edge = element_1.GetHalfEdge();
        HalfEdge<2>* walker = edge->GetNextHalfEdge();
        //Find common halfedge
        HalfEdge<2>* common_edge = nullptr;
        bool edge_found = false;
        while(!edge_found)
        {
            HalfEdge<2>* twin_walker_edge = walker->GetTwinHalfEdge();
            HalfEdge<2>* neigh_edge = element_2.GetHalfEdge();
            HalfEdge<2>* next_neigh_edge = neigh_edge->GetNextHalfEdge();
            while(next_neigh_edge!=neigh_edge)
            {
                if (twin_walker_edge==next_neigh_edge)
                {
                    edge_found = true;
                    common_edge = walker;
                }
                next_neigh_edge = next_neigh_edge->GetNextHalfEdge();
            }
            walker = walker->GetNextHalfEdge();
        }

        TS_ASSERT(edge_found);
        TS_ASSERT(common_edge);
        c_vector<double, 2> common_edge_node = common_edge->GetTargetNode()->rGetLocation();
        c_vector<double, 2> common_edge_twin_node = common_edge->GetTwinHalfEdge()->GetTargetNode()->rGetLocation();
        c_vector<double, 2> common_edge_prev_edge_node
        = common_edge->GetPreviousHalfEdge()->GetTargetNode()->rGetLocation();
        c_vector<double, 2> common_edge_twin_prev_edge_node
        = common_edge->GetTwinHalfEdge()->GetPreviousHalfEdge()->GetTargetNode()->rGetLocation();
        TS_ASSERT_DELTA(common_edge_node(0), common_edge_twin_prev_edge_node(0),1e-12);
        TS_ASSERT_DELTA(common_edge_node(1), common_edge_twin_prev_edge_node(1),1e-12);
        TS_ASSERT_DELTA(common_edge_prev_edge_node(0), common_edge_twin_node(0),1e-12);
        TS_ASSERT_DELTA(common_edge_prev_edge_node(1), common_edge_twin_node(1),1e-12);

        c_vector<double, 2> element_1_centroid = mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(element_1_centroid(0),0.5,1e-12);
        TS_ASSERT_DELTA(element_1_centroid(1),0.5,1e-12);
    }

    void TestBasic2DMesh()
    {
        // Create a 2D mesh comprising six nodes and two elements
        std::vector<HENode<2>*> nodes_2d;
        nodes_2d.push_back(new HENode<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new HENode<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new HENode<2>(2, false, 1.5, 1.0));
        nodes_2d.push_back(new HENode<2>(3, false, 1.0, 2.0));
        nodes_2d.push_back(new HENode<2>(4, false, 0.0, 1.0));
        nodes_2d.push_back(new HENode<2>(5, false, 2.0, 0.0));

        std::vector<std::vector<HENode<2>*> > nodes_elements_2d(2);
        for (unsigned i=0; i<5; i++)
        {
            nodes_elements_2d[0].push_back(nodes_2d[i]);
        }
        nodes_elements_2d[1].push_back(nodes_2d[1]);
        nodes_elements_2d[1].push_back(nodes_2d[5]);
        nodes_elements_2d[1].push_back(nodes_2d[2]);

        std::vector<HEElement<2>*> elements_2d;
        elements_2d.push_back(new HEElement<2>(0, nodes_elements_2d[0]));
        elements_2d.push_back(new HEElement<2>(1, nodes_elements_2d[1]));

        HEVertexMesh<2> mesh_2d(nodes_2d, elements_2d);

        // Test the mesh has the correct numbers of nodes and elements
        TS_ASSERT_EQUALS(mesh_2d.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_2d.GetNumNodes(), 6u);

        // Further tests of the mesh
        TS_ASSERT_DELTA(mesh_2d.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
        TS_ASSERT_EQUALS(mesh_2d.GetElement(1)->GetHalfEdge()->GetNextHalfEdge()->GetTargetNode()->GetIndex(), 2u);

        // Test that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);

        // Nodes 0 and 4 are only in element 0
        nodes_2d[0]->UpdateElementIndices();
        nodes_2d[4]->UpdateElementIndices();
        TS_ASSERT_EQUALS(nodes_2d[0]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes_2d[4]->rGetContainingElementIndices(), temp_list1);

        // Nodes 1 and 2 are in element 0 and 1
        temp_list1.insert(1u);
        nodes_2d[1]->UpdateElementIndices();
        nodes_2d[2]->UpdateElementIndices();
        TS_ASSERT_EQUALS(nodes_2d[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes_2d[2]->rGetContainingElementIndices(), temp_list1);

        // Node 5 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        nodes_2d[5]->UpdateElementIndices();
        TS_ASSERT_EQUALS(nodes_2d[5]->rGetContainingElementIndices(), temp_list2);

        //Test full edges
        TS_ASSERT_EQUALS(mesh_2d.GetNumFullEdges(), 7u);
        //Traverse edges
        for(unsigned int i=0; i<mesh_2d.GetNumFullEdges(); ++i)
        {
            FullEdge<2>* full_edge = mesh_2d.GetFullEdge(i);
            //Check to see if full edges contain all the vertices
            TS_ASSERT_EQUALS(std::count(nodes_2d.begin(),nodes_2d.end(), full_edge->first->GetTargetNode()), 1u);
            TS_ASSERT_EQUALS(std::count(nodes_2d.begin(),nodes_2d.end(), full_edge->second->GetTargetNode()), 1u);
        }
    }

    void TestGetCentroidOfElement()
    {
        /**
         * Tests here are simply adapted from of VertexMesh tests
         * Also test for conversion of VertexMesh into HEVertexMesh
         */
        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        HEVertexMesh<2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 2> triangle_centroid = triangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(triangle_centroid[0], 2.0/3.0, 1e-4);
        TS_ASSERT_DELTA(triangle_centroid[1], 1.0/3.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        HEVertexMesh<2> square_mesh(square_nodes, square_elements);

        c_vector<double, 2> square_centroid = square_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(square_centroid[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(square_centroid[1], 0.5, 1e-6);

        // Test method with a single rectangular element away from the origin
        std::vector<Node<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new Node<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new Node<2>(3, false, 10.0, 14.0));
        std::vector<VertexElement<2,2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new VertexElement<2,2>(0, far_rectangle_nodes));
        HEVertexMesh<2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 2> far_rectangle_centroid = far_rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_centroid[0], 10.5, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_centroid[1], 12.0, 1e-4);

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<Node<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new Node<2>(0, false,  2.0*0.5*sqrt(3.0) - 1.0*0.5,  2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(1, false, -2.0*0.5*sqrt(3.0) - 1.0*0.5, -2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(2, false, -2.0*0.5*sqrt(3.0) + 1.0*0.5, -2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(3, false,  2.0*0.5*sqrt(3.0) + 1.0*0.5,  2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        std::vector<VertexElement<2,2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new VertexElement<2,2>(0, angled_rectangle_nodes));
        HEVertexMesh<2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 2> angled_rectangle_centroid = angled_rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(angled_rectangle_centroid(1), 0.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> circle_elements;
        circle_elements.push_back(new VertexElement<2,2>(0, circle_nodes));
        HEVertexMesh<2> circle_mesh(circle_nodes, circle_elements);

        c_vector<double, 2> circle_centroid = circle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(circle_centroid[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(circle_centroid[1], 0.0, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> hexagon_elements;
        hexagon_elements.push_back(new VertexElement<2,2>(0, hexagon_nodes));
        HEVertexMesh<2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        c_vector<double, 2> hexagon_centroid = hexagon_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(hexagon_centroid[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(hexagon_centroid[1], 0.0, 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        HEVertexMesh<2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 2> irregular_centroid = irregular_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(irregular_centroid[0], 2.6269, 1e-4);
        TS_ASSERT_DELTA(irregular_centroid[1], 0.8930, 1e-4);

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);
        HEVertexMesh<2> regular_he_mesh(&regular_mesh);

        c_vector<double, 2> regular_centroid_5 = regular_he_mesh.GetCentroidOfElement(5);
        TS_ASSERT_DELTA(regular_centroid_5[0], 2.0,           1e-4);
        TS_ASSERT_DELTA(regular_centroid_5[1], 2.5/sqrt(3.0), 1e-4);

        c_vector<double, 2> regular_centroid_7 = regular_he_mesh.GetCentroidOfElement(7);
        TS_ASSERT_DELTA(regular_centroid_7[0], 4.0,           1e-4);
        TS_ASSERT_DELTA(regular_centroid_7[1], 2.5/sqrt(3.0), 1e-4);
    }

    void TestGetVolumeOfElement()
    {
        /**
         * Tests here are simply adapted from of VertexMesh tests
         */
        // Note that the method GetVolumeOfElement() only works in 2D and 3D

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> triangle_elements;
        triangle_elements.push_back(new VertexElement<2,2>(0, triangle_nodes));
        HEVertexMesh<2> triangle_mesh(triangle_nodes, triangle_elements);

        TS_ASSERT_DELTA(triangle_mesh.GetVolumeOfElement(0), 1.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));
        HEVertexMesh<2> square_mesh(square_nodes, square_elements);

        TS_ASSERT_DELTA(square_mesh.GetVolumeOfElement(0), 1.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> circle_elements;
        circle_elements.push_back(new VertexElement<2,2>(0, circle_nodes));
        HEVertexMesh<2> circle_mesh(circle_nodes, circle_elements);

        TS_ASSERT_DELTA(circle_mesh.GetVolumeOfElement(0), M_PI, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<VertexElement<2,2>*> hexagon_elements;
        hexagon_elements.push_back(new VertexElement<2,2>(0, hexagon_nodes));
        HEVertexMesh<2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        TS_ASSERT_DELTA(hexagon_mesh.GetVolumeOfElement(0), 1.5*sqrt(3.0), 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158,  1.5588));
        std::vector<VertexElement<2,2>*> irregular_elements;
        irregular_elements.push_back(new VertexElement<2,2>(0, irregular_nodes));
        HEVertexMesh<2> irregular_mesh(irregular_nodes, irregular_elements);

        TS_ASSERT_DELTA(irregular_mesh.GetVolumeOfElement(0), 1.4684, 1e-3);

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);
        HEVertexMesh<2> regular_he_mesh(&regular_mesh);

        for (HEVertexMesh<2>::HEElementIterator iter = regular_he_mesh.GetElementIteratorBegin();
                iter != regular_he_mesh.GetElementIteratorEnd();
                ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TS_ASSERT_DELTA(regular_he_mesh.GetVolumeOfElement(elem_index), 0.5*sqrt(3.0), 1e-4);
        }
    }
};
#endif /*TESTDCELMESH_HPP_*/


