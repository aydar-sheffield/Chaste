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

    void TestNodeIterator()
    {
        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        VertexMesh<2,2> v_mesh;
        v_mesh.ConstructFromMeshReader(mesh_reader);

        HEVertexMesh<2> mesh;
        mesh.ConvertFromVertexMesh(&v_mesh);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 30u);

        unsigned counter = 0;
        for (HEVertexMesh<2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
                iter != mesh.GetNodeIteratorEnd();
                ++iter)
        {
            unsigned node_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, node_index); // assumes the iterator will give nodes 0,1..,N in that order
            counter++;
        }
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), counter);

        // Check that the node iterator correctly handles deleted nodes
        mesh.GetNode(0)->MarkAsDeleted();

        counter = 0;
        for (HEVertexMesh<2>::NodeIterator iter = mesh.GetNodeIteratorBegin();
                iter != mesh.GetNodeIteratorEnd();
                ++iter)
        {
            unsigned node_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter+1, node_index); // assumes the iterator will give nodes 1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), counter+1);

        // For coverage, test with an empty mesh
        HEVertexMesh<2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        HEVertexMesh<2>::NodeIterator iter = empty_mesh.GetNodeIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (iter != empty_mesh.GetNodeIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);

        // Coverage of AbstractMesh::SetElementOwnerships()
        TS_ASSERT_THROWS_NOTHING(empty_mesh.SetElementOwnerships());
    }

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
        TS_ASSERT_EQUALS(nodes_2d[0]->GetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes_2d[4]->GetContainingElementIndices(), temp_list1);

        // Nodes 1 and 2 are in element 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(nodes_2d[1]->GetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes_2d[2]->GetContainingElementIndices(), temp_list1);

        // Node 5 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(nodes_2d[5]->GetContainingElementIndices(), temp_list2);

        //Test full edges
        TS_ASSERT_EQUALS(mesh_2d.GetNumFullEdges(), 7u);
        //Traverse edges
        for(unsigned int i=0; i<mesh_2d.GetNumFullEdges(); ++i)
        {
            FullEdge<2>* full_edge = mesh_2d.GetFullEdge(i);
            //Check to see if full edges contain all the vertices
            TS_ASSERT_EQUALS(std::count(nodes_2d.begin(),nodes_2d.end(), (*full_edge)(0)->GetTargetNode()), 1u);
            TS_ASSERT_EQUALS(std::count(nodes_2d.begin(),nodes_2d.end(), (*full_edge)(1)->GetTargetNode()), 1u);
        }
    }

    void TestGetRosetteRankOfElement()
    {
        // Test method on a honeycomb mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> v_regular_mesh;
        v_regular_mesh.ConstructFromMeshReader(mesh_reader);

        HEVertexMesh<2> regular_mesh;
        regular_mesh.ConvertFromVertexMesh(&v_regular_mesh);

        // The rosette rank of each element should be 3
        for (unsigned elem_idx = 0 ; elem_idx < regular_mesh.GetNumElements() ; elem_idx++)
        {
            TS_ASSERT_EQUALS(regular_mesh.GetRosetteRankOfElement(elem_idx), 3u);
        }

        // Generate a five-element rosette
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.00000, 0.00000));
        nodes.push_back(new HENode<2>(2, true, 0.80902, 0.58779));
        nodes.push_back(new HENode<2>(3, true, 0.30902, 0.95106));
        nodes.push_back(new HENode<2>(4, true, -0.30902, 0.95106));
        nodes.push_back(new HENode<2>(5, true, -0.80902, 0.58779));
        nodes.push_back(new HENode<2>(6, true, -1.00000, 0.00000));
        nodes.push_back(new HENode<2>(7, true, -0.80902, -0.58779));
        nodes.push_back(new HENode<2>(8, true, -0.30902, -0.95106));
        nodes.push_back(new HENode<2>(9, true, 0.30902, -0.95106));
        nodes.push_back(new HENode<2>(10, true, 0.80902, -0.58779));

        // Make 5 quadrangular cells
        std::vector<HENode<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        std::vector<HENode<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        std::vector<HENode<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[5]);
        nodes_elem_3.push_back(nodes[6]);
        nodes_elem_3.push_back(nodes[7]);
        std::vector<HENode<2>*> nodes_elem_4;
        nodes_elem_4.push_back(nodes[0]);
        nodes_elem_4.push_back(nodes[7]);
        nodes_elem_4.push_back(nodes[8]);
        nodes_elem_4.push_back(nodes[9]);
        std::vector<HENode<2>*> nodes_elem_5;
        nodes_elem_5.push_back(nodes[0]);
        nodes_elem_5.push_back(nodes[9]);
        nodes_elem_5.push_back(nodes[10]);
        nodes_elem_5.push_back(nodes[1]);

        // Make 5 vertex elements
        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_3));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_4));
        vertex_elements.push_back(new HEElement<2>(4, nodes_elem_5));

        HEVertexMesh<2> rosette_mesh(nodes, vertex_elements);

        // The rosette rank of each element should be five
        for (unsigned elem_idx = 0 ; elem_idx < rosette_mesh.GetNumElements() ; elem_idx++)
        {
            TS_ASSERT_EQUALS(rosette_mesh.GetRosetteRankOfElement(elem_idx), 5u);
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

        triangle_mesh.GetElement(0)->ComputeVolume();
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

        square_mesh.GetElement(0)->ComputeVolume();
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

        circle_mesh.GetElement(0)->ComputeVolume();
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

        hexagon_mesh.GetElement(0)->ComputeVolume();
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

        irregular_mesh.GetElement(0)->ComputeVolume();
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
            regular_he_mesh.GetElement(elem_index)->ComputeVolume();
            TS_ASSERT_DELTA(regular_he_mesh.GetVolumeOfElement(elem_index), 0.5*sqrt(3.0), 1e-4);
        }
    }

    void TestGetSurfaceAreaOfElement()
    {
        // Note that the method GetSurfaceAreaOfElement() only works in 2D and 3D

        // Test method with a single triangular element
        std::vector<HENode<2>*> triangle_nodes;
        triangle_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new HENode<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new HENode<2>(2, false, 0.0, 1.0));
        std::vector<HEElement<2>*> triangle_elements;
        triangle_elements.push_back(new HEElement<2>(0, triangle_nodes));
        HEVertexMesh<2> triangle_mesh(triangle_nodes, triangle_elements);

        TS_ASSERT_DELTA(triangle_mesh.GetSurfaceAreaOfElement(0), 3.0 + sqrt(5.0), 1e-4);

        // Test method with a single square element
        std::vector<HENode<2>*> square_nodes;
        square_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new HENode<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new HENode<2>(3, false, 0.0, 1.0));
        std::vector<HEElement<2>*> square_elements;
        square_elements.push_back(new HEElement<2>(0, square_nodes));
        HEVertexMesh<2> square_mesh(square_nodes, square_elements);

        TS_ASSERT_DELTA(square_mesh.GetSurfaceAreaOfElement(0), 4.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<HENode<2>*> circle_nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            circle_nodes.push_back(new HENode<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<HEElement<2>*> circle_elements;
        circle_elements.push_back(new HEElement<2>(0, circle_nodes));
        HEVertexMesh<2> circle_mesh(circle_nodes, circle_elements);

        TS_ASSERT_DELTA(circle_mesh.GetSurfaceAreaOfElement(0), 2.0*M_PI, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<HENode<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new HENode<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<HEElement<2>*> hexagon_elements;
        hexagon_elements.push_back(new HEElement<2>(0, hexagon_nodes));
        HEVertexMesh<2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        TS_ASSERT_DELTA(hexagon_mesh.GetSurfaceAreaOfElement(0), 6.0, 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<HENode<2>*> irregular_nodes;
        irregular_nodes.push_back(new HENode<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new HENode<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new HENode<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new HENode<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new HENode<2>(4, false, 3.3158,  1.5588));
        std::vector<HEElement<2>*> irregular_elements;
        irregular_elements.push_back(new HEElement<2>(0, irregular_nodes));
        HEVertexMesh<2> irregular_mesh(irregular_nodes, irregular_elements);

        TS_ASSERT_DELTA(irregular_mesh.GetSurfaceAreaOfElement(0), 5.1263, 1e-3);

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> v_regular_mesh;
        v_regular_mesh.ConstructFromMeshReader(mesh_reader);

        HEVertexMesh<2> regular_mesh;
        regular_mesh.ConvertFromVertexMesh(&v_regular_mesh);

        for (HEVertexMesh<2>::HEElementIterator iter = regular_mesh.GetElementIteratorBegin();
                iter != regular_mesh.GetElementIteratorEnd();
                ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TS_ASSERT_DELTA(regular_mesh.GetSurfaceAreaOfElement(elem_index), 2*sqrt(3.0), 1e-4);
        }
    }

    void TestGetAreaGradientOfElementAtNode()
    {
        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new HENode<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new HENode<2>(3, false, 0.0, 1.0));

        std::vector<VertexElement<2,2>*> square_elements;
        square_elements.push_back(new VertexElement<2,2>(0, square_nodes));

        VertexMesh<2,2> square_mesh(square_nodes, square_elements);

        VertexElement<2,2>* p_element = square_mesh.GetElement(0);

        c_vector<double, 2> element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 0);
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);

        element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 1);
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);

        element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 2);
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);

        element_area_gradient = square_mesh.GetAreaGradientOfElementAtNode(p_element, 3);
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
    }

    void TestGetPerimeterGradientAtNode()
    {
        // Test method with a single square element
        std::vector<HENode<2>*> square_nodes;
        square_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new HENode<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new HENode<2>(3, false, 0.0, 1.0));

        std::vector<HEElement<2>*> square_elements;
        square_elements.push_back(new HEElement<2>(0, square_nodes));

        HEVertexMesh<2> square_mesh(square_nodes, square_elements);

        HEElement<2>* p_element = square_mesh.GetElement(0);

        c_vector<double, 2> element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 0);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 1);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 2);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);

        element_perimeter_gradient = square_mesh.GetPerimeterGradientOfElementAtNode(p_element, 3);
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);
    }



    void TestCalculateMomentsOfElement()
    {
        // Test method with a single triangular element
        std::vector<HENode<2>*> isos_triangle_nodes;
        isos_triangle_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        isos_triangle_nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        isos_triangle_nodes.push_back(new HENode<2>(2, false, 0.0, 1.0));
        std::vector<HEElement<2>*> isos_triangle_elements;
        isos_triangle_elements.push_back(new HEElement<2>(0, isos_triangle_nodes));
        HEVertexMesh<2> isos_triangle_mesh(isos_triangle_nodes, isos_triangle_elements);

        c_vector<double, 3> isos_triangle_moments = isos_triangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(isos_triangle_moments[0], 0.0277, 1e-4);
        TS_ASSERT_DELTA(isos_triangle_moments[1], 0.0277, 1e-4);
        TS_ASSERT_DELTA(isos_triangle_moments[2], -0.0138, 1e-4);

        // Test method with a single triangular element
        std::vector<HENode<2>*> triangle_nodes;
        triangle_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new HENode<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new HENode<2>(2, false, 0.0, 1.0));
        std::vector<HEElement<2>*> triangle_elements;
        triangle_elements.push_back(new HEElement<2>(0, triangle_nodes));
        HEVertexMesh<2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 3> triangle_moments = triangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(triangle_moments(0), 1.0/18.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(triangle_moments(1), 2.0/9.0, 1e-6);   // Iyy
        TS_ASSERT_DELTA(triangle_moments(2), -5.0/90.0, 1e-6); // Ixy

        // Test method with a single rectangular element parallel to the x-axis
        std::vector<HENode<2>*> horizontal_rectangle_nodes;
        horizontal_rectangle_nodes.push_back(new HENode<2>(0, false,  2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new HENode<2>(1, false, -2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new HENode<2>(2, false, -2.0, -1.0));
        horizontal_rectangle_nodes.push_back(new HENode<2>(3, false,  2.0, -1.0));
        std::vector<HEElement<2>*> horizontal_rectangle_elements;
        horizontal_rectangle_elements.push_back(new HEElement<2>(0, horizontal_rectangle_nodes));
        HEVertexMesh<2> horizontal_rectangle_mesh(horizontal_rectangle_nodes, horizontal_rectangle_elements);

        c_vector<double, 3> horizontal_rectangle_moments = horizontal_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(horizontal_rectangle_moments(0), 8.0/3.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(horizontal_rectangle_moments(1), 32.0/3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(horizontal_rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        // Test method with the same shape, but supply nodes in clockwise manner
        std::vector<HENode<2>*> clockwise_rectangle_nodes;
        clockwise_rectangle_nodes.push_back(new HENode<2>(0, false, -2.0, -1.0));
        clockwise_rectangle_nodes.push_back(new HENode<2>(1, false, -2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new HENode<2>(2, false,  2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new HENode<2>(3, false,  2.0, -1.0));
        std::vector<HEElement<2>*> clockwise_rectangle_elements;
        clockwise_rectangle_elements.push_back(new HEElement<2>(0, clockwise_rectangle_nodes));
        HEVertexMesh<2> clockwise_rectangle_mesh(clockwise_rectangle_nodes, clockwise_rectangle_elements);

        c_vector<double, 3> clockwise_rectangle_moments = clockwise_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(clockwise_rectangle_moments(0), 8.0/3.0, 1e-6);  // Ixx
        TS_ASSERT_DELTA(clockwise_rectangle_moments(1), 32.0/3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(clockwise_rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        // Test method with a single rectangular element parallel to the y-axis
        std::vector<HENode<2>*> vertical_rectangle_nodes;
        vertical_rectangle_nodes.push_back(new HENode<2>(0, false,  1.0,  2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(1, false, -1.0,  2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(2, false, -1.0, -2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(3, false,  1.0, -2.0));
        std::vector<HEElement<2>*> vertical_rectangle_elements;
        vertical_rectangle_elements.push_back(new HEElement<2>(0, vertical_rectangle_nodes));
        HEVertexMesh<2> vertical_rectangle_mesh(vertical_rectangle_nodes, vertical_rectangle_elements);

        c_vector<double, 3> vertical_rectangle_moments = vertical_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_moments(0), 32.0/3.0, 1e-6); // Ixx
        TS_ASSERT_DELTA(vertical_rectangle_moments(1), 8.0/3.0, 1e-6);  // Iyy
        TS_ASSERT_DELTA(vertical_rectangle_moments(2), 0.0, 1e-6);      // Ixy = 0 by symmetry

        // Test method with a single rectangular element away from the origin
        std::vector<HENode<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new HENode<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new HENode<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new HENode<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new HENode<2>(3, false, 10.0, 14.0));
        std::vector<HEElement<2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new HEElement<2>(0, far_rectangle_nodes));
        HEVertexMesh<2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 3> far_rectangle_moments = far_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_moments[0], 16.0/3.0, 1e-4); // Ixx
        TS_ASSERT_DELTA(far_rectangle_moments[1], 1.0/3.0, 1e-4);  // Iyy
        TS_ASSERT_DELTA(far_rectangle_moments[2], 0.0, 1e-4);      // Ixy = 0 by symmetry

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<HENode<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new HENode<2>(0, false,  2.0*0.5*sqrt(3.0) - 1.0*0.5,  2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new HENode<2>(1, false, -2.0*0.5*sqrt(3.0) - 1.0*0.5, -2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new HENode<2>(2, false, -2.0*0.5*sqrt(3.0) + 1.0*0.5, -2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new HENode<2>(3, false,  2.0*0.5*sqrt(3.0) + 1.0*0.5,  2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        std::vector<HEElement<2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new HEElement<2>(0, angled_rectangle_nodes));
        HEVertexMesh<2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 3> angled_rectangle_moments = angled_rectangle_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_moments[0], 14.0/3.0, 1e-4);    // Ixx
        TS_ASSERT_DELTA(angled_rectangle_moments[1], 26.0/3.0, 1e-4);    // Iyy
        TS_ASSERT_DELTA(angled_rectangle_moments[2], 2.0*sqrt(3.0), 1e-4); // Ixy

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<HENode<2>*> irregular_nodes;
        irregular_nodes.push_back(new HENode<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new HENode<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new HENode<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new HENode<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new HENode<2>(4, false, 3.3158,  1.5588));
        std::vector<HEElement<2>*> irregular_elements;
        irregular_elements.push_back(new HEElement<2>(0, irregular_nodes));
        HEVertexMesh<2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 3> irregular_moments = irregular_mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(irregular_moments[0], 0.3521, 1e-4); // Ixx
        TS_ASSERT_DELTA(irregular_moments[1], 0.1264, 1e-4); // Iyy
        TS_ASSERT_DELTA(irregular_moments[2], 0.1162, 1e-4); // Ixy

        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        VertexMesh<2,2> v_regular_mesh;
        v_regular_mesh.ConstructFromMeshReader(mesh_reader);

        HEVertexMesh<2> regular_mesh;
        regular_mesh.ConvertFromVertexMesh(&v_regular_mesh);

        for (unsigned i=0; i<regular_mesh.GetNumElements(); i++)
        {
            c_vector<double, 3> regular_moments = regular_mesh.CalculateMomentsOfElement(i);
            TS_ASSERT_DELTA(regular_moments(0), 5*sqrt(3.0)/16/9, 1e-6); // Ixx
            TS_ASSERT_DELTA(regular_moments(1), 5*sqrt(3.0)/16/9, 1e-6); // Iyy
            TS_ASSERT_DELTA(regular_moments(2), 0.0,            1e-6); // Ixy = 0 by symmetry
        }
    }

    void TestGetShortAxisOfElement()
    {
        // Test method with a single triangular element
        std::vector<HENode<2>*> triangle_nodes;
        triangle_nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        triangle_nodes.push_back(new HENode<2>(2, false, 0.0, 1.0));
        std::vector<HEElement<2>*> triangle_elements;
        triangle_elements.push_back(new HEElement<2>(0, triangle_nodes));
        HEVertexMesh<2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 2> triangle_short_axis = triangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(triangle_short_axis[0], 1/sqrt(2.0), 1e-4);
        TS_ASSERT_DELTA(triangle_short_axis[1], 1/sqrt(2.0), 1e-4);

        // Test method with a single rectangular element parallel to the x-axis
        std::vector<HENode<2>*> horizontal_rectangle_nodes;
        horizontal_rectangle_nodes.push_back(new HENode<2>(0, false,  2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new HENode<2>(1, false, -2.0,  1.0));
        horizontal_rectangle_nodes.push_back(new HENode<2>(2, false, -2.0, -1.0));
        horizontal_rectangle_nodes.push_back(new HENode<2>(3, false,  2.0, -1.0));
        std::vector<HEElement<2>*> horizontal_rectangle_elements;
        horizontal_rectangle_elements.push_back(new HEElement<2>(0, horizontal_rectangle_nodes));
        HEVertexMesh<2> horizontal_rectangle_mesh(horizontal_rectangle_nodes, horizontal_rectangle_elements);

        c_vector<double, 2> horizontal_rectangle_short_axis = horizontal_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(horizontal_rectangle_short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(horizontal_rectangle_short_axis(1), 1.0, 1e-6);

        // Test method with the same shape, but supply nodes in clockwise manner
        std::vector<HENode<2>*> clockwise_rectangle_nodes;
        clockwise_rectangle_nodes.push_back(new HENode<2>(0, false, -2.0, -1.0));
        clockwise_rectangle_nodes.push_back(new HENode<2>(1, false, -2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new HENode<2>(2, false,  2.0,  1.0));
        clockwise_rectangle_nodes.push_back(new HENode<2>(3, false,  2.0, -1.0));
        std::vector<HEElement<2>*> clockwise_rectangle_elements;
        clockwise_rectangle_elements.push_back(new HEElement<2>(0, clockwise_rectangle_nodes));
        HEVertexMesh<2> clockwise_rectangle_mesh(clockwise_rectangle_nodes, clockwise_rectangle_elements);

        c_vector<double, 3> clockwise_rectangle_short_axis = clockwise_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(clockwise_rectangle_short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(clockwise_rectangle_short_axis(1), 1.0, 1e-6);

        // Test method with a single rectangular element parallel to the y-axis
        std::vector<HENode<2>*> vertical_rectangle_nodes;
        vertical_rectangle_nodes.push_back(new HENode<2>(0, false,  1.0,  2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(1, false, -1.0,  2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(2, false, -1.0, -2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(3, false,  1.0, -2.0));
        std::vector<HEElement<2>*> vertical_rectangle_elements;
        vertical_rectangle_elements.push_back(new HEElement<2>(0, vertical_rectangle_nodes));
        HEVertexMesh<2> vertical_rectangle_mesh(vertical_rectangle_nodes, vertical_rectangle_elements);

        c_vector<double, 2> vertical_rectangle_short_axis = vertical_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_short_axis(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertical_rectangle_short_axis(1), 0.0, 1e-6);

        // Test method with a single rectangular element parallel to the y-axis, centred away from the origin
        std::vector<HENode<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new HENode<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new HENode<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new HENode<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new HENode<2>(3, false, 10.0, 14.0));
        std::vector<HEElement<2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new HEElement<2>(0, far_rectangle_nodes));
        HEVertexMesh<2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 2> far_rectangle_short_axis = far_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_short_axis[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_short_axis[1], 0.0, 1e-4);

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<HENode<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new HENode<2>(0, false,  2.0*0.5*sqrt(3.0) - 1.0*0.5,  2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new HENode<2>(1, false, -2.0*0.5*sqrt(3.0) - 1.0*0.5, -2.0*0.5 + 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new HENode<2>(2, false, -2.0*0.5*sqrt(3.0) + 1.0*0.5, -2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        angled_rectangle_nodes.push_back(new HENode<2>(3, false,  2.0*0.5*sqrt(3.0) + 1.0*0.5,  2.0*0.5 - 1.0*0.5*sqrt(3.0)));
        std::vector<HEElement<2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new HEElement<2>(0, angled_rectangle_nodes));
        HEVertexMesh<2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 2> angled_rectangle_short_axis = angled_rectangle_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(angled_rectangle_short_axis(1), -0.5*sqrt(3.0), 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<HENode<2>*> irregular_nodes;
        irregular_nodes.push_back(new HENode<2>(0, false, 2.7663,  0.1033));
        irregular_nodes.push_back(new HENode<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new HENode<2>(2, false, 1.9433,  0.5313));
        irregular_nodes.push_back(new HENode<2>(3, false, 2.7425,  1.9461));
        irregular_nodes.push_back(new HENode<2>(4, false, 3.3158,  1.5588));
        std::vector<HEElement<2>*> irregular_elements;
        irregular_elements.push_back(new HEElement<2>(0, irregular_nodes));
        HEVertexMesh<2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 2> irregular_short_axis = irregular_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(irregular_short_axis[0], 0.9210, 1e-4);
        TS_ASSERT_DELTA(irregular_short_axis[1], -0.3894, 1e-4);

        // Test with a single trapezoidal element of width 1, lengths 3*sqrt(3.0) and sqrt(3.0), rotated 30 degrees anticlockwise
        std::vector<HENode<2>*> trapezium_nodes;
        trapezium_nodes.push_back(new HENode<2>(0, false,  1.0, 0.0));
        trapezium_nodes.push_back(new HENode<2>(1, false,  2.0, sqrt(3.0)));
        trapezium_nodes.push_back(new HENode<2>(2, false, -2.5, -sqrt(3.0)/2.0));
        trapezium_nodes.push_back(new HENode<2>(3, false, -0.5, -sqrt(3.0)/2.0));
        std::vector<HEElement<2>*> trapezium_elements;
        trapezium_elements.push_back(new HEElement<2>(0, trapezium_nodes));
        HEVertexMesh<2> trapezium_mesh(trapezium_nodes, trapezium_elements);

        c_vector<double, 2> trapezium_short_axis = trapezium_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(trapezium_short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(trapezium_short_axis(1), -0.5*sqrt(3.0), 1e-6);

        // Test method with a single hexagonal element centred at the origin
        std::vector<HENode<2>*> hexagon_nodes;
        for (unsigned i=0; i<6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new HENode<2>(i, false, cos(theta), sin(theta)));
        }
        std::vector<HEElement<2>*> hexagon_elements;
        hexagon_elements.push_back(new HEElement<2>(0, hexagon_nodes));
        HEVertexMesh<2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        // Since the element is a regular polygon, the short axis is a random vector, so test the random seed
        c_vector<double, 2> hexagon_short_axis = hexagon_mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(norm_2(hexagon_short_axis), 1.0, 1e-6);
        TS_ASSERT_DELTA(hexagon_short_axis(0), 0.5488, 1e-4);
        TS_ASSERT_DELTA(hexagon_short_axis(1), 0.8359, 1e-4);
    }

    void TestGetElongationShapeFactorOfElement()
    {
        // Test method with a single square element
        std::vector<HENode<2>*> square_nodes;
        square_nodes.push_back(new HENode<2>(0, false,  1.0,  1.0));
        square_nodes.push_back(new HENode<2>(1, false, -1.0,  1.0));
        square_nodes.push_back(new HENode<2>(2, false, -1.0, -1.0));
        square_nodes.push_back(new HENode<2>(3, false,  1.0, -1.0));
        std::vector<HEElement<2>*> square_elements;
        square_elements.push_back(new HEElement<2>(0, square_nodes));
        HEVertexMesh<2> square_mesh(square_nodes, square_elements);

        // By symmetry, the two eigenvalues are equal, so the elongation shape factor equals one
        double square_elongation_shape_factor = square_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(square_elongation_shape_factor, 1.0, 1e-6);

        // We should obtain the same result regardless of the square's size
        square_mesh.Scale(5.0, 5.0);
        square_elongation_shape_factor = square_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(square_elongation_shape_factor, 1.0, 1e-6);

        // Test method with a single rectangular element
        std::vector<HENode<2>*> vertical_rectangle_nodes;
        vertical_rectangle_nodes.push_back(new HENode<2>(0, false,  1.0,  2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(1, false, -1.0,  2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(2, false, -1.0, -2.0));
        vertical_rectangle_nodes.push_back(new HENode<2>(3, false,  1.0, -2.0));
        std::vector<HEElement<2>*> vertical_rectangle_elements;
        vertical_rectangle_elements.push_back(new HEElement<2>(0, vertical_rectangle_nodes));
        HEVertexMesh<2> vertical_rectangle_mesh(vertical_rectangle_nodes, vertical_rectangle_elements);

        /*
         * For a rectangle of width a (parallel to the x axis) and height b (parallel to the
         * y axis), the second moments of area are given by J_xx = a*b^3/12, J_yy = b*a^3/12
         * and J_xy = 0. Therefore the two eigenvalues are given by J_xx and J_yy and if b>a,
         * the elongation shape factor should be equal to sqrt(a*b^3/b*a^3) = b/a.
         */
        double vertical_rectangle_elongation_shape_factor = vertical_rectangle_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_elongation_shape_factor, 4.0/2.0, 1e-6);

        // We should obtain the same result regardless of the orientation of the rectangle
        vertical_rectangle_mesh.Scale(2.0, 0.5);
        vertical_rectangle_elongation_shape_factor = vertical_rectangle_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_elongation_shape_factor, 4.0/2.0, 1e-6);

        // Check the formula is correct for a different aspect ratio
        vertical_rectangle_mesh.Scale(3.0, 1.0);
        vertical_rectangle_elongation_shape_factor = vertical_rectangle_mesh.GetElongationShapeFactorOfElement(0);
        TS_ASSERT_DELTA(vertical_rectangle_elongation_shape_factor, 12.0/2.0, 1e-6);
    }

};
#endif /*TESTDCELMESH_HPP_*/


