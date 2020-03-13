#ifndef TESTDCELMESH_HPP_
#define TESTDCELMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "HEElement.hpp"
#include "HalfEdge.hpp"
#include "HEVertex.hpp"
#include "VertexMesh.hpp"
#include "HEVertexMesh.hpp"
#include "FakePetscSetup.hpp"
class TestHEVertexMesh : public CxxTest::TestSuite
{
public:
    void TestBasicMesh()
    {
        std::vector<HEVertex<2>* > vertices;
        vertices.push_back(new HEVertex<2>(0,false,0.0,0.0));
        vertices.push_back(new HEVertex<2>(1,false,1.0,0.0));
        vertices.push_back(new HEVertex<2>(2,false,1.0,1.0));
        vertices.push_back(new HEVertex<2>(3,false,0.0,1.0));

        vertices.push_back(new HEVertex<2>(4,false,2.0,0.0));
        vertices.push_back(new HEVertex<2>(5,false,2.0,1.0));

        std::vector<HEVertex<2>* > element_1_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        std::vector<HEVertex<2>* > element_2_vertices = {vertices[1], vertices[4], vertices[5], vertices[2]};
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
        c_vector<double, 2> common_edge_vertex = common_edge->GetTargetVertex()->rGetLocation();
        c_vector<double, 2> common_edge_twin_vertex = common_edge->GetTwinHalfEdge()->GetTargetVertex()->rGetLocation();
        c_vector<double, 2> common_edge_prev_edge_vertex
        = common_edge->GetPreviousHalfEdge()->GetTargetVertex()->rGetLocation();
        c_vector<double, 2> common_edge_twin_prev_edge_vertex
        = common_edge->GetTwinHalfEdge()->GetPreviousHalfEdge()->GetTargetVertex()->rGetLocation();
        TS_ASSERT_DELTA(common_edge_vertex(0), common_edge_twin_prev_edge_vertex(0),1e-12);
        TS_ASSERT_DELTA(common_edge_vertex(1), common_edge_twin_prev_edge_vertex(1),1e-12);
        TS_ASSERT_DELTA(common_edge_prev_edge_vertex(0), common_edge_twin_vertex(0),1e-12);
        TS_ASSERT_DELTA(common_edge_prev_edge_vertex(1), common_edge_twin_vertex(1),1e-12);

        c_vector<double, 2> element_1_centroid = mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(element_1_centroid(0),0.5,1e-12);
        TS_ASSERT_DELTA(element_1_centroid(1),0.5,1e-12);
    }

    void TestBasic2DMesh()
    {
        // Create a 2D mesh comprising seven nodes and two elements
        std::vector<HEVertex<2>*> nodes_2d;
        nodes_2d.push_back(new HEVertex<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new HEVertex<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new HEVertex<2>(2, false, 1.5, 1.0));
        nodes_2d.push_back(new HEVertex<2>(3, false, 1.0, 2.0));
        nodes_2d.push_back(new HEVertex<2>(4, false, 0.0, 1.0));
        nodes_2d.push_back(new HEVertex<2>(5, false, 2.0, 0.0));

        std::vector<std::vector<HEVertex<2>*> > nodes_elements_2d(2);
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
        TS_ASSERT_EQUALS(mesh_2d.GetElement(1)->GetHalfEdge()->GetNextHalfEdge()->GetTargetVertex()->GetIndex(), 2u);

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


    }
};
#endif /*TESTDCELMESH_HPP_*/


