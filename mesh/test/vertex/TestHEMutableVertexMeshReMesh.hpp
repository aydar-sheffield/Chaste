/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTHEMUTABLEVERTEXMESHREMESH_HPP_
#define TESTHEMUTABLEVERTEXMESHREMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"
#include "FileComparison.hpp"
#include "Warnings.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#include "HEMutableVertexMesh.hpp"
class TestHEMutableVertexMeshReMesh : public CxxTest::TestSuite
{
public:
    void TestCollapseEdge()
    {
        /*
         * Create a mesh comprising a single triangular element, as shown below.
         * We will test that the nodes marked with an x are merged correctly.
         *
         *      /|
         *     / |
         *    /  |
         *   /   |
         *  /    |
         *  --xx-
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true, 1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true, 0.4, 0.0));
        nodes.push_back(new HENode<2>(4, true, 0.6, 0.0));

        unsigned node_indices_elem_0[5] = {0, 3, 4, 1, 2};
        std::vector<HENode<2>*> nodes_elem_0;
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Merge nodes 3 and 4
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdge(1);

        HalfEdge<2>* half_edge = full_edge->at(0)->GetElement()==nullptr ? full_edge->at(1) : full_edge->at(0);
        assert(half_edge->GetTargetNode()->GetIndex()==4);
        vertex_mesh.IdentifySwapType(*full_edge);

        // Test the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        // Test the correct nodes are boundary nodes
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Test the merged node is in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[1], 0.0, 1e-3);

        // Test the elements own the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 2u);

        // Test the element's area and perimeter are computed correctly
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 2+sqrt(2.0), 1e-6);
    }

    void TestCollapseEdgeWhenLowIndexNodeMustBeAddedToElement()
    {
        /**
         * Create a mesh comprising two square elements, as shown below. We will test that the
         * nodes marked with an x are merged correctly. We will test node merging in the case
         * where, when the elements previously containing the high-index node are updated to
         * contain the low-index node, at least one of these elements did not already contain
         * the low-index node.
         *
         *   -----x-x---
         *  |     |     |
         *  |     |     |
         *   ----- -----
         *
         * \todo I think this should be a T1 swap (see #1263)
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.00, 0.00));
        nodes.push_back(new HENode<2>(1, true, 1.00, 0.00));
        nodes.push_back(new HENode<2>(2, true, 2.00, 0.00));
        nodes.push_back(new HENode<2>(3, true, 2.00, 1.00));
        nodes.push_back(new HENode<2>(4, true, 1.01, 1.00));
        nodes.push_back(new HENode<2>(5, true, 1.00, 1.00));
        nodes.push_back(new HENode<2>(6, true, 0.00, 2.00));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 5, 6};
        unsigned node_indices_elem_1[5] = {1, 2, 3, 4, 5};
        for (unsigned i=0; i<5; i++)
        {
            if (i < 4)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Merge nodes 4 and 5
        assert(vertex_elements[1]->GetHalfEdge(3)->GetTargetNode()->GetIndex()==5);
        assert(vertex_elements[1]->GetHalfEdge(3)->GetOriginNode()->GetIndex()==4);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[1]->GetHalfEdge(3));

        vertex_mesh.IdentifySwapType(*full_edge);

        // Test the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Test the correct nodes are boundary nodes
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Test that the moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 1.005, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 1.0, 1e-8);

        // Test the elements own the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        unsigned node_indices_element_0[4] = {0, 1, 4, 5};
        unsigned node_indices_element_1[4] = {1, 2, 3, 4};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }
    }

    void TestPerformT1SwapAndIdentifySwapType()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *  _____
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true,  1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true,  1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true,  0.0, 1.0));
        nodes.push_back(new HENode<2>(4, false, 0.5, 0.4));
        nodes.push_back(new HENode<2>(5, false, 0.5, 0.6));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {2, 5, 4, 1};
        unsigned node_indices_elem_2[3] = {1, 4, 0};
        unsigned node_indices_elem_3[4] = {0, 4, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        // Perform a T1 swap on nodes 4 and 5
        assert(vertex_elements[3]->GetHalfEdge(1)->GetTargetNode()->GetIndex()==5);
        assert(vertex_elements[3]->GetHalfEdge(1)->GetOriginNode()->GetIndex()==4);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[3]->GetHalfEdge(1));

        vertex_mesh.IdentifySwapType(*full_edge);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[4] = {2, 3, 5, 4};
        unsigned node_indices_element_1[3] = {2, 4, 1};
        unsigned node_indices_element_2[4] = {1, 4, 5, 0};
        unsigned node_indices_element_3[3] = {0, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test T1 swap location tracking
        std::vector< c_vector<double, 2> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestPerformT1SwapOnBoundary()
    {
        /*
         * Create a mesh comprising six nodes contained in three elements such that all nodes are
         * boundary nodes, as shown below. We will test that that a T1 swap is correctly implemented.
         *  _____
         * |\   /
         * | \ /
         * |  |
         * | / \
         * |/___\
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true, 1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true, 0.0, 1.0));
        nodes.push_back(new HENode<2>(4, true, 0.5, 0.4));
        nodes.push_back(new HENode<2>(5, true, 0.5, 0.6));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[3] = {1, 4, 0};
        unsigned node_indices_elem_2[4] = {0, 4, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        // Perform a T1 swap on nodes 5 and 4 (this way round to ensure coverage of boundary node tracking)
        // Perform a T1 swap on nodes 4 and 5
        assert(vertex_elements[2]->GetHalfEdge(1)->GetTargetNode()->GetIndex()==5);
        assert(vertex_elements[2]->GetHalfEdge(1)->GetOriginNode()->GetIndex()==4);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[2]->GetHalfEdge(1));

        vertex_mesh.IdentifySwapType(*full_edge);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test that each element contains the correct number nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[4] = {2, 3, 5, 4};
        unsigned node_indices_element_1[4] = {1, 4, 5, 0};
        unsigned node_indices_element_2[4] = {0, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = (i!=5);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }
    }

    void TestPerformT1SwapOnBoundary2()
    {
        /*
         * Create a mesh comprising six nodes contained in three elements such that all but one node
         * are boundary nodes, as shown below. We will test that that a T1 swap is correctly implemented.
         *
         * |\   /|
         * | \ / |
         * |  |  |
         * | / \ |
         * |/___\|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true,  1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true,  1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true,  0.0, 1.0));
        nodes.push_back(new HENode<2>(4, false, 0.5, 0.4));
        nodes.push_back(new HENode<2>(5, true,  0.5, 0.6));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_1[3] = {1, 4, 0};
        unsigned node_indices_elem_2[4] = {0, 4, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Perform a T1 swap on nodes 5 and 4 (this way round to ensure coverage of boundary node tracking)
        assert(vertex_elements[2]->GetHalfEdge(1)->GetTargetNode()->GetIndex()==5);
        assert(vertex_elements[2]->GetHalfEdge(1)->GetOriginNode()->GetIndex()==4);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[2]->GetHalfEdge(1));
        vertex_mesh.IdentifySwapType(*full_edge);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test that each element contains the correct number nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = {1, 2, 4};
        unsigned node_indices_element_1[4] = {1, 4, 5, 0};
        unsigned node_indices_element_2[3] = {0, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestPerformT1SwapWhenVoidForms()
    {
        /*
         * Create a mesh containing six nodes containing in two elements. We will test that
         * a T1 swap is correctly performed in the case where a void forms as a result of
         * the rearrangement, as shown below.
         *
         * |\   /|     |\      /|
         * | \ / |     | \    / |
         * |  |  |  => | /    \ |
         * | / \ |     |/      \|
         * |/   \|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true, 1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true, 0.0, 1.0));
        nodes.push_back(new HENode<2>(4, true, 0.5, 0.4));
        nodes.push_back(new HENode<2>(5, true, 0.5, 0.6));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 4, 5, 3};
        unsigned node_indices_elem_1[4] = {4, 1, 2, 5};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        // Perform a T1 swap on nodes 5 and 4.
        assert(vertex_elements[0]->GetHalfEdge(1)->GetTargetNode()->GetIndex()==5);
        assert(vertex_elements[0]->GetHalfEdge(1)->GetOriginNode()->GetIndex()==4);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[0]->GetHalfEdge(1));
        vertex_mesh.IdentifySwapType(*full_edge);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test that each element contains the correct number of nodes following the rearrangement
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[3] = {0, 5, 3};
        unsigned node_indices_element_1[3] = {4, 1, 2};
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);

        // Test that the correct nodes are labelled as boundary nodes following the rearrangement
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestPerformT1SwapExceptions()
    {
        /*
         * Create a mesh comprising six nodes containing in two triangle and two rhomboid elements,
         * where two nodes (those with indices 4 and 5) have the same location. We will test that
         * trying to perform a T1 swap on these nodes throws the correct exception.
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, false, 1.0, 1.0));
        nodes.push_back(new HENode<2>(3, false, 0.0, 1.0));
        nodes.push_back(new HENode<2>(4, false, 0.5, 0.5));
        nodes.push_back(new HENode<2>(5, false, 0.5, 0.5));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {2, 5, 4, 1};
        unsigned node_indices_elem_2[3] = {1, 4, 0};
        unsigned node_indices_elem_3[4] = {0, 4, 5, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1*2.0/1.5);

        // Test that trying to perform a T1 swap on nodes 4 and 5 throws the correct exception
        assert(vertex_elements[3]->GetHalfEdge(1)->GetTargetNode()->GetIndex()==5);
        assert(vertex_elements[3]->GetHalfEdge(1)->GetOriginNode()->GetIndex()==4);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[3]->GetHalfEdge(1));
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(*full_edge), "Nodes are too close together, this shouldn't happen");
    }

    void TestDoNotPerforT1SwapWithRemovingEdgeFromTriangularElement()
    {
        /**
         * In this test we check that a T1 swap does not occur if one of the elements is triangular
         * and would loose an edge by swapping nodes. The mesh looks like this
         *
         *       ______________
         *      |\             |
         *      | \ _________  |
         *      |  |          \| ...where the funny shaped element in the middle is supposed to be
         *      |  |_________ /|    a very long triangle that has the third vertex on the right hand boundary.
         *      | /            |
         *      |/_____________|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true,  2.0, 0.0));
        nodes.push_back(new HENode<2>(2, true,  2.0, 2.0));
        nodes.push_back(new HENode<2>(3, true,  0.0, 2.0));
        nodes.push_back(new HENode<2>(4, false, 0.2, 0.95));
        nodes.push_back(new HENode<2>(5, true, 2.0, 1.0));
        nodes.push_back(new HENode<2>(6, false, 0.2, 1.05));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[4] = {0, 1, 5, 4};
        unsigned node_indices_elem_1[4] = {5, 2, 3, 6};
        unsigned node_indices_elem_2[4] = {0, 4, 6, 3};
        unsigned node_indices_elem_3[3] = { 4, 5, 6};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
            }
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Ensure that the inner edge will be considered for a swap
        vertex_mesh.SetCellRearrangementThreshold(0.11);

        // Check for T1 swaps and carry them out if allowed - the short edge should not swap!
        vertex_mesh.CheckForSwapsFromShortEdges();

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);

        // Test that each element still contains the correct nodes following the rearrangement
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_elem_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_elem_1[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_elem_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_elem_3[i]);
            }
        }
    }

    void TestExceptionForVoidRemovalWithRemovingEdgeFromTriangularElement()
    {
        /**
         * In this test we check that void removal does not occur if one of the adjacent elements is triangular
         * and would loose an edge by swapping nodes. The code should throw and exception in this case.
         * The mesh looks like this
         *
         *       ______________./This corner is not a node.
         *      |\      1      |
         *      | \ _________  |
         *      |  |   void   \| ...where elements 1, and 2 are triangles that share the right hand vertex
         *      |  |_________ /|    with the triangular void in the middle.
         *      | /     2      |
         *      |/_____________|.This corner is not a node either.
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true,  0.0, 2.0));
        nodes.push_back(new HENode<2>(2, true, 0.2, 0.95));
        nodes.push_back(new HENode<2>(3, true, 2.0, 1.0));
        nodes.push_back(new HENode<2>(4, true, 0.2, 1.05));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[4] = {0, 2, 4, 1};
        unsigned node_indices_elem_1[3] = {1, 4, 3};
        unsigned node_indices_elem_2[3] = {0, 3, 2};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            }
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Ensure that the inner edge will be considered for a swap
        vertex_mesh.SetCellRearrangementThreshold(0.11);

        // Check for possible swaps and carry them out if allowed - the short edge should not swap and
        // the void should not be removed!
        TS_ASSERT_THROWS_THIS(vertex_mesh.CheckForSwapsFromShortEdges(),
                              "Triangular element next to triangular void, not implemented yet.");
    }

    void TestPerformT2Swap()
    {
        /*
         * Create a mesh comprising six nodes contained in three trapezium element and
         * a central triangle element, as shown below. We will test that a T2 swap
         * correctly removes the triangle element from the mesh.
         *
         *      /|\
         *     / | \
         *    /  |  \    (the triangular element has index zero)
         *   /2 /_\ 1\
         *  /  /   \  \
         * /__/__3__\__\
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true, 0.5, 0.5));
        nodes.push_back(new HENode<2>(3, false, 0.4, 0.25));
        nodes.push_back(new HENode<2>(4, false, 0.6, 0.25));
        nodes.push_back(new HENode<2>(5, false, 0.5, 0.3));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        unsigned node_indices_elem_3[4] = {0, 1, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Perform a T2 swap on the central triangle element
        HEElement<2>* p_element_0 = vertex_mesh.GetElement(0);
        c_vector<double, 2> centroid_of_element_0_before_swap = vertex_mesh.GetCentroidOfElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFullEdges(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFullEdges(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 7u);

        for (unsigned j=1; j<4; j++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 3u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(0)->GetIndex(), j%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(1)->GetIndex(), (j+1)%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(2)->GetIndex(), 6u);
        }

        // Test boundary property of nodes. All are boundary nodes except node 3.
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            bool expected_boundary_node = (i!=3);
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), expected_boundary_node);
        }

        // Test the location of the new node:
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], centroid_of_element_0_before_swap[0], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], centroid_of_element_0_before_swap[1], 1e-10);

        // Test the tracking of the T2 swap location:
        TS_ASSERT_DELTA(vertex_mesh.GetLastT2SwapLocation()[0], centroid_of_element_0_before_swap[0], 1e-10);
        TS_ASSERT_DELTA(vertex_mesh.GetLastT2SwapLocation()[1], centroid_of_element_0_before_swap[1], 1e-10);
    }

    void TestPerformT2SwapWithBoundaryNodes()
    {
        /*
         * Create a mesh comprising six nodes contained in two trapezium elements
         * and one triangle element, as shown below. We will test that a T2 swap
         * is performed correctly when boundary nodes are involved.
         *
         *       /|\
         *      / | \
         *     /  |  \    (the triangular element has index zero)
         *    /2 /_\ 1\
         *   /  /   \  \
         *  /__/     \__\
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true, 0.5, 0.5));
        nodes.push_back(new HENode<2>(3, true, 0.4, 0.25));
        nodes.push_back(new HENode<2>(4, true, 0.6, 0.25));
        nodes.push_back(new HENode<2>(5, false, 0.5, 0.3));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFullEdges(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Perform a T2 swap on the central triangle element
        HEElement<2>* p_element_0 = vertex_mesh.GetElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFullEdges(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFullEdges(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 7u);

        for (unsigned j=1; j<3; j++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNumNodes(), 3u);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(0)->GetIndex(), j%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(1)->GetIndex(), (j+1)%3);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(j)->GetNode(2)->GetIndex(), 6u);
        }
        // Test boundary property of nodes. All are boundary nodes.
        for (unsigned i=0; i<vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Make five nodes to assign to two elements
        std::vector<HENode<2>*> nodes2;
        nodes2.push_back(new HENode<2>(0, true, 1.0, 0.0));
        nodes2.push_back(new HENode<2>(1, true, 0.5, 0.5));
        nodes2.push_back(new HENode<2>(2, true, 0.4, 0.25));
        nodes2.push_back(new HENode<2>(3, true, 0.6, 0.25));
        nodes2.push_back(new HENode<2>(4, true, 0.5, 0.3));

        /*
         *  Make one trapezium element with a central triangular element out of these nodes
         *
         *       |\
         *       | \
         *       |  \
         *      /_\ 1\   Triangular element is element zero
         *         \_ \
         *           \_\
         *
         */

        // Triangle element
        std::vector<HENode<2>*> nodes2_elem_0;
        nodes2_elem_0.push_back(nodes2[2]);
        nodes2_elem_0.push_back(nodes2[3]);
        nodes2_elem_0.push_back(nodes2[4]);

        // Trapezium
        std::vector<HENode<2>*> nodes2_elem_1;
        nodes2_elem_1.push_back(nodes2[0]);
        nodes2_elem_1.push_back(nodes2[1]);
        nodes2_elem_1.push_back(nodes2[4]);
        nodes2_elem_1.push_back(nodes2[3]);

        std::vector<HEElement<2>*> vertex_elements2;
        vertex_elements2.push_back(new HEElement<2>(0, nodes2_elem_0));
        vertex_elements2.push_back(new HEElement<2>(1, nodes2_elem_1));

        // Make a vertex mesh
        HEMutableVertexMesh<2> vertex_mesh2(nodes2, vertex_elements2);

        TS_ASSERT_EQUALS(vertex_mesh2.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumFullEdges(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumNodes(), 5u);

        // Perform a T2 swap on the middle triangle element
        p_element_0 = vertex_mesh2.GetElement(0);
        vertex_mesh2.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh2.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumFullEdges(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumNodes(), 3u);

        TS_ASSERT_EQUALS(vertex_mesh2.GetNumAllElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumAllFullEdges(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumAllNodes(), 6u);

        TS_ASSERT_EQUALS(vertex_mesh2.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetElement(1)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetElement(1)->GetNode(2)->GetIndex(), 5u);

        // Test boundary property of nodes. All are boundary nodes.
        for (unsigned i=0; i<vertex_mesh2.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh2.GetNode(i)->IsBoundaryNode(), true);
        }
    }

    void TestT2SwapsDontOccurWithTriangularNeighbours()
    {
        // Make 6 nodes to assign to four elements
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, false, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, false, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, false, 0.5, 0.5));
        nodes.push_back(new HENode<2>(3, false, 0.4, 0.25));
        nodes.push_back(new HENode<2>(4, false, 0.6, 0.25));
        nodes.push_back(new HENode<2>(5, false, 0.5, 0.3));

        // Make two triangles and two trapezium elements out of these nodes
        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[3] = {2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        unsigned node_indices_elem_3[4] = {0, 1, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements, 0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);

        // Attempt to perform a T2 swap on the middle triangle element
        HEElement<2>* p_element_0 = vertex_mesh.GetElement(0);
        TS_ASSERT_THROWS_THIS( vertex_mesh.PerformT2Swap(*p_element_0),
                               "One of the neighbours of a small triangular element is also a triangle - "
                               "dealing with this has not been implemented yet" );
    }

    void TestPerformT2SwapWithRosettes()
    {
        /* Create a mesh containing a smaller triangular element, each of whose nodes are
         * 'rosette' nodes. Test that a T2 swap correctly removes the triangular element
         * from the mesh.
         *  _________
         *  |\      |
         *  | \     |
         *  | |_\___|
         *  | / |   |
         *  |/__|___|
         */

        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0,  0.0));
        nodes.push_back(new HENode<2>(1, true,  1.0,  0.0));
        nodes.push_back(new HENode<2>(2, true,  2.0,  0.0));
        nodes.push_back(new HENode<2>(3, false, 0.99, 1.0));
        nodes.push_back(new HENode<2>(4, false, 1.0,  1.0));
        nodes.push_back(new HENode<2>(5, false, 2.0,  1.0));
        nodes.push_back(new HENode<2>(6, false, 0.99, 1.01));
        nodes.push_back(new HENode<2>(7, false, 0.0,  2.0));
        nodes.push_back(new HENode<2>(8, false, 2.0,  2.0));

        std::vector<HENode<2>*> nodes_elem_0;
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[3]);

        std::vector<HENode<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);

        std::vector<HENode<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[6]);
        nodes_elem_2.push_back(nodes[7]);

        std::vector<HENode<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[3]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[6]);

        std::vector<HENode<2>*> nodes_elem_4;
        nodes_elem_4.push_back(nodes[4]);
        nodes_elem_4.push_back(nodes[5]);
        nodes_elem_4.push_back(nodes[8]);
        nodes_elem_4.push_back(nodes[7]);
        nodes_elem_4.push_back(nodes[6]);

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));
        vertex_elements.push_back(new HEElement<2>(4, nodes_elem_4));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Perform a T2 swap on the central triangle element
        HEElement<2>* p_element_3 = vertex_mesh.GetElement(3);
        c_vector<double, 2> centroid_of_element_0_before_swap = vertex_mesh.GetCentroidOfElement(3);
        vertex_mesh.PerformT2Swap(*p_element_3);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFullEdges(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllFullEdges(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllNodes(), 10u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 9u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 9u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(0), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(2), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(0), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(2), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNodeGlobalIndex(3), 7u);
    }

    void TestReMeshForT1Swaps()
    {
        /*
         * Read in a vertex mesh that contains several pairs of nodes that are close enough for
         * T1 swaps to be performed, as shown below. The mesh consists of six elements and all
         * T1 swaps are performed on all horizontal edges. We will test that the ReMesh() method
         * correctly performs T1 swaps for internal and boundary elements, and correctly updates
         * which nodes are labelled as boundary nodes.
         *
         *      /\    /\
         *     /  \__/  \
         *    /   /  \   \
         *    \__/\__/\__/
         *    /  \/  \/  \
         *    \   \__/   /
         *     \  /  \  /
         *      \/    \/
         */
        VertexMeshReader<2,2> mesh_reader("cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1");
        MutableVertexMesh<2,2> vertex_mesh;

        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        HEMutableVertexMesh<2> he_mesh;
        he_mesh.ConvertFromVertexMesh(&vertex_mesh);
        he_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(he_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(he_mesh.GetNumNodes(), 22u);

        //Test whether Vertex mesh and HE Mesh are the same...
        //Tests whether HEMutableVertexMesh::ConvertFromVertexMesh() performs as intended
        for (unsigned int i=0; i<8; ++i)
        {
            VertexElement<2,2>* v_element = vertex_mesh.GetElement(i);
            HEElement<2>* he_element = he_mesh.GetElement(i);
            const unsigned int n_v_nodes = v_element->GetNumNodes();
            const unsigned int n_he_nodes = he_element->GetNumNodes();
            assert(n_v_nodes == n_he_nodes);
            for (unsigned int j=0; j<n_v_nodes; ++j)
            {
                TS_ASSERT_EQUALS(v_element->GetNode(j)->GetIndex(), he_element->GetNode(j)->GetIndex());
                TS_ASSERT_DELTA(v_element->GetNode(j)->rGetLocation()[0], he_element->GetNode(j)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(v_element->GetNode(j)->rGetLocation()[1], he_element->GetNode(j)->rGetLocation()[1],1e-6);
            }
        }

        // Calls ReMesh() to identify and perform any T1 swaps
        he_mesh.ReMesh();
        vertex_mesh.ReMesh();

        //Test whether ReMesh methods for HEMutableVertexMesh and MutableVertexMesh perform exactly the same
        for (unsigned int i=0; i<8; ++i)
        {
            VertexElement<2,2>* v_element = vertex_mesh.GetElement(i);
            HEElement<2>* he_element = he_mesh.GetElement(i);
            const unsigned int n_v_nodes = v_element->GetNumNodes();
            const unsigned int n_he_nodes = he_element->GetNumNodes();
            assert(n_v_nodes == n_he_nodes);

            //Shift vertex ordering if necessary:
            //i.e. if Vertex Element node order is 1 2 3 4
            //and HEElement order is 2 3 4 1, then elements are the same.
            for (unsigned int j=0; j<n_v_nodes; ++j)
            {
                if (he_element->GetNode(j)->GetIndex()==v_element->GetNode((j+1)%n_v_nodes)->GetIndex())
                {
                    he_element->SetHalfEdge(he_element->GetHalfEdge()->GetPreviousHalfEdge());
                    break;
                }
                if (he_element->GetNode(j)->GetIndex()==v_element->GetNode((j-1+n_v_nodes)%n_v_nodes)->GetIndex())
                {
                    he_element->SetHalfEdge(he_element->GetHalfEdge()->GetNextHalfEdge());
                    break;
                }
            }

            for (unsigned int j=0; j<n_v_nodes; ++j)
            {
                TS_ASSERT_EQUALS(v_element->GetNode(j)->GetIndex(),he_element->GetNode(j)->GetIndex());
                TS_ASSERT_DELTA(v_element->GetNode(j)->rGetLocation()[0], he_element->GetNode(j)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(v_element->GetNode(j)->rGetLocation()[1], he_element->GetNode(j)->rGetLocation()[1],1e-6);
            }
        }

        TS_ASSERT_EQUALS(he_mesh.GetNumElements(), 8u);
        TS_ASSERT_EQUALS(he_mesh.GetNumNodes(), 22u);

        std::string dirname = "TestVertexMeshReMesh";
        std::string mesh_filename = "vertex_remesh_T1";

        MutableVertexMesh<2,2>* vertex_mesh_to_write = he_mesh.ConvertToMutableVertexMesh();

        //Test if mesh conversion from HEMesh to VertexMesh works as intended
        for (unsigned int i=0; i<8; ++i)
        {
            VertexElement<2,2>* v_element = vertex_mesh.GetElement(i);
            VertexElement<2,2>* v2_element = vertex_mesh_to_write->GetElement(i);
            const unsigned int n_v_nodes = v_element->GetNumNodes();
            const unsigned int n_v2_nodes = v2_element->GetNumNodes();
            assert(n_v_nodes == n_v2_nodes);

            for (unsigned int j=0; j<n_v_nodes; ++j)
            {
                TS_ASSERT_EQUALS(v_element->GetNode(j)->GetIndex(),v2_element->GetNode(j)->GetIndex());
                TS_ASSERT_DELTA(v_element->GetNode(j)->rGetLocation()[0], v2_element->GetNode(j)->rGetLocation()[0],1e-6);
                TS_ASSERT_DELTA(v_element->GetNode(j)->rGetLocation()[1], v2_element->GetNode(j)->rGetLocation()[1],1e-6);
            }
        }

        // Save the mesh data using mesh writers
        VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
        mesh_writer.WriteFilesUsingMesh(*vertex_mesh_to_write);

        // Check the positions are updated correctly
        OutputFileHandler handler("TestVertexMeshReMesh", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T1.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_remesh_T1.cell";

        FileComparison comparer1(results_file1, "cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1_after_remesh.node");
        TS_ASSERT(comparer1.CompareFiles());
        FileComparison comparer2(results_file2, "cell_based/test/data/TestMutableVertexMesh/vertex_remesh_T1_after_remesh.cell");
        TS_ASSERT(comparer2.CompareFiles());
    }

    void TestReMeshExceptionWhenNonBoundaryNodesAreContainedOnlyInTwoElements()
    {
        /*
         * Create a mesh comprising six nodes contained in two elements, as shown below. We will test
         * that when we attempt to mesh the nodes marked x, the correct exception is thrown.
         *   ______
         *  |     /|
         *  |    / |
         *  |   x  |
         *  |  x   |
         *  | /    |
         *  |/_____|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.00, 0.00));
        nodes.push_back(new HENode<2>(1, true,  1.00, 0.00));
        nodes.push_back(new HENode<2>(2, true,  1.00, 1.00));
        nodes.push_back(new HENode<2>(3, true,  0.00, 1.00));
        nodes.push_back(new HENode<2>(4, false, 0.49, 0.49));
        nodes.push_back(new HENode<2>(5, false, 0.51, 0.51));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = {0, 1, 2, 5, 4};
        unsigned node_indices_elem_1[5] = {0, 4, 5, 2, 3};
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "There are non-boundary nodes contained only in two elements; something has gone wrong.");
    }

    void TestIdentifySwapTypeExceptionWhenBoundaryNodeIsContainedInThreeElements()
    {
        /*
         * Create a mesh as shown below, where the two nodes marked with an x are to be merged.
         * We will test that the correct exception is thrown when IdentifySwapType() is called.
         *             _____
         *           /      |
         *          /       |
         *   -----xx--------|
         *  |       \       |
         *  |________\______|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.00, 0.0));
        nodes.push_back(new HENode<2>(1, true, 2.00, 0.0));
        nodes.push_back(new HENode<2>(2, true, 3.00, 0.0));
        nodes.push_back(new HENode<2>(3, true, 3.00, 1.0));
        nodes.push_back(new HENode<2>(4, true, 3.00, 2.0));
        nodes.push_back(new HENode<2>(5, true, 2.00, 2.0));
        nodes.push_back(new HENode<2>(6, true, 1.00, 1.0));
        nodes.push_back(new HENode<2>(7, true, 0.99, 1.0));
        nodes.push_back(new HENode<2>(8, true, 0.00, 1.0));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[5] = {0, 1, 6, 7, 8};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 6};
        unsigned node_indices_elem_2[4] = {3, 4, 5, 6};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
        }
        nodes_elem_0.push_back(nodes[node_indices_elem_0[4]]);

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 9u);

        assert(vertex_elements[0]->GetHalfEdge(2)->GetTargetNode()->GetIndex()==7);
        assert(vertex_elements[0]->GetHalfEdge(2)->GetOriginNode()->GetIndex()==6);
        FullEdge<2>* full_edge = vertex_mesh.GetFullEdgeFromHalfEdge(vertex_elements[0]->GetHalfEdge(2));
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(*full_edge), "There is a boundary node contained in three elements something has gone wrong.");
    }

    void TestReMeshExceptionWhenNonBoundaryNodeIsContainedOnlyInTwoElements()
    {
        /*
         * Create a mesh as shown below, where the two nodes marked with an x are to be merged.
         * We will test that the correct exception is thrown when ReMesh() is called.
         *
         * |\   /|
         * | \ / |
         * |  x  |
         * |  x  |
         * |__|__|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true, 1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true, 0.0, 1.0));
        nodes.push_back(new HENode<2>(4, false, 0.5, 0.49));
        nodes.push_back(new HENode<2>(5, true, 0.5, 0.51));
        nodes.push_back(new HENode<2>(6, true, 0.5, 0.0));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = {1, 2, 5, 4, 6};
        unsigned node_indices_elem_1[5] = {0, 6, 4, 5, 3};
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "There is a non-boundary node contained only in two elements; something has gone wrong.");
    }

    void TestAnotherReMeshExceptionWhenNonBoundaryNodesAreContainedOnlyInTwoElements()
    {
        /*
         * Create a mesh as shown below, where the two central nodes marked with an x are to be merged. We will
         * test that the correct exception is thrown when ReMesh() is called. Note that the extra node at the
         * top of the mesh is required to stop the element from containing only three nodes, otherwise the
         * ReMesh() method would not call IdentifySwapType().
         *
         *  __x__
         * |\   /|
         * | \ / |
         * |  x  |
         * |  x  |
         * |__|__|
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.00));
        nodes.push_back(new HENode<2>(1, true,  1.0, 0.00));
        nodes.push_back(new HENode<2>(2, true,  1.0, 1.00));
        nodes.push_back(new HENode<2>(3, true,  0.0, 1.00));
        nodes.push_back(new HENode<2>(4, false, 0.5, 0.49));
        nodes.push_back(new HENode<2>(5, false, 0.5, 0.51));
        nodes.push_back(new HENode<2>(6, true,  0.5, 0.00));
        nodes.push_back(new HENode<2>(7, true,  0.5, 1.00));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2;
        unsigned node_indices_elem_0[5] = {1, 2, 5, 4, 6};
        unsigned node_indices_elem_1[5] = {0, 6, 4, 5, 3};
        unsigned node_indices_elem_2[4] = {2, 3, 5, 7};
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
        }
        nodes_elem_0.push_back(nodes[node_indices_elem_0[4]]);
        nodes_elem_1.push_back(nodes[node_indices_elem_1[4]]);

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);
        vertex_mesh.SetCellRearrangementThreshold(0.1);

        TS_ASSERT_THROWS_THIS(vertex_mesh.ReMesh(), "There are non-boundary nodes contained only in two elements; something has gone wrong.");
    }

    void TestPerformIntersectionSwap()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements.
         * We will test that when a node is moved to overlap with an element, it is correctly
         * found and dealt with by the CheckForIntersections() method.
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true,  1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true,  1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true,  0.0, 1.0));
        nodes.push_back(new HENode<2>(4, false, 0.4, 0.5));
        nodes.push_back(new HENode<2>(5, false, 0.6, 0.5));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {2, 5, 4, 1};
        unsigned node_indices_elem_2[3] = {1, 4, 0};
        unsigned node_indices_elem_3[4] = {0, 4, 5, 3};
        for (unsigned i=0; i<3; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }
        nodes_elem_1.push_back(nodes[node_indices_elem_1[3]]);
        nodes_elem_3.push_back(nodes[node_indices_elem_3[3]]);

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Move node 4 so that it overlaps element 0
        ChastePoint<2> point = vertex_mesh.GetNode(4)->GetPoint();
        point.SetCoordinate(1u, 0.7);
        vertex_mesh.GetNode(4)->SetPoint(point);

        // Merge intersection to maintain non-overlapping elements
        vertex_mesh.SetCheckForInternalIntersections(true);
        TS_ASSERT_EQUALS(vertex_mesh.GetCheckForInternalIntersections(), true);
        vertex_mesh.CheckForIntersections();

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.7, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.6, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);

        // Test that each element contains the correct nodes following the rearrangement
        unsigned node_indices_element_0[4] = {2, 3, 4, 5};
        unsigned node_indices_element_1[3] = {2, 5, 1};
        unsigned node_indices_element_2[4] = {1, 5, 4, 0};
        unsigned node_indices_element_3[3] = {0, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.24, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.20, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.36, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.20, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 2.4232, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 2.2806, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 2.7294, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 2.3062, 1e-4);
    }

    void TestPerformIntersectionSwapOtherWayRound()
    {
        /*
         * This test is very similar to TestPerformIntersectionSwap() but with a different ordering
         * of nodes and elements, to ensure full coverage of the CheckForIntersections() method.
         */
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true,  0.0, 0.0));
        nodes.push_back(new HENode<2>(1, true,  1.0, 0.0));
        nodes.push_back(new HENode<2>(2, true,  1.0, 1.0));
        nodes.push_back(new HENode<2>(3, true,  0.0, 1.0));
        nodes.push_back(new HENode<2>(4, false, 0.4, 0.5));
        nodes.push_back(new HENode<2>(5, false, 0.6, 0.5));

        std::vector<HENode<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {2, 3, 4};
        unsigned node_indices_elem_1[4] = {0, 5, 4, 3};
        unsigned node_indices_elem_2[3] = {1, 5, 0};
        unsigned node_indices_elem_3[4] = {2, 4, 5, 1};
        for (unsigned i=0; i<3; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }
        nodes_elem_1.push_back(nodes[node_indices_elem_1[3]]);
        nodes_elem_3.push_back(nodes[node_indices_elem_3[3]]);

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));
        vertex_elements.push_back(new HEElement<2>(2, nodes_elem_2));
        vertex_elements.push_back(new HEElement<2>(3, nodes_elem_3));

        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        // Move node 5 so that it overlaps element 0
        ChastePoint<2> point = vertex_mesh.GetNode(5)->GetPoint();
        point.SetCoordinate(1u, 0.7);
        vertex_mesh.GetNode(5)->SetPoint(point);

        // Merge intersection to maintain non-overlapping elements
        vertex_mesh.SetCheckForInternalIntersections(true);
        vertex_mesh.CheckForIntersections();

        // Test the rearrangement is correct, as in TestPerformIntersectionSwap()
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.6, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.7, 1e-3);

        unsigned node_indices_element_0[4] = {2, 3, 4, 5};
        unsigned node_indices_element_1[3] = {0, 4, 3};
        unsigned node_indices_element_2[4] = {1, 5, 4, 0};
        unsigned node_indices_element_3[3] = {2, 5, 1};
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i), node_indices_element_0[i]);
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(i), node_indices_element_2[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i), node_indices_element_1[i]);
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNodeGlobalIndex(i), node_indices_element_3[i]);
            }
        }

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.24, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.20, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.36, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.20, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 2.4232, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 2.2806, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 2.7294, 1e-4);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 2.3062, 1e-4);
    }
};

#endif /*TestHEMutableVertexMeshReMesh_HPP_*/
