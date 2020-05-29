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
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 4u);
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
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Perform a T2 swap on the central triangle element
        HEElement<2>* p_element_0 = vertex_mesh.GetElement(0);
        vertex_mesh.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 3u);
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
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumNodes(), 5u);

        // Perform a T2 swap on the middle triangle element
        p_element_0 = vertex_mesh2.GetElement(0);
        vertex_mesh2.PerformT2Swap(*p_element_0);

        TS_ASSERT_EQUALS(vertex_mesh2.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh2.GetNumNodes(), 3u);

        TS_ASSERT_EQUALS(vertex_mesh2.GetNumAllElements(), 2u);
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
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumAllElements(), 5u);
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
};

#endif /*TestHEMutableVertexMeshReMesh_HPP_*/
