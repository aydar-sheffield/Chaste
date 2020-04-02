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

#ifndef TESTHEMUTABLEVERTEXMESH_HPP_
#define TESTHEMUTABLEVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <set>

#include "HEMutableVertexMesh.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestMutableVertexMesh : public CxxTest::TestSuite
{
public:

    // This also tests that boundary nodes are updated on element division
    void TestDivideVertexElementGivenAxisOfDivision()
    {
        // Make five nodes, 0, 1 and 2 are boundary nodes
        std::vector<HENode<2>*> nodes;
        nodes.push_back(new HENode<2>(0, true, 1.0, -2.0));
        nodes.push_back(new HENode<2>(1, true, 1.0, 2.0));
        nodes.push_back(new HENode<2>(2, true, -1.0, 2.0));
        nodes.push_back(new HENode<2>(3, false, -1.0, -2.0));
        nodes.push_back(new HENode<2>(4, false, 0.0, 3.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<HENode<2>*> nodes_elem_0;
        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[i]);
        }

        // Make a triangular element out of nodes 1,4,2
        std::vector<HENode<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[2]);

        std::vector<HEElement<2>*> vertex_elements;
        vertex_elements.push_back(new HEElement<2>(0, nodes_elem_0));
        vertex_elements.push_back(new HEElement<2>(1, nodes_elem_1));

        // Make a vertex mesh
        HEMutableVertexMesh<2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

        c_vector<double, 2> axis_of_division;
        axis_of_division(0) = 1.0;
        axis_of_division(1) = 0.0;

        // Divide element 0 along given axis
        unsigned new_element_index = vertex_mesh.DivideElementAlongGivenAxis(vertex_mesh.GetElement(0), axis_of_division, true);

        TS_ASSERT_EQUALS(new_element_index, vertex_mesh.GetNumElements()-1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Now test the position of new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.0, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], -1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.0, 1e-8);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);


        // Now test the nodes in each element
        unsigned expected_node_indices_element_0[4] = {5, 6, 3,0};
        unsigned expected_node_indices_element_1[3] = {1, 4, 2};
        unsigned expected_node_indices_element_2[4] = {6, 5, 1, 2};

        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(i)->GetIndex(), expected_node_indices_element_0[i]);
            if (i < 3)
            {
                TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(i)->GetIndex(), expected_node_indices_element_1[i]);
            }
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(i)->GetIndex(), expected_node_indices_element_2[i]);
        }

        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->GetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);

        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->GetContainingElementIndices(), expected_elements_containing_node_6);
    }

};

#endif /*TESTMUTABLEVERTEXMESH_HPP_*/
