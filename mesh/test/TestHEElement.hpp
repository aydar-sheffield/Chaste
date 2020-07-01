#ifndef TESTHEELEMENT_HPP_
#define TESTHEELEMENT_HPP_

#include <cxxtest/TestSuite.h>
#include "HEElement.hpp"

#include "FakePetscSetup.hpp"
class TestHEElement : public CxxTest::TestSuite
{
public:

    void TestBasicElement()
    {
        /**
         * 3_______2
         * |       |
         * |       |
         * |       |
         * |_______|
         * 0        1
         */
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        std::vector<HENode<2>* > element_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        HEElement<2> element(0, element_vertices);

        //Correct number of nodes
        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);

        //Testing (inner) edge traversal in CCW order
        HalfEdge<2>* edge = element.GetHalfEdge();
        HalfEdge<2>* twin_edge = edge->GetTwinHalfEdge();
        HalfEdge<2>* next_edge = edge->GetNextHalfEdge();
        HalfEdge<2>* next_edge_twin = next_edge->GetTwinHalfEdge();
        HalfEdge<2>* previous_edge = edge->GetPreviousHalfEdge();
        HalfEdge<2>* previous_edge_twin = previous_edge->GetTwinHalfEdge();

        c_vector<double, 2> edge_target = edge->GetTargetNode()->rGetLocation();
        c_vector<double, 2> next_edge_twin_target = next_edge_twin->GetTargetNode()->rGetLocation();
        TS_ASSERT_DELTA(edge_target(0), 1.0,1e-12);
        TS_ASSERT_DELTA(edge_target(1), 0.0,1e-12);
        TS_ASSERT_DELTA(edge_target(0), next_edge_twin_target(0),1e-12);
        TS_ASSERT_DELTA(edge_target(1), next_edge_twin_target(1),1e-12);

        c_vector<double, 2> edge_twin_target = twin_edge->GetTargetNode()->rGetLocation();
        c_vector<double, 2> previous_edge_target = previous_edge->GetTargetNode()->rGetLocation();
        TS_ASSERT_DELTA(edge_twin_target(0), previous_edge_target(0), 1e-12);
        TS_ASSERT_DELTA(edge_twin_target(1), previous_edge_target(1), 1e-12);

        //Should be vertex 2
        c_vector<double, 2> previous_edge_twin_target = previous_edge_twin->GetTargetNode()->rGetLocation();
        TS_ASSERT_DELTA(previous_edge_twin_target(0), 0.0,1e-12);
        TS_ASSERT_DELTA(previous_edge_twin_target(1), 1.0,1e-12);
        //Testing edge target traversal in CCW order
        while(next_edge!=edge)
        {
            TS_ASSERT(next_edge->GetElement()== &element);
            previous_edge = next_edge->GetPreviousHalfEdge();
            edge_target = previous_edge->GetTargetNode()->rGetLocation();
            next_edge_twin = next_edge->GetTwinHalfEdge();
            next_edge_twin_target = next_edge_twin->GetTargetNode()->rGetLocation();
            TS_ASSERT_DELTA(edge_target(0), next_edge_twin_target(0),1e-12);
            TS_ASSERT_DELTA(edge_target(1), next_edge_twin_target(1),1e-12);
            next_edge = next_edge->GetNextHalfEdge();
        }

        //Testing (outer) edge traversal in CW order
        TS_ASSERT(!twin_edge->GetElement());
        HalfEdge<2>* next_twin_edge = twin_edge->GetNextHalfEdge();
        while (next_twin_edge!=twin_edge)
        {
            TS_ASSERT(!next_twin_edge->GetElement());
            HalfEdge<2>* next_twin_edge_twin = next_twin_edge->GetTwinHalfEdge();
            HalfEdge<2>* previous_twin_edge = next_twin_edge->GetPreviousHalfEdge();
            c_vector<double, 2> next_twin_edge_twin_target = next_twin_edge_twin->GetTargetNode()->rGetLocation();
            c_vector<double, 2> previous_twin_edge_target = previous_twin_edge->GetTargetNode()->rGetLocation();
            TS_ASSERT_DELTA(next_twin_edge_twin_target(0), previous_twin_edge_target(0),1e-12);
            TS_ASSERT_DELTA(next_twin_edge_twin_target(1), previous_twin_edge_target(1),1e-12);
            next_twin_edge = next_twin_edge->GetNextHalfEdge();
        }

    }

    void TestHEElementNodeIterator()
    {
        /**
         * 3_______2
         * |       |
         * |       |
         * |       |
         * |_______|
         * 0        1
         */
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        std::vector<HENode<2>* > element_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        HEElement<2> element(0, element_vertices);

        //Correct number of nodes
        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);

        //Checking if default node iteration works
        unsigned int n_nodes= 0;
        typename HEElement<2>::NodeIterator start_iter = element.GetNodeIteratorBegin();
        typename HEElement<2>::NodeIterator end_iter = element.GetNodeIteratorEnd();
        TS_ASSERT(start_iter!=end_iter);
        //Last iterator points to the first node due to cyclicality
        TS_ASSERT(*end_iter == vertices[1]);
        for (; start_iter != end_iter; ++start_iter, ++n_nodes)
        {
            TS_ASSERT(*start_iter==element_vertices[n_nodes]);
        }
        TS_ASSERT_EQUALS(n_nodes, 4u);
        TS_ASSERT(!(start_iter !=end_iter));
        TS_ASSERT(*start_iter== vertices[1]);

        //Mark nodes 1 and 2 as deleted
        vertices[1]->MarkAsDeleted();
        vertices[2]->MarkAsDeleted();
        start_iter = element.GetNodeIteratorBegin();
        end_iter = element.GetNodeIteratorEnd();
        n_nodes= 0;
        TS_ASSERT(*start_iter == element_vertices[0]);
        ++start_iter;
        TS_ASSERT(*start_iter == element_vertices[3]);
        ++start_iter;
        TS_ASSERT(!(start_iter != end_iter));
    }

    void TestHEElementEdgeIterator()
    {
        /**
         * 3_______2
         * |       |
         * |       |
         * |       |
         * |_______|
         * 0        1
         */
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        std::vector<HENode<2>* > element_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        HEElement<2> element(0, element_vertices);

        //Correct number of nodes
        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);

        //Checking if default node iteration works
        unsigned int n_edges= 0;
        typename HEElement<2>::EdgeIterator start_iter = element.GetEdgeIteratorBegin();
        typename HEElement<2>::EdgeIterator end_iter = element.GetEdgeIteratorEnd();
        TS_ASSERT(start_iter!=end_iter);

        for (; start_iter!=end_iter; ++start_iter, n_edges++)
        {
            TS_ASSERT(*start_iter==element.GetHalfEdge(n_edges));
        }
        TS_ASSERT_EQUALS(n_edges, 4u);
        TS_ASSERT(!(start_iter !=end_iter));
        TS_ASSERT(*start_iter == element.GetHalfEdge());

        n_edges= 0;
        //Mark edges 1 and 2 as deleted
        element.GetHalfEdge(1)->SetDeletedStatus(true);
        element.GetHalfEdge(2)->SetDeletedStatus(true);
        start_iter = element.GetEdgeIteratorBegin();
        end_iter = element.GetEdgeIteratorEnd();

        TS_ASSERT(*start_iter == element.GetHalfEdge());
        ++start_iter;
        TS_ASSERT(*start_iter == element.GetHalfEdge(3));
        ++start_iter;
        TS_ASSERT(!(start_iter != end_iter));
    }

    void TestAddNode()
    {
        /**
         * 3_______2
         * |       |
         * |       |
         * |       |
         * |_______|
         * 0        1
         */
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        std::vector<HENode<2>* > element_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        HEElement<2> element(0, element_vertices);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);

        //Add node between nodes 0 and 1 and below the edge connecting the two nodes
        HENode<2>* added_node = new HENode<2>(4, false, 0.5, -0.5);
        std::vector<HENode<2>* > new_element_vertices = {vertices[0], added_node, vertices[1], vertices[2], vertices[3]};
        HalfEdge<2>* node_inserted_edge = element.GetHalfEdge(vertices[1]);
        element.AddNode(node_inserted_edge, added_node);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 5u);

        unsigned int n_nodes= 0;
        typename HEElement<2>::NodeIterator start_iter = element.GetNodeIteratorBegin();
        typename HEElement<2>::NodeIterator end_iter = element.GetNodeIteratorEnd();
        for (; start_iter != end_iter; ++start_iter, ++n_nodes)
        {
            TS_ASSERT(*start_iter==new_element_vertices[n_nodes]);
        }

        n_nodes= 0;
        HalfEdge<2>* next_edge = element.GetHalfEdge();
        do
        {
            TS_ASSERT(next_edge->GetOriginNode() == new_element_vertices[n_nodes]);
            n_nodes++;
            next_edge = next_edge->GetNextHalfEdge();
        }while(next_edge != element.GetHalfEdge());
    }

    void TestDeleteNodeContainedInOneElement()
    {
        /**
         * 3_______2
         * |       |
         * |       |
         * |       |
         * |_______|
         * 0        1
         */
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        std::vector<HENode<2>* > element_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        HEElement<2> element(0, element_vertices);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);

        //Delete node 1
        element.DeleteNode(vertices[1]);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 3u);
        TS_ASSERT(element.GetHalfEdge()->GetOriginNode()==vertices[0]);
        TS_ASSERT(element.GetHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[2]);
        TS_ASSERT(element.GetHalfEdge()->GetPreviousHalfEdge()->GetOriginNode()==vertices[3]);
    }

    void TestDeleteNodeContainedInTwoElements()
    {
        //Two square elements sharing an edge
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        vertices.push_back(new HENode<2>(4,false,2.0,0.0));
        vertices.push_back(new HENode<2>(5,false,2.0,1.0));

        std::vector<HENode<2>* > element_0_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        std::vector<HENode<2>* > element_1_vertices = {vertices[1], vertices[4], vertices[5], vertices[2]};

        std::vector<HEElement<2>*> elements_2d;
        elements_2d.push_back(new HEElement<2>(0, element_0_vertices));
        elements_2d.push_back(new HEElement<2>(1, element_1_vertices));

        //Delete common node 1 and check if the elements have updated their nodes accordingly
        elements_2d[0]->DeleteNode(vertices[1]);
        TS_ASSERT_EQUALS(elements_2d[0]->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(elements_2d[1]->GetNumNodes(), 3u);

        //Checking if nodes have been properly updated in the first element...
        TS_ASSERT(elements_2d[0]->GetHalfEdge()->GetOriginNode()==vertices[0]);
        TS_ASSERT(elements_2d[0]->GetHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[2]);
        TS_ASSERT(elements_2d[0]->GetHalfEdge()->GetPreviousHalfEdge()->GetOriginNode()==vertices[3]);

        //Checking if nodes have been properly updated in the second element...
        TS_ASSERT(elements_2d[1]->GetHalfEdge()->GetOriginNode()==vertices[2]);
        TS_ASSERT(elements_2d[1]->GetHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[4]);
        TS_ASSERT(elements_2d[1]->GetHalfEdge()->GetPreviousHalfEdge()->GetOriginNode()==vertices[5]);
    }

    void TestDeleteNodeContainedInThreeElements()
    {
        //Two square elements sharing an edge and a rectangular element below, sharing common edges with both
        //A void between elements forms, such that the elements contain no common edge
        std::vector<HENode<2>* > vertices;
        vertices.push_back(new HENode<2>(0,false,0.0,0.0));
        vertices.push_back(new HENode<2>(1,false,1.0,0.0));
        vertices.push_back(new HENode<2>(2,false,1.0,1.0));
        vertices.push_back(new HENode<2>(3,false,0.0,1.0));
        vertices.push_back(new HENode<2>(4,false,2.0,0.0));
        vertices.push_back(new HENode<2>(5,false,2.0,1.0));

        vertices.push_back(new HENode<2>(6,false,0.0,-1.0));
        vertices.push_back(new HENode<2>(7,false,2.0,-1.0));

        std::vector<HENode<2>* > element_0_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        std::vector<HENode<2>* > element_1_vertices = {vertices[1], vertices[4], vertices[5], vertices[2]};
        std::vector<HENode<2>* > element_2_vertices = {vertices[0], vertices[6], vertices[7], vertices[4], vertices[1]};

        std::vector<HEElement<2>*> elements_2d;
        elements_2d.push_back(new HEElement<2>(0, element_0_vertices));
        elements_2d.push_back(new HEElement<2>(1, element_1_vertices));
        elements_2d.push_back(new HEElement<2>(2, element_2_vertices));

        //Delete common node 1 and check if the elements have updated their nodes accordingly
        elements_2d[0]->DeleteNode(vertices[1]);
        TS_ASSERT_EQUALS(elements_2d[0]->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(elements_2d[1]->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(elements_2d[2]->GetNumNodes(), 4u);

        //Checking if nodes have been properly updated in the first element...
        TS_ASSERT(elements_2d[0]->GetHalfEdge()->GetOriginNode()==vertices[0]);
        TS_ASSERT(elements_2d[0]->GetHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[2]);
        TS_ASSERT(elements_2d[0]->GetHalfEdge()->GetPreviousHalfEdge()->GetOriginNode()==vertices[3]);

        //Checking if nodes have been properly updated in the second element...
        TS_ASSERT(elements_2d[1]->GetHalfEdge()->GetOriginNode()==vertices[2]);
        TS_ASSERT(elements_2d[1]->GetHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[4]);
        TS_ASSERT(elements_2d[1]->GetHalfEdge()->GetPreviousHalfEdge()->GetOriginNode()==vertices[5]);

        //Checking if nodes have been properly updated in the third element...
        TS_ASSERT(elements_2d[2]->GetHalfEdge()->GetOriginNode()==vertices[0]);
        TS_ASSERT(elements_2d[2]->GetHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[6]);
        TS_ASSERT(elements_2d[2]->GetHalfEdge()->GetNextHalfEdge()->GetNextHalfEdge()->GetOriginNode()==vertices[7]);
        TS_ASSERT(elements_2d[2]->GetHalfEdge()->GetPreviousHalfEdge()->GetOriginNode()==vertices[4]);
    }
};
#endif /*TESTHEELEMENT_HPP_*/


