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
};
#endif /*TESTHEELEMENT_HPP_*/


