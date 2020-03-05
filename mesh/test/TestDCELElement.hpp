#ifndef TESTDCELELEMENT_HPP_
#define TESTDCELELEMENT_HPP_

#include <cxxtest/TestSuite.h>
#include "DCELElement.hpp"

#include "FakePetscSetup.hpp"
class TestDCELElement : public CxxTest::TestSuite
{
public:

    void TestBasicElement()
    {
        std::vector<DCELVertex<2>* > vertices;
        vertices.push_back(new DCELVertex<2>(0.0,0.0));
        vertices.push_back(new DCELVertex<2>(1.0,0.0));
        vertices.push_back(new DCELVertex<2>(1.0,1.0));
        vertices.push_back(new DCELVertex<2>(0.0,1.0));
        std::vector<DCELVertex<2>* > element_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        DCELElement<2> element(element_vertices);

        //Testing edge traversal in CCW order
        DCELHalfEdge<2>* edge = element.GetHalfEdge();
        DCELHalfEdge<2>* twin_edge = edge->GetTwinHalfEdge();
        DCELHalfEdge<2>* next_edge = edge->GetNextHalfEdge();
        DCELHalfEdge<2>* next_edge_twin = next_edge->GetTwinHalfEdge();
        DCELHalfEdge<2>* previous_edge = edge->GetPreviousHalfEdge();
        DCELHalfEdge<2>* previous_edge_twin = previous_edge->GetTwinHalfEdge();

        c_vector<double, 2> edge_target = edge->GetTargetVertex()->rGetLocation();
        c_vector<double, 2> next_edge_twin_target = next_edge_twin->GetTargetVertex()->rGetLocation();
        TS_ASSERT_DELTA(edge_target(0), next_edge_twin_target(0),1e-12);
        TS_ASSERT_DELTA(edge_target(1), next_edge_twin_target(1),1e-12);
        c_vector<double, 2> edge_twin_target = twin_edge->GetTargetVertex()->rGetLocation();
        c_vector<double, 2> previous_edge_target = previous_edge->GetTargetVertex()->rGetLocation();
        TS_ASSERT_DELTA(edge_twin_target(0), previous_edge_target(0), 1e-12);
        TS_ASSERT_DELTA(edge_twin_target(1), previous_edge_target(1), 1e-12);
        c_vector<double, 2> previous_edge_twin_target = previous_edge_twin->GetTargetVertex()->rGetLocation();
        TS_ASSERT_DELTA(previous_edge_twin_target(0), 0,1e-12);
        TS_ASSERT_DELTA(previous_edge_twin_target(1), 1.0,1e-12);
        while(next_edge!=edge)
        {
            TS_ASSERT(next_edge->GetElement()== &element);
            previous_edge = next_edge->GetPreviousHalfEdge();
            edge_target = previous_edge->GetTargetVertex()->rGetLocation();
            next_edge_twin = next_edge->GetTwinHalfEdge();
            next_edge_twin_target = next_edge_twin->GetTargetVertex()->rGetLocation();
            TS_ASSERT_DELTA(edge_target(0), next_edge_twin_target(0),1e-12);
            TS_ASSERT_DELTA(edge_target(1), next_edge_twin_target(1),1e-12);
            next_edge = next_edge->GetNextHalfEdge();
        }

    }
};
#endif /*TESTDCELELEMENT_HPP_*/


