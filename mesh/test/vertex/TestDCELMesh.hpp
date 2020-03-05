#ifndef TESTDCELMESH_HPP_
#define TESTDCELMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "DCELElement.hpp"
#include "DCELHalfEdge.hpp"
#include "DCELVertex.hpp"

#include "FakePetscSetup.hpp"
class TestDCELMesh : public CxxTest::TestSuite
{
public:
    void TestBasicMesh()
    {
        std::vector<DCELVertex<2>* > vertices;
        vertices.push_back(new DCELVertex<2>(0.0,0.0));
        vertices.push_back(new DCELVertex<2>(1.0,0.0));
        vertices.push_back(new DCELVertex<2>(1.0,1.0));
        vertices.push_back(new DCELVertex<2>(0.0,1.0));

        vertices.push_back(new DCELVertex<2>(2.0,0.0));
        vertices.push_back(new DCELVertex<2>(2.0,1.0));

        std::vector<DCELVertex<2>* > element_1_vertices = {vertices[0], vertices[1], vertices[2], vertices[3]};
        std::vector<DCELVertex<2>* > element_2_vertices = {vertices[1], vertices[4], vertices[5], vertices[2]};
        DCELElement<2> element_1(element_1_vertices);
        DCELElement<2> element_2(element_2_vertices);

        DCELHalfEdge<2>* edge_1 = element_1.GetHalfEdge();
        c_vector<double, 2> vert = edge_1->GetTargetVertex()->rGetLocation();
        TS_ASSERT_DELTA(vert(0), 1.0,1e-12);
        TS_ASSERT_DELTA(vert(1), 0.0,1e-12);
        DCELHalfEdge<2>* edge_1_twin = edge_1->GetTwinHalfEdge();
        c_vector<double, 2> vert_twin = edge_1_twin->GetTargetVertex()->rGetLocation();
        TS_ASSERT_DELTA(vert_twin(0), 0.0,1e-12);
        TS_ASSERT_DELTA(vert_twin(1), 0.0,1e-12);
        DCELHalfEdge<2>* walker = edge_1->GetNextHalfEdge();
        while(walker!=edge_1)
        {
            c_vector<double, 2> target_vertex = walker->GetTargetVertex()->rGetLocation();
            DCELHalfEdge<2>* walker_twin = walker->GetTwinHalfEdge();
            walker=walker->GetNextHalfEdge();
        }
        std::cout<<"Backwards: "<<std::endl;
        walker = edge_1->GetPreviousHalfEdge();
        while(walker!=edge_1)
        {
            c_vector<double, 2> target_vertex = walker->GetTargetVertex()->rGetLocation();
            std::cout<<target_vertex(0)<<" "<<target_vertex(1)<<std::endl;
            walker=walker->GetPreviousHalfEdge();
        }

        std::cout<<"Second element"<<std::endl;
        DCELHalfEdge<2>* edge_2 = element_2.GetHalfEdge();
        walker = edge_2->GetNextHalfEdge();
        vert = edge_2->GetTargetVertex()->rGetLocation();
        std::cout<<vert(0)<<" "<<vert(1)<<std::endl;
        while(walker!=edge_2)
        {
            c_vector<double, 2> target_vertex = walker->GetTargetVertex()->rGetLocation();
            std::cout<<target_vertex(0)<<" "<<target_vertex(1)<<std::endl;
            walker=walker->GetNextHalfEdge();
        }
        std::cout<<"Backwards: "<<std::endl;
        walker = edge_2->GetPreviousHalfEdge();
        while(walker!=edge_2)
        {
            c_vector<double, 2> target_vertex = walker->GetTargetVertex()->rGetLocation();
            std::cout<<target_vertex(0)<<" "<<target_vertex(1)<<std::endl;
            walker=walker->GetPreviousHalfEdge();
        }
    }
};
#endif /*TESTDCELMESH_HPP_*/


