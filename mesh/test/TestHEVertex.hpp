#ifndef TESTHEVERTEX_HPP_
#define TESTHEVERTEX_HPP_

#include <cxxtest/TestSuite.h>
#include "HEVertex.hpp"

#include "FakePetscSetup.hpp"
class TestHEVertex : public CxxTest::TestSuite
{
public:

    void TestBasicVertex()
    {
        //Testing constructors
        HEVertex<2> vertex0(3.0,3.0);
        c_vector<double, 2> vertex_location = vertex0.rGetLocation();
        TS_ASSERT_DELTA(vertex_location(0),3.0,1e-12);
        TS_ASSERT_DELTA(vertex_location(1),3.0,1e-12);
        c_vector<double, 2> location1;
        location1(0) = 1.0;
        location1(1) = 1.0;
        HEVertex<2> vertex1(location1);
        c_vector<double, 2> vertex_location1 = vertex1.rGetLocation();
        TS_ASSERT_DELTA(vertex_location1(0), location1(0),1e-12);
        TS_ASSERT_DELTA(vertex_location1(1), location1(1),1e-12);
    }
};
#endif /*TESTHEVERTEX_HPP_*/


