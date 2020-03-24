#ifndef TESTHENODE_HPP_
#define TESTHENODE_HPP_

#include <cxxtest/TestSuite.h>
#include "HENode.hpp"

#include "FakePetscSetup.hpp"
class TestHENode : public CxxTest::TestSuite
{
public:

    void TestBasicNode()
    {
        //Testing constructors
        HENode<2> node0(0, false, 3.0,3.0);
        c_vector<double, 2> node_location = node0.rGetLocation();
        TS_ASSERT_DELTA(node_location(0),3.0,1e-12);
        TS_ASSERT_DELTA(node_location(1),3.0,1e-12);
        c_vector<double, 2> location1;
        location1(0) = 1.0;
        location1(1) = 1.0;
        HENode<2> node1(1,location1);
        c_vector<double, 2> node_location1 = node1.rGetLocation();
        TS_ASSERT_DELTA(node_location1(0), location1(0),1e-12);
        TS_ASSERT_DELTA(node_location1(1), location1(1),1e-12);
    }
};
#endif /*TESTHENODE_HPP_*/


