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

    void TestEdgeIterators()
    {
        //Here we create a node with three outgoing (incoming) edges
        HENode<2>* node = new HENode<2>(0, false, 0.0, 0.0);

        //Make edges
        HalfEdge<2>* out_edge_0 = new HalfEdge<2>(); //outgoing edge
        HalfEdge<2>* twin_out_edge_0 = new HalfEdge<2>(); //its twin, also incoming edge
        out_edge_0->SetTwinHalfEdge(twin_out_edge_0,true);
        node->SetOutgoingEdge(out_edge_0);
        twin_out_edge_0->SetTargetNode(node);

        //Another edge outgoing from the node
        HalfEdge<2>* out_edge_1 = new HalfEdge<2>(); //outgoing edge
        HalfEdge<2>* twin_out_edge_1 = new HalfEdge<2>(); //its twin, also incoming edge
        out_edge_1->SetTwinHalfEdge(twin_out_edge_1, true);
        twin_out_edge_1->SetTargetNode(node);

        //Another edge outgoing from the node
        HalfEdge<2>* out_edge_2 = new HalfEdge<2>(); //outgoing edge
        HalfEdge<2>* twin_out_edge_2 = new HalfEdge<2>(); //its twin, also incoming edge
        out_edge_2->SetTwinHalfEdge(twin_out_edge_2, true);
        twin_out_edge_2->SetTargetNode(node);

        twin_out_edge_0->SetNextHalfEdge(out_edge_1,true);
        twin_out_edge_1->SetNextHalfEdge(out_edge_2,true);
        twin_out_edge_2->SetNextHalfEdge(out_edge_0, true);

        //Testing Outgoing edge iterators
        typename HENode<2>::OutgoingEdgeIterator start_iter = node->GetOutgoingEdgeIteratorBegin();
        typename HENode<2>::OutgoingEdgeIterator end_iter = node->GetOutgoingEdgeIteratorEnd();

        TS_ASSERT(start_iter!=end_iter);
        TS_ASSERT(start_iter->GetOriginNode()==node);
        TS_ASSERT(end_iter->GetOriginNode()==node);

        TS_ASSERT(*start_iter==out_edge_0);
        TS_ASSERT(*end_iter==out_edge_0);
        ++start_iter;

        TS_ASSERT(*start_iter==out_edge_1);
        TS_ASSERT(start_iter->GetOriginNode()==node);
        ++start_iter;
        TS_ASSERT(*start_iter==out_edge_2);
        TS_ASSERT(start_iter->GetOriginNode()==node);
        ++start_iter;
        TS_ASSERT(! (start_iter!=end_iter));
        TS_ASSERT(start_iter->GetOriginNode()==node);

        unsigned int n_edges= 0;
        start_iter = node->GetOutgoingEdgeIteratorBegin();
        for (; start_iter!=end_iter; ++start_iter)
        {
            n_edges++;
        }
        TS_ASSERT_EQUALS(n_edges, 3u);

        //Testing Incoming edge iterators
        typename HENode<2>::IncomingEdgeIterator start_iter_in = node->GetIncomingEdgeIteratorBegin();
        typename HENode<2>::IncomingEdgeIterator end_iter_in = node->GetIncomingEdgeIteratorEnd();

        TS_ASSERT(start_iter_in!=end_iter_in);
        TS_ASSERT(start_iter_in->GetTargetNode()==node);
        TS_ASSERT(end_iter_in->GetTargetNode()==node);

        TS_ASSERT(*start_iter_in==twin_out_edge_0);
        TS_ASSERT(*end_iter_in==twin_out_edge_0);
        ++start_iter_in;

        TS_ASSERT(*start_iter_in==twin_out_edge_1);
        TS_ASSERT(start_iter_in->GetTargetNode()==node);
        ++start_iter_in;
        TS_ASSERT(*start_iter_in==twin_out_edge_2);
        TS_ASSERT(start_iter_in->GetTargetNode()==node);
        ++start_iter_in;
        TS_ASSERT(! (start_iter_in!=end_iter_in));
        TS_ASSERT(start_iter_in->GetTargetNode()==node);

        n_edges= 0;
        start_iter_in = node->GetIncomingEdgeIteratorBegin();
        for (; start_iter_in!=end_iter_in; ++start_iter_in)
        {
            n_edges++;
        }
        TS_ASSERT_EQUALS(n_edges, 3u);
    }


};
#endif /*TESTHENODE_HPP_*/


