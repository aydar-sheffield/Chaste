/*
 * HENode.hpp
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef HENODE_HPP_
#define HENODE_HPP_

#include "HalfEdge.hpp"
#include "UblasVectorInclude.hpp"
#include "Node.hpp"
template<unsigned int SPACE_DIM>
class HENode : public Node<SPACE_DIM>
{
public:
    /**
     * Constructor that takes in the node's location as a std::vector.
     *
     * @param index  the index of the node in the mesh
     * @param coords  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param edge the pointer to the outgoing edge
     */
    HENode(unsigned index, std::vector<double> coords,
             bool isBoundaryNode=false, HalfEdge<SPACE_DIM>* pEdge = nullptr);

    /**
     * Constructor that takes in the node's location as a c_vector.
     *
     * @param index  the index of the node in the mesh
     * @param location  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param edge the pointer to the outgoing edge
     */
    HENode(unsigned index, c_vector<double, SPACE_DIM> location,
             bool isBoundaryNode=false, HalfEdge<SPACE_DIM>* pEdge = nullptr);

    /**
     * Constructor that takes the coordinates of the node's location as separate input arguments.
     *
     * @param index  the index of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param v1 the x-coordinate of the node in the mesh (defaults to 0)
     * @param v2 the y-coordinate of the node in the mesh (defaults to 0)
     * @param v3 the z-coordinate of the node in the mesh (defaults to 0)
     * @param edge the pointer to the outgoing edge
     */
    HENode(unsigned index, bool isBoundaryNode=false,
             double v1=0, double v2=0, double v3=0,
             HalfEdge<SPACE_DIM>* pEdge = nullptr);
    /**
     * Constructor to convert Node object to HENode object
     * @param node
     */
    HENode(const Node<SPACE_DIM> &rNode, HalfEdge<SPACE_DIM>* pEdge = nullptr);

    ~HENode();
    HalfEdge<SPACE_DIM>* GetOutgoingEdge() const override;
    void SetOutgoingEdge(HalfEdge<SPACE_DIM>* pEdge) override;

    HalfEdge<SPACE_DIM>* GetIncomingEdge() const;
    /** Get Node the outgoing edge points to **/
    HENode<SPACE_DIM>* GetNextNode() const;

    std::set<unsigned> GetContainingElementIndices() const;

    virtual unsigned int GetNumContainingElements() const override;
    //Forward declaration
    class OutgoingEdgeIterator;

    /**
     * @return an iterator to the first outgoing edge of the node
     *
     * @param skipDeletedElements whether to include deleted edge
     */
    inline OutgoingEdgeIterator GetOutgoingEdgeIteratorBegin(bool skipDeletedEdges = true);

    /**
     * @return an iterator to one past the last Node in the element.
     */
    inline OutgoingEdgeIterator GetOutgoingEdgeIteratorEnd();

    /**
     * A smart iterator over the outgoing edges of the node.
     */
    class OutgoingEdgeIterator
    {
    public:
        /**
         * Dereference the iterator giving you a pointer to the current edge.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline HalfEdge<SPACE_DIM>* operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline HalfEdge<SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename HENode<SPACE_DIM>::OutgoingEdgeIterator& rOther);

        OutgoingEdgeIterator& operator=(const typename HENode<SPACE_DIM>::OutgoingEdgeIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline OutgoingEdgeIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * HENode::GetElementIteratorBegin and HENode::GetElementIteratorEnd instead.
         *
         * @param rElement the element to iterator over
         * @param edgeIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        OutgoingEdgeIterator(HENode<SPACE_DIM>& rNode,
                             bool isLastOutEdge = false,
                             bool skipDeletedEdges = true);

    private:
        /** The node over which we're iterating outgoing edges. */
        HENode<SPACE_DIM>& mrNode;

        /** The actual edge iterator. */
        HalfEdge<SPACE_DIM>* mpEdgeIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedEdges;

        bool mPastFirstIterator;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedEdge();
    };

    //Forward declaration
    class IncomingEdgeIterator;

    /**
     * @return an iterator to the first outgoing edge of the node
     *
     * @param skipDeletedElements whether to include deleted edge
     */
    inline IncomingEdgeIterator GetIncomingEdgeIteratorBegin(bool skipDeletedEdges = true);

    /**
     * @return an iterator to one past the last Node in the element.
     */
    inline IncomingEdgeIterator GetIncomingEdgeIteratorEnd();

    /**
     * A smart iterator over the Incoming edges of the node.
     */
    class IncomingEdgeIterator
    {
    public:
        /**
         * Dereference the iterator giving you a pointer to the current edge.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline HalfEdge<SPACE_DIM>* operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline HalfEdge<SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename HENode<SPACE_DIM>::IncomingEdgeIterator& rOther);

        IncomingEdgeIterator& operator=(const typename HENode<SPACE_DIM>::IncomingEdgeIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline IncomingEdgeIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * HENode::GetElementIteratorBegin and HENode::GetElementIteratorEnd instead.
         *
         * @param rElement the element to iterator over
         * @param edgeIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        IncomingEdgeIterator(HENode<SPACE_DIM>& rNode,
                             bool isLastOutEdge = false,
                             bool skipDeletedEdges = true);

    private:
        /** The node over which we're iterating Incoming edges. */
        HENode<SPACE_DIM>& mrNode;

        /** The actual edge iterator. */
        HalfEdge<SPACE_DIM>* mpEdgeIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedEdges;

        bool mPastFirstIterator;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedEdge();
    };

private:
    /** Outgoing edge **/
    HalfEdge<SPACE_DIM>* mpEdge;
};

//////////////////////////////////////////////////////////////////////////////
// OutgoingEdgeIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::OutgoingEdgeIterator HENode<SPACE_DIM>::GetOutgoingEdgeIteratorBegin(bool skipDeletedEdges)
{
    return OutgoingEdgeIterator(*this, false, skipDeletedEdges);
}

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::OutgoingEdgeIterator HENode<SPACE_DIM>::GetOutgoingEdgeIteratorEnd()
{
    return OutgoingEdgeIterator(*this, true);
}

template <unsigned SPACE_DIM>
HalfEdge<SPACE_DIM>* HENode<SPACE_DIM>::OutgoingEdgeIterator::operator*()
{
    return mpEdgeIter;
}

template <unsigned SPACE_DIM>
HalfEdge<SPACE_DIM>* HENode<SPACE_DIM>::OutgoingEdgeIterator::operator->()
{
    return mpEdgeIter;
}

template <unsigned SPACE_DIM>
bool HENode<SPACE_DIM>::OutgoingEdgeIterator::operator!=(const typename HENode<SPACE_DIM>::OutgoingEdgeIterator& rOther)
{
    return (mpEdgeIter != rOther.mpEdgeIter) || (mPastFirstIterator != rOther.mPastFirstIterator);
}

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::OutgoingEdgeIterator& HENode<SPACE_DIM>::OutgoingEdgeIterator::operator=(const typename HENode<SPACE_DIM>::OutgoingEdgeIterator& rOther)
{
    mpEdgeIter = rOther.mpEdgeIter;
    mrNode = rOther.mrNode;
    mPastFirstIterator = rOther.mPastFirstIterator;
    return *this;
}

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::OutgoingEdgeIterator& HENode<SPACE_DIM>::OutgoingEdgeIterator::operator++()
{
    mPastFirstIterator = true;
    do
    {
        mpEdgeIter = mpEdgeIter->GetTwinHalfEdge()->GetNextHalfEdge();
    }while (!IsAtEnd() && !IsAllowedEdge());
    return (*this);
}

template <unsigned SPACE_DIM>
HENode<SPACE_DIM>::OutgoingEdgeIterator::OutgoingEdgeIterator(HENode<SPACE_DIM>& rNode,
                                                              bool isLastOutEdge,
                                                              bool skipDeletedEdges)
: mrNode(rNode),
  mpEdgeIter(rNode.GetOutgoingEdge()),
  mSkipDeletedEdges(skipDeletedEdges),
  mPastFirstIterator(isLastOutEdge)
{
    assert(mpEdgeIter);
    if (!IsAllowedEdge())
        ++(*this);
}

template <unsigned SPACE_DIM>
bool HENode<SPACE_DIM>::OutgoingEdgeIterator::IsAtEnd()
{
    return mrNode.GetOutgoingEdge()==mpEdgeIter && mPastFirstIterator;
}

template <unsigned int SPACE_DIM>
bool HENode<SPACE_DIM>::OutgoingEdgeIterator::IsAllowedEdge()
{
    return !(mSkipDeletedEdges&&mpEdgeIter->IsDeleted());
}

//////////////////////////////////////////////////////////////////////////////
// IncomingEdgeIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::IncomingEdgeIterator HENode<SPACE_DIM>::GetIncomingEdgeIteratorBegin(bool skipDeletedEdges)
{
    return IncomingEdgeIterator(*this, false, skipDeletedEdges);
}

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::IncomingEdgeIterator HENode<SPACE_DIM>::GetIncomingEdgeIteratorEnd()
{
    return IncomingEdgeIterator(*this, true);
}

template <unsigned SPACE_DIM>
HalfEdge<SPACE_DIM>* HENode<SPACE_DIM>::IncomingEdgeIterator::operator*()
{
    return mpEdgeIter;
}

template <unsigned SPACE_DIM>
HalfEdge<SPACE_DIM>* HENode<SPACE_DIM>::IncomingEdgeIterator::operator->()
{
    return mpEdgeIter;
}

template <unsigned SPACE_DIM>
bool HENode<SPACE_DIM>::IncomingEdgeIterator::operator!=(const typename HENode<SPACE_DIM>::IncomingEdgeIterator& rOther)
{
    return (mpEdgeIter != rOther.mpEdgeIter) || (mPastFirstIterator != rOther.mPastFirstIterator);
}

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::IncomingEdgeIterator& HENode<SPACE_DIM>::IncomingEdgeIterator::operator=(const typename HENode<SPACE_DIM>::IncomingEdgeIterator& rOther)
{
    mpEdgeIter = rOther.mpEdgeIter;
    mrNode = rOther.mrNode;
    mPastFirstIterator = rOther.mPastFirstIterator;
    return *this;
}

template <unsigned SPACE_DIM>
typename HENode<SPACE_DIM>::IncomingEdgeIterator& HENode<SPACE_DIM>::IncomingEdgeIterator::operator++()
{
    mPastFirstIterator = true;
    do
    {
        mpEdgeIter = mpEdgeIter->GetNextHalfEdge()->GetTwinHalfEdge();
    }while (!IsAtEnd() && !IsAllowedEdge());
    return (*this);
}

template <unsigned SPACE_DIM>
HENode<SPACE_DIM>::IncomingEdgeIterator::IncomingEdgeIterator(HENode<SPACE_DIM>& rNode,
                                                              bool isLastOutEdge,
                                                              bool skipDeletedEdges)
: mrNode(rNode),
  mpEdgeIter(rNode.GetOutgoingEdge()->GetTwinHalfEdge()),
  mSkipDeletedEdges(skipDeletedEdges),
  mPastFirstIterator(isLastOutEdge)
{
    assert(mpEdgeIter);
    if (!IsAllowedEdge())
        ++(*this);
}

template <unsigned SPACE_DIM>
bool HENode<SPACE_DIM>::IncomingEdgeIterator::IsAtEnd()
{
    return mrNode.GetOutgoingEdge()->GetTwinHalfEdge()==mpEdgeIter && mPastFirstIterator;
}

template <unsigned int SPACE_DIM>
bool HENode<SPACE_DIM>::IncomingEdgeIterator::IsAllowedEdge()
{
    return !(mSkipDeletedEdges&&mpEdgeIter->IsDeleted());
}

#endif /* HENODE_HPP_ */
