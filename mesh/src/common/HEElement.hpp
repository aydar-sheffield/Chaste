/*
 * HEElement.h
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef HEELEMENT_HPP_
#define HEELEMENT_HPP_
#include "HalfEdge.hpp"
#include "HENode.hpp"
#include "AbstractElement.hpp"
#include "VertexElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template<unsigned int SPACE_DIM>
class HEElement : public AbstractElement<SPACE_DIM, SPACE_DIM>
{
private:
    HalfEdge<SPACE_DIM>* mpHalfEdge;
    unsigned int mNumNodes;
    double mVolume;
    double mSurfaceArea;
    void CommonConstructor(const std::vector<HENode<SPACE_DIM>* > node_list);
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractElement<SPACE_DIM, SPACE_DIM> >(*this);
        archive & mpHalfEdge;
    }

public:

    class NodeIterator;
    /**
     * @return an iterator to the first Node in the element.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline NodeIterator GetNodeIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last Node in the element.
     */
    inline NodeIterator GetNodeIteratorEnd();

    /**
     * Default constructor
     * @param index of the element.
     */
    HEElement(unsigned int index = 0);

    /**
     * Constructs the element given the vector of HE vertices
     * Node list must be given in CCW order
     * @param index of the element
     * @param node_list
     */
    HEElement(unsigned int index, const std::vector<HENode<SPACE_DIM>* > node_list);

    /**
     * Constructs the element given an edge. It is assumed that the edge is
     * contained in a doubly-linked list with fully initalized edges.
     * @param index of the element
     * @param edge the element points to
     */
    HEElement(unsigned int index, HalfEdge<SPACE_DIM>* edge);


    HEElement(unsigned int index, const std::vector<Node<SPACE_DIM>* > node_list);
    HEElement(const VertexElement<SPACE_DIM, SPACE_DIM> &rElement);

    ~HEElement();

    /**
     * Returns the halfedge this element associated with. By default its the edge with
     * local index 0
     * @param local_index
     * @return first halfedge
     */
    HalfEdge<SPACE_DIM>* GetHalfEdge(const unsigned int local_index = 0) const;

    /**
     * Finds halfedge by its target node
     * @param pTarget node target
     * @return HalfEdge whose target node is pTarget
     */
    HalfEdge<SPACE_DIM>* GetHalfEdge(HENode<SPACE_DIM>* pTarget) const;

    void SetHalfEdge(HalfEdge<SPACE_DIM>* pEdge);

    /**
     * Get node with specified index. This random access is an expensive operation
     * @param local_index
     * @return Node with local_index
     */
    virtual HENode<SPACE_DIM>* GetNode(unsigned int local_index) const;

    /**
     * Given the local index of a node owned by this element, return the
     * global index of the node in the mesh.
     *
     * @param localIndex the node's local index in this element
     * @return the global index
     */
    virtual unsigned GetNodeGlobalIndex(unsigned localIndex) const;

    virtual void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode) override
    {
        UpdateNode(rIndex, static_cast<HENode<SPACE_DIM>* >(pNode));
    }
    void UpdateNode(const unsigned &rIndex, HENode<SPACE_DIM>* pNode);

    /**
     * Replace one of the nodes in this element with another.
     *
     * @param pOldNode  pointer to the current node
     * @param pNewNode  pointer to the replacement node
     */
    virtual void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode) override
    {
        ReplaceNode(static_cast<HENode<SPACE_DIM>* >(pOldNode), static_cast<HENode<SPACE_DIM>* >(pNewNode));
    }
    void ReplaceNode(HENode<SPACE_DIM>* pOldNode, HENode<SPACE_DIM>* pNewNode);

    virtual void MarkAsDeleted();

    /**
     * Inherited method, which is not necessary in halfedge based data structures
     */
    virtual void RegisterWithNodes()
    {}

    /**
     * Informs edges of this element that they belong to this element
     */
    void RegisterWithHalfEdges();

    /**
     * Inserts a node into the edge, which splits the edge into two.
     * A new edge is created between pNode and the origin of pEdge.
     * @param pEdge edge, where pNode is inserted
     * @param pNode node to be inserted
     * @return newly created halfedge
     */
    HalfEdge<SPACE_DIM>* AddNode(HalfEdge<SPACE_DIM>* pEdge, HENode<SPACE_DIM>* pNode);

    void DeleteNode(HENode<SPACE_DIM>* pNode);

    /**
     * @return the number of nodes
     */
    virtual unsigned int GetNumNodes() const override;

    /**
     * Returns number of nodes, or computes the number and stores it
     * @param compute if true, computes the number of nodes and stores
     * @return number of nodes
     */
    unsigned int GetNumNodes(const bool compute);

    bool IsElementOnBoundary() const;

    /**
     * Fetch volume of this element. This returns element volume since last update.
     * To compute current volume use ComputeVolume()
     * @return volume of this element
     */
    double GetVolume() const;

    /**
     * Fetch surface area of this element. This returns element surface area since last update.
     * To compute current surface area use ComputeSurfaceArea()
     * @return surface area of this element
     */
    double GetSurfaceArea() const;

    /**
     * Computes, returns, and updates volume of this element
     * @return computed volume of the element
     */
    double ComputeVolume();

    /**
     * Computes, returns, and updates surface area of this element
     * @return computed surface area of the element
     */
    double ComputeSurfaceArea();
    /**
     * A smart iterator over the vertices in the element.
     */
    class NodeIterator
    {
    public:
        /**
         * Dereference the iterator giving you a pointer to the current Node.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline HENode<SPACE_DIM>* operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline HENode<SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename HEElement<SPACE_DIM>::NodeIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline NodeIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * HEElement::GetElementIteratorBegin and HEElement::GetElementIteratorEnd instead.
         *
         * @param rElement the element to iterator over
         * @param nodeIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        NodeIterator(HEElement<SPACE_DIM>& rElement,
                       HENode<SPACE_DIM>* nodeIter,
                       bool skipDeletedElements = true);

    private:
        /** The element we're iterating over. */
        HEElement& mrElement;

        /** The actual node iterator. */
        HENode<SPACE_DIM>* mpNodeIter;

        /**
         * We are iterating over half edges
         */
        HalfEdge<SPACE_DIM>* mpEdge;

        /**
         * Local HalfEdge index
         */
        unsigned int local_index;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedVertices;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedNode();
    };
};

//////////////////////////////////////////////////////////////////////////////
// HENodeIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template <unsigned SPACE_DIM>
typename HEElement<SPACE_DIM>::NodeIterator HEElement<SPACE_DIM>::GetNodeIteratorBegin(
    bool skipDeletedElements)
{
    return NodeIterator(*this, this->GetHalfEdge()->GetOriginNode(), skipDeletedElements);
}

template <unsigned SPACE_DIM>
typename HEElement<SPACE_DIM>::NodeIterator HEElement<SPACE_DIM>::GetNodeIteratorEnd()
{
    return NodeIterator(*this, nullptr);
}

template <unsigned SPACE_DIM>
HENode<SPACE_DIM>* HEElement<SPACE_DIM>::NodeIterator::operator*()
{
    assert(!IsAtEnd());
    return mpNodeIter;
}

template <unsigned SPACE_DIM>
HENode<SPACE_DIM>* HEElement<SPACE_DIM>::NodeIterator::operator->()
{
    assert(!IsAtEnd());
    return mpNodeIter;
}

template <unsigned SPACE_DIM>
bool HEElement<SPACE_DIM>::NodeIterator::operator!=(const typename HEElement<SPACE_DIM>::NodeIterator& rOther)
{
    return mpNodeIter != rOther.mpNodeIter;
}

template <unsigned SPACE_DIM>
typename HEElement<SPACE_DIM>::NodeIterator& HEElement<SPACE_DIM>::NodeIterator::operator++()
{
    do
    {
        local_index++;
        mpEdge = mpEdge->GetNextHalfEdge();
        mpNodeIter = mpEdge->GetOriginNode();
    } while (!IsAtEnd() && !IsAllowedNode());
    if (IsAtEnd())
        mpNodeIter = nullptr;
    return (*this);
}

template <unsigned SPACE_DIM>
HEElement<SPACE_DIM>::NodeIterator::NodeIterator(HEElement<SPACE_DIM>& Element,
                                                     HENode<SPACE_DIM>* nodeIter,
                                                     bool skipDeletedVertices)
        : mrElement(Element),
          mpNodeIter(nodeIter),
          mSkipDeletedVertices(skipDeletedVertices)
{
    if (Element.mNumNodes==0||!nodeIter)
    {
        local_index = UINT_MAX;
    }
    else
    {
        HalfEdge<SPACE_DIM>* edge = Element.GetHalfEdge();
        mpEdge = edge;
        local_index = 0;
        while(mpEdge->GetOriginNode()!=mpNodeIter&&local_index<Element.GetNumNodes())
        {
            local_index++;
            mpEdge = mpEdge->GetNextHalfEdge();
        }
        if (!IsAllowedNode())
            ++(*this);
    }
}

template <unsigned SPACE_DIM>
bool HEElement<SPACE_DIM>::NodeIterator::IsAtEnd()
{
    return local_index>=mrElement.GetNumNodes();
}

template <unsigned int SPACE_DIM>
bool HEElement<SPACE_DIM>::NodeIterator::IsAllowedNode()
{
    return !(mSkipDeletedVertices&&mpEdge->GetOriginNode()->IsDeleted());
}

#endif /* _HEElement_HPP_ */
