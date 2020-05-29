/*
 * HEVertexMesh.hpp
 *
 *  Created on: 11 Mar 2020
 *      Author: aydar
 */

#ifndef HEVERTEXMESH_HPP_
#define HEVERTEXMESH_HPP_

#include "AbstractMesh.hpp"
#include "HEElement.hpp"
#include "VertexMesh.hpp"
#include "FullEdge.hpp"

template <unsigned int SPACE_DIM>
class HEVertexMesh: public AbstractMesh<SPACE_DIM, SPACE_DIM>
{
protected:
    /** List of elements */
    std::vector<HEElement<SPACE_DIM>* > mElements;

    /** List of edge */
    std::vector<FullEdge<SPACE_DIM>* > mFullEdges;

    /**
     * Constructs full edges for easier edge traversal across the mesh
     */
    void ConstructFullEdges();

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     * @return local index
     */
    virtual unsigned SolveNodeMapping(unsigned index) const override;

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a winding number test, which counts the number of times the
     * polygon associated with the element winds around the given point.
     * The point is outside only when this "winding number" vanishes;
     * otherwise, the point is inside.
     *
     * One must decide whether a point on the polygon's boundary is inside
     * or outside: we adopt the standard convention that a point on a left
     * or bottom edge is inside, and a point on a right or top edge is outside.
     * This way, if two distinct polygons share a common boundary segment,
     * then a point on that segment will be in one polygon or the other, but
     * not both at the same time.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return if the point is included in the element.
     */
    bool ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /**
     * Get the local index of a given element which is the start vertex of the edge
     * of the element that the overlapping point rTestPoint is closest to.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return the local index
     */
    unsigned GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

public:
    /** Forward declaration of element iterator. */
    class HEElementIterator;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline HEElementIterator GetElementIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline HEElementIterator GetElementIteratorEnd();

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to HEVertices
     * @param elements vector of pointers to HEElements
     */
    HEVertexMesh(std::vector<HENode<SPACE_DIM>* > vertices,
                 std::vector<HEElement<SPACE_DIM>* > elements);

    /**
     * Constructor that converts VertexElements into HEVertexMesh
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     */
    HEVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                 std::vector<VertexElement<SPACE_DIM, SPACE_DIM>*> vertexElements);

    /**
     * Constructor that creates HEVertexMesh from VertexMesh
     * @param vertex_mesh mesh to be copied into HEVertexMesh
     */
    HEVertexMesh(VertexMesh<SPACE_DIM, SPACE_DIM>* vertex_mesh);

    /**
     * Default
     */
    HEVertexMesh();

    /**
     * Destructor
     */
    ~HEVertexMesh();

    /**
     * Convert Nodes and VertexElements into HEVertices and HEElements
     * @param nodes
     * @param vertexElements
     */
    void ConvertFromVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                 std::vector<VertexElement<SPACE_DIM, SPACE_DIM>*> vertexElements);

    /**
     * Delete mNodes, mEdges and mElements.
     */
    virtual void Clear();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of full edges
     */
    virtual unsigned int GetNumFullEdges() const;

    /**
     * @return the number of full edges
     */
    unsigned int GetNumAllFullEdges() const;

    /**
     * @return the number of HEElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned int GetNumAllElements() const;

    /**
     * @param index  the global index of a specified halfedge element.
     *
     * @return a pointer to the vertex element
     */
    HEElement<SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Get the node with a given index in the mesh.
     *
     * @param index the global index of the node
     * @return a pointer to the node.
     */
    virtual HENode<SPACE_DIM>* GetNode(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * A formula for the centroid of a plane polygon may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x, centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Get the volume (or area in 2D, or length in 1D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the surfacearea of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    double ComputeVolumeOfElement(unsigned index);
    double ComputeSurfaceAreaOfElement(unsigned index);

    FullEdge<SPACE_DIM>* GetFullEdge(const unsigned int index) const;

    /**
     * A smart iterator over the elements in the mesh.
     *
     * \todo This is the same as in AbstractTetrahedralMesh and PottsMesh - merge? (#1379)
     */
    class HEElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline HEElement<SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline HEElement<SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename HEVertexMesh<SPACE_DIM>::HEElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline HEElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * VertexMesh::GetElementIteratorBegin and VertexMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        HEElementIterator(HEVertexMesh<SPACE_DIM>& rMesh,
                          typename std::vector<HEElement<SPACE_DIM>*>::iterator elementIter,
                          bool skipDeletedElements = true);

    private:
        /** The mesh we're iterating over. */
        HEVertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<HEElement<SPACE_DIM>*>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedElement();
    };
};

//////////////////////////////////////////////////////////////////////////////
// HEElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template <unsigned SPACE_DIM>
typename HEVertexMesh<SPACE_DIM>::HEElementIterator HEVertexMesh<SPACE_DIM>::GetElementIteratorBegin(
    bool skipDeletedElements)
{
    return HEElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template <unsigned SPACE_DIM>
typename HEVertexMesh<SPACE_DIM>::HEElementIterator HEVertexMesh<SPACE_DIM>::GetElementIteratorEnd()
{
    return HEElementIterator(*this, mElements.end());
}

template <unsigned SPACE_DIM>
HEElement<SPACE_DIM>& HEVertexMesh<SPACE_DIM>::HEElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template <unsigned SPACE_DIM>
HEElement<SPACE_DIM>* HEVertexMesh<SPACE_DIM>::HEElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template <unsigned SPACE_DIM>
bool HEVertexMesh<SPACE_DIM>::HEElementIterator::operator!=(const typename HEVertexMesh<SPACE_DIM>::HEElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template <unsigned SPACE_DIM>
typename HEVertexMesh<SPACE_DIM>::HEElementIterator& HEVertexMesh<SPACE_DIM>::HEElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    } while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template <unsigned SPACE_DIM>
HEVertexMesh<SPACE_DIM>::HEElementIterator::HEElementIterator(
    HEVertexMesh<SPACE_DIM>& rMesh,
    typename std::vector<HEElement<SPACE_DIM>*>::iterator elementIter,
    bool skipDeletedElements)
        : mrMesh(rMesh),
          mElementIter(elementIter),
          mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.empty())
    {
        // Cope with empty meshes
        mElementIter = mrMesh.mElements.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mElementIter == mrMesh.mElements.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template <unsigned SPACE_DIM>
bool HEVertexMesh<SPACE_DIM>::HEElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template <unsigned SPACE_DIM>
bool HEVertexMesh<SPACE_DIM>::HEElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /* HEVERTEXMESH_HPP_ */
