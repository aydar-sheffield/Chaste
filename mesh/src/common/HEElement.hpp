/*
 * HEElement.h
 *
 *  Created on: 3 Mar 2020
 *      Author: aydar
 */

#ifndef HEELEMENT_HPP_
#define HEELEMENT_HPP_
#include "HalfEdge.hpp"
#include "HEVertex.hpp"
#include "AbstractElement.hpp"
#include "VertexElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template<unsigned int SPACE_DIM>
class HEElement : public AbstractElement<SPACE_DIM, SPACE_DIM>
{
private:
    HalfEdge<SPACE_DIM>* mHalfEdge;
    unsigned int mNumNodes;
    void CommonConstructor(const std::vector<HEVertex<SPACE_DIM>* > vertex_list);
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
        archive & mHalfEdge;
    }

public:
    HEElement(unsigned int index = 0);
    /**
     * Constructs the element given the vector of HE vertices
     * Vertex list must be given in CCW order
     * @param index of the element
     * @param vertex_list
     */
    HEElement(unsigned int index, const std::vector<HEVertex<SPACE_DIM>* > vertex_list);

    /**
     * Constructs the element given the vector of halfedges
     * Edge list must be given in CCW order
     * @param index
     * @param edge_list
     */
    HEElement(unsigned int index, const std::vector<HalfEdge<SPACE_DIM>* > edge_list);

    HEElement(unsigned int index, const std::vector<Node<SPACE_DIM>* > node_list);
    HEElement(const VertexElement<SPACE_DIM, SPACE_DIM> &element);

    ~HEElement();
    HalfEdge<SPACE_DIM>* GetHalfEdge() const;
    void SetHalfEdge(HalfEdge<SPACE_DIM>* edge);

    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode) override
    {
        UpdateNode(rIndex, static_cast<HEVertex<SPACE_DIM>* >(pNode));
    }
    void UpdateNode(const unsigned &rIndex, HEVertex<SPACE_DIM>* pVertex);
    void MarkAsDeleted();
    void RegisterWithNodes();
    void AddNode(const unsigned int &rIndex, HEVertex<SPACE_DIM>* pNode);

    void UpdateNumNodes();
    virtual unsigned int GetNumNodes() const override;
};

#endif /* _HEElement_HPP_ */
