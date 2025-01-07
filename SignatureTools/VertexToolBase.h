#ifndef VERTEX_TOOLBASE_H
#define VERTEX_TOOLBASE_H

#include "TVector3.h"
#include "art/Framework/Principal/Event.h"

namespace signature
{

class VertexToolBase 
{
public:
    virtual ~VertexToolBase() = default;
    virtual TVector3 findVertex(const art::Event& evt) const = 0;
};

}

#endif