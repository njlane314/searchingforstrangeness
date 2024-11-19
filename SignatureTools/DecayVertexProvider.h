#ifndef DECAYVERTEXPROVIDER_H
#define DECAYVERTEXPROVIDER_H

#include "TVector3.h"
#include "art/Framework/Principal/Event.h"

class DecayVertexProvider 
{
public:
    virtual ~DecayVertexProvider() = default;
    virtual TVector3 getDecayVertex(const art::Event& evt) const = 0;
};

#endif