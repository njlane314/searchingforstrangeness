#ifndef SIGNATURE_INTEGRITY_H
#define SIGNATURE_INTEGRITY_H

#include "ClarityToolBase.h" 

namespace claritytools {

class SignatureIntegrity : ClarityToolBase {

public:
    explicit SignatureIntegrity(const fhicl::ParameterSet& pset) :
     _chan_act_reg{pset.get<int>("ChannelActiveRegion", 3)}
    {
        configure(pset);
    }

    ~SignatureIntegrity() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    bool filter(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);

private:


   const int _chan_act_reg;

};

bool SignatureIntegrity::filter(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
{
    auto isChannelRegionActive = [&](const TVector3& point) -> bool {
        for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
            try {
                geo::WireID wire = _geo->NearestWireID(point, plane);
                raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

                for (int offset = -_chan_act_reg; offset <= _chan_act_reg; ++offset) {
                    raw::ChannelID_t neighboring_channel = central_channel + offset;

                    if (neighboring_channel < 0 || static_cast<size_t>(neighboring_channel) >= _geo->Nchannels())
                        continue; 

                    if (_bad_channel_mask[neighboring_channel])
                        return false; 
                }
            } catch (const cet::exception&) {
                return false; 
            }
        }
        return true;
    };

    for (const auto& sig : patt) {
        for (const auto& mcp_s : sig) {
            TVector3 start(mcp_s->Vx(), mcp_s->Vy(), mcp_s->Vz());
            if (!isChannelRegionActive(start))
                return false;

            if (std::abs(mcp_s->PdgCode()) != 13) { 
                TVector3 end(mcp_s->EndX(), mcp_s->EndY(), mcp_s->EndZ());
                if (!isChannelRegionActive(end))
                    return false;
            }
        }
    }

    return true;
}


DEFINE_ART_CLASS_TOOL(SignatureIntegrity)

}

#endif
