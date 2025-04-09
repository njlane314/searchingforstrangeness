#ifndef SIGNATURE_INTEGRITY_H
#define SIGNATURE_INTEGRITY_H

#include "ClarityToolBase.h"
#include "CommonFunctions/Corrections.h"

namespace signature 
{
    class SignatureIntegrity : public ClarityToolBase {
    public:
        explicit SignatureIntegrity(const fhicl::ParameterSet& pset)
            : ClarityToolBase(pset),
            _channelActiveRegion{pset.get<int>("ChannelActiveRegion", 3)},
            _bad_channel_mask(_geoService->Nchannels(), false) {}

        bool filter(const art::Event& e, const Signature& sig,
                    SignatureType type, common::PandoraView view) const override { 
            bool allActive = true;
            _metrics.start_active.clear();
            _metrics.end_active.clear();
            if (type == signature::kPrimaryMuonSignature) { 
                for (const auto& mcp_s : sig) {
                    bool start = checkStart(mcp_s, view);
                    _metrics.start_active.push_back(start);
                    if (!start) allActive = false;
                }
            } else {
                for (const auto& mcp_s : sig) {
                    bool start = checkStart(mcp_s, view);
                    bool end = checkEnd(mcp_s, view);
                    _metrics.start_active.push_back(start);
                    _metrics.end_active.push_back(end);
                    if (!start || !end) allActive = false;
                }
            }
            return allActive;
        }

        std::string getToolName() const override {
            return "SignatureIntegrity";
        }

        std::unique_ptr<ClarityMetrics> getMetrics() const override {
            return std::make_unique<IntegrityMetrics>(_metrics);
        }

    private:
        int _channelActiveRegion;
        mutable IntegrityMetrics _metrics;
        std::vector<bool> _bad_channel_mask; 

        bool isChannelRegionActive(const TVector3& point, common::PandoraView view, int actReg) const {
            for (geo::PlaneID const& plane : _geoService->IteratePlaneIDs()) {
                if (static_cast<unsigned int>(plane.Plane) != static_cast<unsigned int>(view)) {
                    continue;
                }
                try {
                    geo::WireID wire = _geoService->NearestWireID(point, plane);
                    raw::ChannelID_t centralChannel = _geoService->PlaneWireToChannel(wire);
                    for (int offset = -actReg; offset <= actReg; ++offset) {
                        raw::ChannelID_t neighboringChannel = centralChannel + offset;
                        if (neighboringChannel < 0 || static_cast<size_t>(neighboringChannel) >= _geoService->Nchannels()) {
                            continue;
                        }
                        if (_bad_channel_mask[neighboringChannel]) {
                            return false;
                        }
                    }
                } catch (const cet::exception&) {
                    return false;
                }
            }
            return true;
        }

        bool checkStart(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const {
            TVector3 mappedStart = {part->Vx(), part->Vy(), part->Vz()};
            common::ApplySCEMapping(mappedStart);
            return isChannelRegionActive(mappedStart, view, _channelActiveRegion);
        }

        bool checkEnd(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const {
            TVector3 mappedEnd = {part->EndX(), part->EndY(), part->EndZ()};
            common::ApplySCEMapping(mappedEnd);
            return isChannelRegionActive(mappedEnd, view, _channelActiveRegion);
        }
    };

    DEFINE_ART_CLASS_TOOL(SignatureIntegrity)
}

#endif 