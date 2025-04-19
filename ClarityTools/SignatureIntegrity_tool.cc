#ifndef SIGNATURE_INTEGRITY_H
#define SIGNATURE_INTEGRITY_H

#include "ClarityToolBase.h"
#include "../CommonFunctions/Corrections.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/WireReadout.h"

namespace signature 
{
    class SignatureIntegrity : public ClarityToolBase {
    public:
        explicit SignatureIntegrity(const fhicl::ParameterSet& pset)
            : ClarityToolBase(pset),
            _channelActiveRegion{pset.get<int>("ChannelActiveRegion", 3)}
        {
            auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
            _bad_channel_mask.resize(wireReadout.Nchannels(), false);
        }

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
            auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
            for (geo::CryostatID::CryostatID_t c = 0; c < _geoService->Ncryostats(); ++c) {
                geo::CryostatID cryoid(c);
                for (geo::TPCID::TPCID_t t = 0; t < _geoService->NTPC(cryoid); ++t) {
                    geo::TPCID tpcid(cryoid, t);
                    for (geo::PlaneID::PlaneID_t p = 0; p < wireReadout.Nplanes(tpcid); ++p) {
                        geo::PlaneID plane(tpcid, p);
                        if (static_cast<unsigned int>(p) != static_cast<unsigned int>(view)) {
                            continue;
                        }
                        try {
                            geo::Point_t geoPoint(point.X(), point.Y(), point.Z());
                            geo::WireID wire = wireReadout.NearestWireID(geoPoint, plane);
                            raw::ChannelID_t centralChannel = wireReadout.PlaneWireToChannel(wire);
                            for (int offset = -actReg; offset <= actReg; ++offset) {
                                raw::ChannelID_t neighboringChannel = centralChannel + offset;
                                if (neighboringChannel < 0 || static_cast<size_t>(neighboringChannel) >= wireReadout.Nchannels()) {
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
                }
            }
            return true;
        }

        bool checkStart(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const {
            TVector3 mappedStart = {part->Vx(), part->Vy(), part->Vz()};
            float x = static_cast<float>(mappedStart.X());
            float y = static_cast<float>(mappedStart.Y());
            float z = static_cast<float>(mappedStart.Z());
            common::ApplySCEMappingXYZ(x, y, z);
            mappedStart.SetXYZ(x, y, z);
            return isChannelRegionActive(mappedStart, view, _channelActiveRegion);
        }

        bool checkEnd(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const {
            TVector3 mappedEnd = {part->EndX(), part->EndY(), part->EndZ()};
            float x = static_cast<float>(mappedEnd.X());
            float y = static_cast<float>(mappedEnd.Y());
            float z = static_cast<float>(mappedEnd.Z());
            common::ApplySCEMappingXYZ(x, y, z);
            mappedEnd.SetXYZ(x, y, z);
            return isChannelRegionActive(mappedEnd, view, _channelActiveRegion);
        }
    };

    DEFINE_ART_CLASS_TOOL(SignatureIntegrity)
}

#endif