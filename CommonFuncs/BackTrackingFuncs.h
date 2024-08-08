#ifndef BACTRACKINGFUNCS_H
#define BACTRACKINGFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"

#include "lardata/Utilities/FindManyInChainP.h"

namespace common
{

void ApplyDetectorOffsets(const float _vtx_t, const float _vtx_x, const float _vtx_y, const float _vtx_z, float &_xtimeoffset, float &_xsceoffset, float &_ysceoffset, float &_zsceoffset)
{
  auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(_vtx_t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(geo::Point_t(_vtx_x, _vtx_y, _vtx_z));
  _xsceoffset = offset.X();
  _ysceoffset = offset.Y();
  _zsceoffset = offset.Z();
}

art::Ptr<simb::MCParticle> getAssocMCParticle(art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &hittruth, const std::vector<art::Ptr<recob::Hit>> &hits, float &purity, float &completeness)
{
    float pfpcharge = 0; // total hit charge from clusters
    float maxcharge = 0; // charge backtracked to best match

    std::unordered_map<int, double> trkide;
    std::unordered_map<int, float> trkq;
    double maxe = -1, tote = 0;
    art::Ptr<simb::MCParticle> maxp_me; 
    
    for (auto h : hits)
    {
        pfpcharge += h->Integral();
        std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth.at(h.key());
        std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth.data(h.key());

        for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p)
        {
            trkide[particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy;                    //store energy per track id
            trkq[particle_vec[i_p]->TrackId()] += h->Integral() * match_vec[i_p]->ideFraction; //store hit integral associated to this hit
            tote += match_vec[i_p]->energy;                                                    //calculate total energy deposited
            if (trkide[particle_vec[i_p]->TrackId()] > maxe)
            { 
                maxe = trkide[particle_vec[i_p]->TrackId()];
                maxp_me = particle_vec[i_p];
                maxcharge = trkq[particle_vec[i_p]->TrackId()];
            }
        } 
    }

    purity = maxcharge / pfpcharge;
    completeness = 0;

    return maxp_me;
}

struct BtPart
{
public:
    BtPart(const int pdg_, const float px_, const float py_, const float pz_, const float e_, const std::vector<unsigned int> &tids_)
        : pdg(pdg_), px(px_), py(py_), pz(pz_), e(e_), tids(tids_) {}
    BtPart(const int pdg_, const float px_, const float py_, const float pz_, const float e_, const unsigned int tid_)
        : pdg(pdg_), px(px_), py(py_), pz(pz_), e(e_) { tids.push_back(tid_); }
    BtPart(const int pdg_,
            const float px_,
            const float py_,
            const float pz_,
            const float e_,
            const std::vector<unsigned int> &tids_,
            const float start_x_,
            const float start_y_,
            const float start_z_,
            const float start_t_) :
        pdg(pdg_),
        px(px_),
        py(py_),
        pz(pz_),
        e(e_),
        tids(tids_),
        start_x(start_x_),
        start_y(start_y_),
        start_z(start_z_),
        start_t(start_t_)
        {}

    BtPart(const int pdg_,
            const float px_,
            const float py_,
            const float pz_,
            const float e_,
            const unsigned int tid_,
            const float start_x_,
            const float start_y_,
            const float start_z_,
            const float start_t_) :
        pdg(pdg_),
        px(px_),
        py(py_),
        pz(pz_),
        e(e_),
        start_x(start_x_),
        start_y(start_y_),
        start_z(start_z_),
        start_t(start_t_)
        { tids.push_back(tid_); }

    int pdg;
    float px, py, pz, e;
    std::vector<unsigned int> tids;
    int nhits = 0;
    float start_x, start_y, start_z, start_t;
};

std::vector<BtPart> initBacktrackingParticleVec(const std::vector<sim::MCShower> &inputMCShower,
                                                const std::vector<sim::MCTrack> &inputMCTrack,
                                                const std::vector<recob::Hit> &inputHits,
                                                const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart)
{
    std::vector<BtPart> btparts_v;
    for (auto mcs : inputMCShower)
    {
        if (mcs.Process() == "primary" || (mcs.MotherPdgCode() == 111 && mcs.Process() == "Decay" && mcs.MotherProcess() == "primary"))
        {
            sim::MCStep mc_step_shower_start = mcs.DetProfile();
            btparts_v.push_back(BtPart(mcs.PdgCode(), mcs.Start().Momentum().Px() * 0.001, mcs.Start().Momentum().Py() * 0.001,
                                        mcs.Start().Momentum().Pz() * 0.001, mcs.Start().Momentum().E() * 0.001, mcs.DaughterTrackID(),
                                        mc_step_shower_start.X(), mc_step_shower_start.Y(), mc_step_shower_start.Z(), mc_step_shower_start.T()));
        }
    }
    for (auto mct : inputMCTrack)
    {
        if (mct.Process() == "primary")
        {
            sim::MCStep mc_step_track_start = mct.Start();
            btparts_v.push_back(BtPart(mct.PdgCode(), mct.Start().Momentum().Px() * 0.001, mct.Start().Momentum().Py() * 0.001,
                                        mct.Start().Momentum().Pz() * 0.001, mct.Start().Momentum().E() * 0.001, mct.TrackID(),
                                        mc_step_track_start.X(), mc_step_track_start.Y(), mc_step_track_start.Z(), mc_step_track_start.T()));
        }
    }

    for (unsigned int ih = 0; ih < inputHits.size(); ih++)
    {
        auto assmcp = assocMCPart->at(ih);
        auto assmdt = assocMCPart->data(ih);
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
        {
            auto mcp = assmcp[ia];
            auto amd = assmdt[ia];
            if (amd->isMaxIDE != 1)
                continue;
            for (auto &btp : btparts_v)
            {
                if (std::find(btp.tids.begin(), btp.tids.end(), mcp->TrackId()) != btp.tids.end())
                {
                    btp.nhits++;
                }
            } 
        }
    }

    return btparts_v;
}

int getAssocBtPart(const std::vector<art::Ptr<recob::Hit>> &hits,
                   const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
                   const std::vector<BtPart> &btpartsv,
                   float &purity,
                   float &completeness,
                   float &overlay_purity)
{
    std::vector<unsigned int> bthitsv(btpartsv.size(), 0);
    
    for (unsigned int ih = 0; ih < hits.size(); ih++)
    {
        art::Ptr<recob::Hit> hitp = hits[ih];
        auto assmcp = assocMCPart->at(hitp.key());
        auto assmdt = assocMCPart->data(hitp.key());
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
        {
            auto mcp = assmcp[ia];
            auto amd = assmdt[ia];
            if (amd->isMaxIDE != 1)
                continue;
            for (unsigned int ib = 0; ib < btpartsv.size(); ++ib)
            {
                auto &btp = btpartsv[ib];
                if (std::find(btp.tids.begin(), btp.tids.end(), mcp->TrackId()) != btp.tids.end())
                {
                    bthitsv[ib]++;
                }
            }
        }
    }

    purity = 0.;
    completeness = 0.;
    unsigned int maxel = (std::max_element(bthitsv.begin(), bthitsv.end()) - bthitsv.begin());
    if (maxel == bthitsv.size())
        return -1;

    if (bthitsv[maxel] == 0)
        return -1;
    
    purity = float(bthitsv[maxel]) / float(hits.size());
    completeness = float(bthitsv[maxel]) / float(btpartsv[maxel].nhits);
    overlay_purity = 1. - std::accumulate(bthitsv.begin(), bthitsv.end(), 0.) / float(hits.size());
    
    return maxel;
}

int getAssocBtPart(const std::vector<art::Ptr<recob::Hit>> &hits,
                   const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
                   const std::vector<BtPart> &btpartsv,
                   float &purity,
                   float &completeness)
{
    std::vector<unsigned int> bthitsv(btpartsv.size(), 0);

    for (unsigned int ih = 0; ih < hits.size(); ih++)
    {
        art::Ptr<recob::Hit> hitp = hits[ih];
        auto assmcp = assocMCPart->at(hitp.key());
        auto assmdt = assocMCPart->data(hitp.key());
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
        {
            auto mcp = assmcp[ia];
            auto amd = assmdt[ia];
            if (amd->isMaxIDE != 1)
                continue;

            for (unsigned int ib = 0; ib < btpartsv.size(); ++ib)
            {
                auto &btp = btpartsv[ib];
                if (std::find(btp.tids.begin(), btp.tids.end(), mcp->TrackId()) != btp.tids.end())
                {
                    bthitsv[ib]++;
                }
            }
        }
    }

    purity = 0.;
    completeness = 0.;
    unsigned int maxel = (std::max_element(bthitsv.begin(), bthitsv.end()) - bthitsv.begin());

    if (maxel == bthitsv.size())
        return -1;

    if (bthitsv[maxel] == 0)
        return -1;
  
    purity = float(bthitsv[maxel]) / float(hits.size());
    completeness = float(bthitsv[maxel]) / float(btpartsv[maxel].nhits);
  
    return maxel;
}

bool isHitBtMonteCarlo(const size_t hit_index,
                       const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
                       float en_threshold)
{
    std::vector<art::Ptr<simb::MCParticle>> particle_vec = assocMCPart->at(hit_index);
    std::vector<anab::BackTrackerHitMatchingData const *> match_vec = assocMCPart->data(hit_index);
 
    bool found_mc_hit = false;
    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p)
    {
        if (match_vec[i_p]->energy > en_threshold)
        {
            found_mc_hit = true;
            break;
        }
    } 

    return found_mc_hit;
}

} 

#endif