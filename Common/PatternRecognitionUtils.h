#ifndef COMMON_PATTERNRECOGNITIONUTILS_H
#define COMMON_PATTERNRECOGNITIONUTILS_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "Common/BacktrackingUtilities.h"

#include <array>
#include <memory>
#include <vector>

namespace common {

struct PRMetrics {
    bool valid = false;
    std::array<int, 3> nshared = {0, 0, 0};
    std::array<int, 3> nhits_O = {0, 0, 0};
    std::array<int, 3> nhits_T = {0, 0, 0};
    std::array<float, 3> purity = {0.f, 0.f, 0.f};
    std::array<float, 3> completeness = {0.f, 0.f, 0.f};
};

inline size_t CountTruthHitsInSlice(
    int trackid, const std::vector<art::Ptr<recob::Hit>> &hits,
    std::unique_ptr<art::FindManyP<simb::MCParticle,
                                   anab::BackTrackerHitMatchingData>>
        &assocMCPart) {
    if (trackid < 0)
        return 0;
    size_t n = 0;
    for (auto const &h : hits) {
        auto parts = assocMCPart->at(h.key());
        for (auto const &mcp : parts) {
            if (mcp->TrackId() == trackid) {
                ++n;
                break;
            }
        }
    }
    return n;
}

inline PRMetrics ComputePRMetrics(
    const art::Event &event, const art::InputTag &cls_tag,
    const art::InputTag &hit_tag, const art::InputTag &bkt_tag,
    std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
    const std::array<int, 3> &tids) {
    PRMetrics out;

    auto const &cluster_h =
        event.getValidHandle<std::vector<recob::Cluster>>(cls_tag);
    art::FindManyP<recob::Hit> assocHits(cluster_h, event, cls_tag);

    std::vector<std::vector<art::Ptr<recob::Hit>>> pfp_hits(slice_pfp_vec.size());
    std::vector<art::Ptr<recob::Hit>> inputHits;
    inputHits.reserve(1024);

    for (size_t ip = 0; ip < slice_pfp_vec.size(); ++ip) {
        auto clusters = slice_pfp_vec[ip].get<recob::Cluster>();
        for (auto const &c : clusters) {
            auto hits = assocHits.at(c.key());
            pfp_hits[ip].insert(pfp_hits[ip].end(), hits.begin(), hits.end());
            inputHits.insert(inputHits.end(), hits.begin(), hits.end());
        }
    }

    auto const &hit_h =
        event.getValidHandle<std::vector<recob::Hit>>(hit_tag);
    auto assocMCPart = std::make_unique<
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
        hit_h, event, bkt_tag);

    out.nhits_T[0] = static_cast<int>(
        CountTruthHitsInSlice(tids[0], inputHits, assocMCPart));
    out.nhits_T[1] = static_cast<int>(
        CountTruthHitsInSlice(tids[1], inputHits, assocMCPart));
    out.nhits_T[2] = static_cast<int>(
        CountTruthHitsInSlice(tids[2], inputHits, assocMCPart));

    const int n_pfp = static_cast<int>(pfp_hits.size());
    if (n_pfp < 3)
        return out;

    std::vector<std::array<int, 3>> shared(n_pfp, {0, 0, 0});
    std::vector<int> n_hits_O(n_pfp, 0);

    for (int i = 0; i < n_pfp; ++i) {
        n_hits_O[i] = static_cast<int>(pfp_hits[i].size());
        for (auto const &h : pfp_hits[i]) {
            auto parts = assocMCPart->at(h.key());
            bool has[3] = {false, false, false};
            for (auto const &mcp : parts) {
                const int tid = mcp->TrackId();
                for (int t = 0; t < 3; ++t) {
                    if (!has[t] && tid == tids[t])
                        has[t] = true;
                }
            }
            for (int t = 0; t < 3; ++t) {
                if (has[t])
                    ++shared[i][t];
            }
        }
    }

    int best_i = -1, best_j = -1, best_k = -1, best_sum = -1;
    for (int i = 0; i < n_pfp; ++i) {
        for (int j = 0; j < n_pfp; ++j) {
            if (j == i)
                continue;
            for (int k = 0; k < n_pfp; ++k) {
                if (k == i || k == j)
                    continue;
                const int sum = shared[i][0] + shared[j][1] + shared[k][2];
                if (sum > best_sum) {
                    best_sum = sum;
                    best_i = i;
                    best_j = j;
                    best_k = k;
                }
            }
        }
    }
    if (best_sum < 0)
        return out;

    out.valid = true;

    auto fill_pair = [&](int i, int t) {
        out.nshared[t] = shared[i][t];
        out.nhits_O[t] = n_hits_O[i];
        out.purity[t] =
            (out.nhits_O[t] > 0)
                ? static_cast<float>(out.nshared[t]) / out.nhits_O[t]
                : 0.f;
        out.completeness[t] =
            (out.nhits_T[t] > 0)
                ? static_cast<float>(out.nshared[t]) / out.nhits_T[t]
                : 0.f;
    };

    fill_pair(best_i, 0);
    fill_pair(best_j, 1);
    fill_pair(best_k, 2);

    return out;
}

} // namespace common

#endif
