#ifndef COMMON_PATTERNRECOGNITIONUTILS_H
#define COMMON_PATTERNRECOGNITIONUTILS_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "Common/BacktrackingUtilities.h"

#include <array>
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

namespace common {

struct PRMetrics {
    bool valid = false;
    std::vector<int> nshared;
    std::vector<int> nhits_O;
    std::vector<int> nhits_T;
    std::vector<float> purity;
    std::vector<float> completeness;
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
    const std::vector<int> &tids) {
    PRMetrics out;
    const size_t n_truth = tids.size();
    if (n_truth == 0)
        return out;

    out.nshared.assign(n_truth, 0);
    out.nhits_O.assign(n_truth, 0);
    out.nhits_T.assign(n_truth, 0);
    out.purity.assign(n_truth, 0.f);
    out.completeness.assign(n_truth, 0.f);

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

    for (size_t t = 0; t < n_truth; ++t) {
        out.nhits_T[t] = static_cast<int>(
            CountTruthHitsInSlice(tids[t], inputHits, assocMCPart));
    }

    const int n_pfp = static_cast<int>(pfp_hits.size());
    if (n_pfp < static_cast<int>(n_truth))
        return out;

    std::vector<std::vector<int>> shared(n_pfp, std::vector<int>(n_truth, 0));
    std::vector<int> n_hits_O(n_pfp, 0);

    for (int i = 0; i < n_pfp; ++i) {
        n_hits_O[i] = static_cast<int>(pfp_hits[i].size());
        for (auto const &h : pfp_hits[i]) {
            auto parts = assocMCPart->at(h.key());
            std::vector<bool> has(n_truth, false);
            for (auto const &mcp : parts) {
                const int tid = mcp->TrackId();
                for (size_t t = 0; t < n_truth; ++t) {
                    if (!has[t] && tid == tids[t]) {
                        has[t] = true;
                    }
                }
            }
            for (size_t t = 0; t < n_truth; ++t) {
                if (has[t]) {
                    ++shared[i][t];
                }
            }
        }
    }

    int best_sum = -1;
    std::vector<int> best_assignment(n_truth, -1);
    std::vector<int> current_assignment(n_truth, -1);
    std::vector<bool> used(n_pfp, false);

    std::function<void(size_t, int)> search = [&](size_t t, int sum) {
        int upper_bound = sum;
        for (size_t tt = t; tt < n_truth; ++tt) {
            int max_val = 0;
            for (int p = 0; p < n_pfp; ++p) {
                if (!used[p]) {
                    max_val = std::max(max_val, shared[p][tt]);
                }
            }
            upper_bound += max_val;
        }
        if (upper_bound <= best_sum)
            return;

        if (t == n_truth) {
            if (sum > best_sum) {
                best_sum = sum;
                best_assignment = current_assignment;
            }
            return;
        }

        for (int p = 0; p < n_pfp; ++p) {
            if (used[p])
                continue;
            used[p] = true;
            current_assignment[t] = p;
            search(t + 1, sum + shared[p][t]);
            current_assignment[t] = -1;
            used[p] = false;
        }
    };

    search(0, 0);
    if (best_sum < 0)
        return out;

    out.valid = true;

    auto fill_pair = [&](int p, size_t t) {
        if (p < 0)
            return;
        out.nshared[t] = shared[p][t];
        out.nhits_O[t] = n_hits_O[p];
        out.purity[t] = (out.nhits_O[t] > 0)
                            ? static_cast<float>(out.nshared[t]) /
                                  out.nhits_O[t]
                            : 0.f;
        out.completeness[t] = (out.nhits_T[t] > 0)
                                  ? static_cast<float>(out.nshared[t]) /
                                        out.nhits_T[t]
                                  : 0.f;
    };

    for (size_t t = 0; t < n_truth; ++t) {
        fill_pair(best_assignment[t], t);
    }

    return out;
}

} // namespace common

#endif
