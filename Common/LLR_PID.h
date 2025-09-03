#ifndef LLRPID_H
#define LLRPID_H

#include "BacktrackingFuncs.h"
#include "Geometry.h"

namespace common {

class LLRPID {
    public:
    LLRPID() {}

    void set_dedx_binning(size_t plane, std::vector<float> bin_edges) {
        dedx_bin_edges[plane] = bin_edges;
        dedx_num_bins[plane] = (bin_edges.size() - 1);
    }

    void set_par_binning(size_t plane, std::vector<std::vector<float>> bin_edges) {
        parameters_bin_edges[plane] = bin_edges;
        for (size_t i = 0; i < bin_edges.size(); i++) {
            parameters_num_bins[plane].push_back(parameters_bin_edges[plane][i].size() - 1);
        }
    }

    void set_corr_par_binning(size_t plane, std::vector<std::vector<float>> bin_edges) {
        corr_parameters_bin_edges[plane] = bin_edges;
        parameters_bin_edges[plane] = bin_edges;
        for (size_t i = 0; i < bin_edges.size(); i++) {
            corr_parameters_num_bins[plane].push_back(corr_parameters_bin_edges[plane][i].size() - 1);
        }
    }

    void set_lookup_tables(size_t plane, std::vector<float> tables) {
        lookup_tables[plane] = tables;
    }

    void set_correction_tables(size_t plane, std::vector<float> tables) {
        correction_tables[plane] = tables;
    }

    size_t digitize(float value, std::vector<float> bin_edges) {
        if (value <= bin_edges[0])
            return 0;
        for (size_t i = 0; i < bin_edges.size(); i++) {
            if (value >= bin_edges[i])
                continue;
            else
                return i - 1;
        }
        return bin_edges.size() - 1;
    }

    size_t findLookupIndex(float dedx_value, std::vector<float> par_value, size_t plane) {

        std::vector<size_t> this_parameters_bins;
        for (size_t i = 0; i < par_value.size(); i++) {
            size_t aux_index = digitize(par_value[i], parameters_bin_edges[plane][i]);
            this_parameters_bins.push_back(aux_index);
        }

        size_t lookup_row = 0, accumulator_par_bins = 1;
        for (size_t i = this_parameters_bins.size(); i-- > 0;) {
            lookup_row += (accumulator_par_bins * this_parameters_bins[i]);
            accumulator_par_bins *= parameters_num_bins[plane][i];
        }

        size_t lookup_row_index;
        lookup_row_index = lookup_row * dedx_num_bins[plane];

        size_t lookup_index = lookup_row_index;
        lookup_index += digitize(dedx_value, dedx_bin_edges[plane]);

        return lookup_index;
    }

    size_t findLookupCorrParameterIndex(std::vector<float> corr_parameter_value, size_t plane) {

        std::vector<size_t> this_corr_parameters_bins;
        for (size_t i = 0; i < corr_parameter_value.size(); i++) {
            size_t aux_index = digitize(corr_parameter_value[i], corr_parameters_bin_edges[plane][i]);
            this_corr_parameters_bins.push_back(aux_index);
        }

        size_t lookup_index = 0, accumulator_par_bins = 1;
        for (size_t i = this_corr_parameters_bins.size(); i-- > 0;) {
            lookup_index += (accumulator_par_bins * this_corr_parameters_bins[i]);
            accumulator_par_bins *= corr_parameters_num_bins[plane][i];
        }

        return lookup_index;
    }

    float LLR_one_hit_one_plane(float dedx_value, std::vector<float> par_value, size_t plane) {
        size_t index = findLookupIndex(dedx_value, par_value, plane);
        return lookup_tables[plane][index];
    }

    float correction_hit_one_plane(std::vector<float> corr_parameter_value, size_t plane) {
        size_t index = findLookupCorrParameterIndex(corr_parameter_value, plane);
        return correction_tables[plane][index];
    }

    float LLR_many_hits_one_plane(std::vector<float> dedx_values, std::vector<std::vector<float>> par_values, size_t plane) {
        float ll_out = 0;
        for (size_t i = 0; i < dedx_values.size(); i++) {
            std::vector<float> aux_par;
            for (std::vector<float> par_value : par_values) {
                aux_par.push_back(par_value[i]);
            }
            ll_out += LLR_one_hit_one_plane(dedx_values[i], aux_par, plane);
        }
        return ll_out;
    }

    std::vector<float> correct_many_hits_one_plane(std::vector<float> dedx_values, std::vector<std::vector<float>> corr_par_values, std::vector<bool> is_to_correct, size_t plane) {
        std::vector<float> dedx_values_corrected;
        for (size_t i = 0; i < dedx_values.size(); i++) {
            float aux_dedx = dedx_values[i];
            if (is_to_correct[i]) {
                aux_dedx *= correction_hit_one_plane(corr_par_values[i], plane);
            }
            dedx_values_corrected.push_back(aux_dedx);
        }
        return dedx_values_corrected;
    }

    std::vector<float> correct_many_hits_one_plane(const art::Ptr<anab::Calorimetry> tkcalo,
                                                   const recob::Track trk,
                                                   const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart,
                                                   const bool fRecalibrateHits,
                                                   const float fEnergyThresholdForMCHits,
                                                   const bool fLocaldEdx) {
        int plane = tkcalo->PlaneID().Plane;

        std::vector<float> dqdx_values, dqdx_values_corrected;
        if (!fLocaldEdx)
            dqdx_values = tkcalo->dQdx();
        else
            dqdx_values = tkcalo->dEdx();
        if (!fRecalibrateHits)
            dqdx_values_corrected = dqdx_values;
        else {
            auto const &pitch = tkcalo->TrkPitchVec();
            auto const &xyz_v = tkcalo->XYZ();

            std::vector<std::vector<float>> corr_par_values;
            for (size_t i = 0; i < xyz_v.size(); i++) {
                auto xyz = xyz_v[i];
                float _dir[3];
                common::TrkDirectionAtXYZ(trk, xyz.X(), xyz.Y(), xyz.Z(), _dir);
                std::vector<float> corr_par_value = common::polarAngles(_dir[0], _dir[1], _dir[2], 2, plane);
                corr_par_value[0] = pitch[i];
                corr_par_values.push_back(corr_par_value);
            }

            std::vector<bool> is_hit_montecarlo;
            const std::vector<size_t> &tp_indices = tkcalo->TpIndices();
            for (size_t i = 0; i < tp_indices.size(); i++) {
                size_t tp_index = tp_indices[i];
                is_hit_montecarlo.push_back(common::isHitBtMonteCarlo(tp_index, assocMCPart, fEnergyThresholdForMCHits));
            }

            dqdx_values_corrected = correct_many_hits_one_plane(dqdx_values, corr_par_values, is_hit_montecarlo, plane);
        }
        return dqdx_values_corrected;
    }

    private:
    size_t dedx_num_bins[3];
    std::vector<float> dedx_bin_edges[3];

    std::vector<size_t> parameters_num_bins[3];
    std::vector<std::vector<float>> parameters_bin_edges[3];

    std::vector<size_t> corr_parameters_num_bins[3];
    std::vector<std::vector<float>> corr_parameters_bin_edges[3];

    std::vector<float> lookup_tables[3];
    std::vector<float> correction_tables[3];
};
} // namespace common

#endif