#ifndef PLOTTER_H
#define PLOTTER_H

#include <functional>
#include <vector>
#include <memory>
#include "TChain.h"

#include "Sample.h" 
#include "Event.h"  

namespace plot {
    template <typename HistType>
    class Plotter {
    private:
        const std::vector<std::unique_ptr<analysis::Sample>>& samples_;
        double target_pot_;

    public:
        Plotter(const std::vector<std::unique_ptr<analysis::Sample>>& samples, double target_pot)
            : samples_(samples), target_pot_(target_pot) {}

        void FillHistogram(std::function<void(const analysis::Event&, double)> filler) {
            for (const auto& sample : samples_) {
                double scale_factor = (sample->GetTotalPOT() > 0.0) ? target_pot_ / sample->GetTotalPOT() : 0.0;
                if (scale_factor == 0.0) 
                    continue;

                TChain* chain = sample->GetEventChain();
                analysis::Event event;
                event.SetBranches(chain);

                Long64_t entries = chain->GetEntries();
                for (Long64_t i = 0; i < entries; ++i) {
                    chain->GetEntry(i);
                    filler(event, scale_factor);
                }
            }
        }
    };
}

#endif