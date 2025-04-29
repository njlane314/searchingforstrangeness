#include "Display.h"
#include "Config.h"
#include "Sample.h"
#include "Event.h"
#include "Plotter.h"

#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <utility>
#include <functional>
#include <TRandom3.h>

int main() {
    std::string config_name = "default";
    auto config = analysis::Config::GetConfig(config_name);
    std::vector<std::unique_ptr<analysis::Sample>> samples;

    if (config.empty())
        return 1;

    for (const auto& [type, files] : config) {
        auto sample = std::make_unique<analysis::Sample>(type);
        for (const auto& file : files) {
            try {
                sample->AddFile(file);
            } catch (const std::runtime_error& e) {
                continue;
            }
        }
        samples.push_back(std::move(sample));
    }

    int img_size = 512;
    std::string output_dir = "./displays";
    plot::Display display(img_size, output_dir);

    if (!samples.empty()) {
        const analysis::Sample& sample = *samples[0];
        TChain* chain = sample.GetEventChain();
        analysis::Event event;
        event.SetBranches(chain);

        TRandom3 rand(0);
        int random_indice = rand.Integer(chain->GetEntries());
        display.VisualiseInput(sample, random_indice);
        display.VisualiseLabels(sample, random_indice);
        display.VisualiseTruth(sample, random_indice);
    }

    display.PlotLabelsLegend();
    display.PlotTruthLegend();

    return 0;
}