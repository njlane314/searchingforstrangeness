#ifndef SAMPLE_H
#define SAMPLE_H

#include <memory>
#include <string>
#include <stdexcept>
#include <TChain.h>
#include <TFile.h>

namespace analysis {
    enum class SampleType {
        GENIE, 
        DIRT, 
        COMSIC, 
        SIGNAL
    };

    class Sample {
    public:
        explicit Sample(SampleType type) 
            : type_(type), 
            event_chain_(std::make_unique<TChain>("selectionfilter/SelectionFilter")), 
            subrun_chain_(std::make_unique<TChain>("selectionfilter/SubRun")), 
            total_pot(0.0) {}

        void AddFile(const std::string& filename) {
            if (!TFile::Open(filename.c_str())) {
                throw std::runtime_error("-- File not found: " + filename);
            }
            event_chain_->Add(filename.c_str());
            subrun_chain_->Add(filename.c_str());
        }

        double GetTotalPOT() const {
            if (total_pot == 0.0) {
                this->AccumulatePOT();
            }
            return total_pot;
        }

        SampleType GetType() const {
            return type_;
        }

        TChain* GetEventChain() const { 
            return event_chain_.get(); 
        }

    private:
        SampleType type_; 
        std::unique_ptr<TChain> event_chain_;  
        std::unique_ptr<TChain> subrun_chain_; 
        mutable double total_pot;            

        void AccumulatePOT() const {
            float pot;
            subrun_chain_->SetBranchAddress("pot", &pot);
            Long64_t nSubruns = subrun_chain_->GetEntries();
            total_pot = 0.0;
            for (Long64_t i = 0; i < nSubruns; ++i) {
                subrun_chain_->GetEntry(i);
                total_pot += pot;
            }
        }
    };
}

#endif