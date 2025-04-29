#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <iostream>

int main() {
    // Open input file
    TFile* inputFile = TFile::Open("/exp/uboone/data/users/nlane/analysis/nl_prodgeni_numi_uboone_overlay_rhc_training_1600_topologytrainingfilter_400_new_analysis.root");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file." << std::endl;
        return 1;
    }

    // Get the tree
    TTree* tree = (TTree*)inputFile->Get("topologyTrainingFilter/EventSelectionFilter");
    if (!tree) {
        std::cerr << "Error: Cannot find tree in input file." << std::endl;
        return 1;
    }

    // Create output file
    TFile* outputFile = TFile::Open("tmva_output.root", "RECREATE");
    if (!outputFile) {
        std::cerr << "Error: Cannot create output file." << std::endl;
        return 1;
    }

    // Initialize TMVA factory
    TMVA::Factory* factory = new TMVA::Factory("TMVAClassification", outputFile,
        "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

    // Initialize DataLoader
    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

    // Add variables
    dataloader->AddVariable("nhits_u", 'I');
    dataloader->AddVariable("nhit_v", 'I');
    dataloader->AddVariable("nhit_w", 'I');
    dataloader->AddVariable("charge_u", 'F');
    dataloader->AddVariable("charge_v", 'F');
    dataloader->AddVariable("chare_w", 'F');
    dataloader->AddVariable("wirerange_u", 'F');
    dataloader->AddVariable("wirerange_v", 'F');
    dataloader->AddVariable("wirerange_w", 'F');
    dataloader->AddVariable("timerange_u", 'F');
    dataloader->AddVariable("timerange_v", 'F');
    dataloader->AddVariable("timerange_w", 'F');

    // Define signal and background cuts
    TCut signalCut = "event_category >= 2 && event_category <= 8";
    TCut backgroundCut = "event_category == 1";

    // Calculate number of signal and background events
    int nSignal = tree->GetEntries(signalCut);
    int nBackground = tree->GetEntries(backgroundCut);

    // Add the tree for signal and background classes (without cuts here)
    dataloader->AddTree(tree, "Signal");
    dataloader->AddTree(tree, "Background");

    // Prepare training and test datasets with cuts
    TString setup = Form("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:NormMode=NumEvents:!V", nSignal/2, nBackground/2);
    dataloader->PrepareTrainingAndTestTree(signalCut, backgroundCut, setup);

    // Book BDT method (since XGBoost is not available in ROOT 6.12.06a)
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:MaxDepth=3");

    // Train, test, and evaluate
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    // Clean up
    outputFile->Close();
    delete factory;
    delete dataloader;
    delete inputFile;

    std::cout << "Training completed. Model saved in tmva_output.root" << std::endl;
    return 0;
}