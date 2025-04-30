#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: tmva_training INPUT_FILE OUTPUT_FILE\n";
        return 1;
    }

    std::string input_file_name(argv[1]);
    std::string output_file_name(argv[2]);

    TFile* input_file = TFile::Open(input_file_name.c_str());
    if (!input_file || input_file->IsZombie()) {
        std::cout << "Error opening input file: " << input_file_name << "\n";
        return 1;
    }

    TTree* tree = dynamic_cast<TTree*>(input_file->Get("strangenessFilter/EventSelectionFilter"));
    if (!tree) {
        std::cout << "Cannot find strangeness filter in " << input_file_name << "\n";
        return 1;
    }

    TFile* output_file = TFile::Open(output_file_name.c_str(), "recreate");
    if (!output_file) {
        std::cout << "Cannot create " << output_file_name << "\n";
        return 1;
    }

    TMVA::Factory* factory = new TMVA::Factory("TMVAClassification", output_file, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

    dataloader->AddVariable("nhit_u", 'I');
    dataloader->AddVariable("nhit_v", 'I');
    dataloader->AddVariable("nhit_w", 'I');
    dataloader->AddVariable("charge_u", 'F');
    dataloader->AddVariable("charge_v", 'F');
    dataloader->AddVariable("charge_w", 'F');
    dataloader->AddVariable("wirerange_u", 'F');
    dataloader->AddVariable("wirerange_v", 'F');
    dataloader->AddVariable("wirerange_w", 'F');
    dataloader->AddVariable("timerange_u", 'F');
    dataloader->AddVariable("timerange_v", 'F');
    dataloader->AddVariable("timerange_w", 'F');

    TCut signal_cut = "in_fiducial && is_mc";
    TCut background_cut = "in_fiducial";

    int n_signal = tree->GetEntries(signal_cut);
    int n_background = tree->GetEntries(background_cut);

    dataloader->AddTree(tree, "Signal");
    dataloader->AddTree(tree, "Background");

    TString setup = Form("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:NormMode=NumEvents:!V", n_signal/2, n_background/2);
    dataloader->PrepareTrainingAndTestTree(signal_cut, background_cut, setup);

    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:MaxDepth=3");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    output_file->Close();

    delete factory;
    delete dataloader;
    delete output_file;
    delete input_file;

    return 0;
}