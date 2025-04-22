#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <iostream>
#include "TMVA/TMVAGui.h"

int main() {
    TFile* inputFile = TFile::Open("input.root");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Cannot open input file." << std::endl;
        return 1;
    }
    TTree* tree = (TTree*)inputFile->Get("selection_tree");
    if (!tree) {
        std::cerr << "Error: Cannot find tree in input file." << std::endl;
        return 1;
    }
    TFile* outputFile = TFile::Open("tmva_output.root", "RECREATE");
    if (!outputFile) {
        std::cerr << "Error: Cannot create output file." << std::endl;
        return 1;
    }
    TMVA::Factory* factory = new TMVA::Factory("TMVAClassification", outputFile,
        "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
    TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");
    dataloader->AddVariable("nHitsU", 'I');
    dataloader->AddVariable("nHitsV", 'I');
    dataloader->AddVariable("nHitsW", 'I');
    dataloader->AddVariable("totalChargeU", 'F');
    dataloader->AddVariable("totalChargeV", 'F');
    dataloader->AddVariable("totalChargeW", 'F');
    dataloader->AddVariable("wireRangeU", 'F');
    dataloader->AddVariable("wireRangeV", 'F');
    dataloader->AddVariable("wireRangeW", 'F');
    dataloader->AddVariable("timeRangeU", 'F');
    dataloader->AddVariable("timeRangeV", 'F');
    dataloader->AddVariable("timeRangeW", 'F');
    TCut signalCut = "event_category >= 2 && event_category <= 8";
    TCut backgroundCut = "event_category == 1";
    int nSignal = tree->GetEntries(signalCut);
    int nBackground = tree->GetEntries(backgroundCut);
    dataloader->AddSignalTree(tree, 1.0, signalCut);
    dataloader->AddBackgroundTree(tree, 1.0, backgroundCut);
    TString setup = Form("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:NormMode=NumEvents:!V", nSignal/2, nBackground/2);
    dataloader->PrepareTrainingAndTestTree("", setup);
    factory->BookMethod(dataloader, TMVA::Types::kXGBoost, "XGBoost", "!H:!V:objective=binary:logistic:eta=0.1:max_depth=6:num_round=100");
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    outputFile->Close();
    delete factory;
    delete dataloader;
    delete inputFile;
    std::cout << "Training completed. Model saved in tmva_output.root" << std::endl;
    if (!gROOT->IsBatch()) TMVA::TMVAGui("tmva_output.root");
    return 0;
}