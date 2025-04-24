//root -l -b -q 'analyse_root.C("/exp/uboone/data/users/nlane/analysis/prod_strange_resample_fhc_run2_fhc_reco2_reco2_selectionfilter_5_new_analysis.root")'
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TList.h>
#include <TObject.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TH1.h>
#include <iostream>

void analyse_directory(TDirectory* dir, int depth = 0) {
    TList* keys = dir->GetListOfKeys();
    if (!keys || keys->GetSize() == 0) {
        std::cout << std::string(depth * 2, ' ') << "  -> Empty directory\n";
        return;
    }

    TIter next(keys);
    TKey* key;
    int objCount = 0;

    while ((key = (TKey*)next())) {
        objCount++;
        TString className = key->GetClassName();
        TString objName = key->GetName();
        std::cout << std::string(depth * 2, ' ') << "Object " << objCount
                  << ": Name = " << objName << ", Type = " << className << "\n";

        if (className == "TTree") {
            TTree* tree = (TTree*)dir->Get(objName);
            if (tree) {
                std::cout << std::string(depth * 2, ' ') << "  -> TTree with "
                          << tree->GetEntries() << " entries\n";
                std::cout << std::string(depth * 2, ' ') << "  -> Branches:\n";

                TObjArray* branches = tree->GetListOfBranches();
                for (int i = 0; i < branches->GetEntries(); ++i) {
                    TBranch* branch = (TBranch*)branches->At(i);
                    std::cout << std::string(depth * 2, ' ') << "    - "
                              << branch->GetName() << " (Type: "
                              << branch->GetClassName() << ")\n";
                }
            }
        }
        else if (className.Contains("TH1")) {
            TH1* hist = (TH1*)dir->Get(objName);
            if (hist) {
                std::cout << std::string(depth * 2, ' ') << "  -> Histogram with "
                          << hist->GetEntries() << " entries\n";
                std::cout << std::string(depth * 2, ' ') << "     X-axis: "
                          << hist->GetNbinsX() << " bins\n";
            }
        }
        else if (className == "TDirectoryFile") {
            TDirectoryFile* subDir = (TDirectoryFile*)dir->Get(objName);
            if (subDir) {
                std::cout << std::string(depth * 2, ' ') << "  -> Entering subdirectory: " << objName << "\n";
                analyse_directory(subDir, depth + 1);
            }
        }
    }

    if (depth > 0) {
        std::cout << std::string(depth * 2, ' ') << "Total objects in directory: " << objCount << "\n";
    }
}

void analyse_root(const char* filename) {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    std::cout << "\n=== Analyzing ROOT file: " << filename << " ===\n";
    std::cout << "File size: " << file->GetSize() << " bytes\n";

    std::cout << "\n=== Objects in the file ===\n";
    analyse_directory(file);

    file->Close();
    delete file;
}

#ifndef __CINT__
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root_file>\n";
        return 1;
    }
    analyse_root(argv[1]);
    return 0;
}
#endif