#include <iostream>
#include <TFile.h>
#include <TTree.h>
//    root -l
//    TFile *file = TFile::Open("<>");
//    file->ls();
//    TDirectory *dir = (TDirectory*)file->Get("<>>");
//    dir->ls();
//    TTree* tree = (TTree*)dir->Get("<>");


// root -l -b -q 'event_dump.C("neutrino_selection_output.root")'
void event_dump(const char* filename) {
 
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Access the nuselection directory
    TDirectory *dir = (TDirectory*)file->Get("nuselection");
    if (!dir) {
        std::cerr << "Error: nuselection directory not found in file: " << filename << std::endl;
        file->Close();
        return;
    }

    // Access the NeutrinoSelectionFilter TTree
    TTree* tree = (TTree*)dir->Get("NeutrinoSelectionFilter");
    if (!tree) {
        std::cerr << "Error: NeutrinoSelectionFilter TTree not found in nuselection directory" << std::endl;
        file->Close();
        return;
    }

    // Set branch addresses
    unsigned int mc_piplus_n_elas;
    unsigned int mc_piplus_n_inelas;
    unsigned int mc_piminus_n_elas;
    unsigned int mc_piminus_n_inelas;

    int mc_kshrt_piplus_pdg;
    float mc_kshrt_piplus_energy;
    float mc_kshrt_piplus_px, mc_kshrt_piplus_py, mc_kshrt_piplus_pz;

    int mc_kshrt_piminus_pdg;
    float mc_kshrt_piminus_energy;
    float mc_kshrt_piminus_px, mc_kshrt_piminus_py, mc_kshrt_piminus_pz;

    float mc_kshrt_total_energy;
    int mc_num_daughters;

    float mc_neutrino_vertex_x, mc_neutrino_vertex_y, mc_neutrino_vertex_z;
    float mc_kaon_decay_x, mc_kaon_decay_y, mc_kaon_decay_z;
    float mc_kaon_decay_distance;

    float mc_piplus_impact_param;
    float mc_piminus_impact_param;

    std::string* mc_piplus_endprocess = nullptr;
    std::string* mc_piminus_endprocess = nullptr;

    // Set branch addresses
    tree->SetBranchAddress("mc_piplus_n_elas", &mc_piplus_n_elas);
    tree->SetBranchAddress("mc_piplus_n_inelas", &mc_piplus_n_inelas);
    tree->SetBranchAddress("mc_piminus_n_elas", &mc_piminus_n_elas);
    tree->SetBranchAddress("mc_piminus_n_inelas", &mc_piminus_n_inelas);

    tree->SetBranchAddress("mc_kshrt_piplus_pdg", &mc_kshrt_piplus_pdg);
    tree->SetBranchAddress("mc_kshrt_piplus_energy", &mc_kshrt_piplus_energy);
    tree->SetBranchAddress("mc_kshrt_piplus_px", &mc_kshrt_piplus_px);
    tree->SetBranchAddress("mc_kshrt_piplus_py", &mc_kshrt_piplus_py);
    tree->SetBranchAddress("mc_kshrt_piplus_pz", &mc_kshrt_piplus_pz);

    tree->SetBranchAddress("mc_kshrt_piminus_pdg", &mc_kshrt_piminus_pdg);
    tree->SetBranchAddress("mc_kshrt_piminus_energy", &mc_kshrt_piminus_energy);
    tree->SetBranchAddress("mc_kshrt_piminus_px", &mc_kshrt_piminus_px);
    tree->SetBranchAddress("mc_kshrt_piminus_py", &mc_kshrt_piminus_py);
    tree->SetBranchAddress("mc_kshrt_piminus_pz", &mc_kshrt_piminus_pz);

    tree->SetBranchAddress("mc_kshrt_total_energy", &mc_kshrt_total_energy);
    tree->SetBranchAddress("mc_num_daughters", &mc_num_daughters);

    tree->SetBranchAddress("mc_neutrino_vertex_x", &mc_neutrino_vertex_x);
    tree->SetBranchAddress("mc_neutrino_vertex_y", &mc_neutrino_vertex_y);
    tree->SetBranchAddress("mc_neutrino_vertex_z", &mc_neutrino_vertex_z);

    tree->SetBranchAddress("mc_kaon_decay_x", &mc_kaon_decay_x);
    tree->SetBranchAddress("mc_kaon_decay_y", &mc_kaon_decay_y);
    tree->SetBranchAddress("mc_kaon_decay_z", &mc_kaon_decay_z);
    tree->SetBranchAddress("mc_kaon_decay_distance", &mc_kaon_decay_distance);

    tree->SetBranchAddress("mc_piplus_impact_param", &mc_piplus_impact_param);
    tree->SetBranchAddress("mc_piminus_impact_param", &mc_piminus_impact_param);

    tree->SetBranchAddress("mc_piplus_endprocess", &mc_piplus_endprocess);
    tree->SetBranchAddress("mc_piminus_endprocess", &mc_piminus_endprocess);

    // Loop over the entries in the TTree
    for (Long64_t i = 0; i < 10; i++) {
        tree->GetEntry(i);

        // Dump the data (this is where you can customize what you want to do with the data)
        std::cout << "Event " << i << std::endl;
        std::cout << "  mc_piplus_n_elas: " << mc_piplus_n_elas << std::endl;
        std::cout << "  mc_piplus_n_inelas: " << mc_piplus_n_inelas << std::endl;
        std::cout << "  mc_piminus_n_elas: " << mc_piminus_n_elas << std::endl;
        std::cout << "  mc_piminus_n_inelas: " << mc_piminus_n_inelas << std::endl;

        std::cout << "  mc_kshrt_piplus_pdg: " << mc_kshrt_piplus_pdg << std::endl;
        std::cout << "  mc_kshrt_piplus_energy: " << mc_kshrt_piplus_energy << std::endl;
        std::cout << "  mc_kshrt_piplus_px: " << mc_kshrt_piplus_px << std::endl;
        std::cout << "  mc_kshrt_piplus_py: " << mc_kshrt_piplus_py << std::endl;
        std::cout << "  mc_kshrt_piplus_pz: " << mc_kshrt_piplus_pz << std::endl;

        std::cout << "  mc_kshrt_piminus_pdg: " << mc_kshrt_piminus_pdg << std::endl;
        std::cout << "  mc_kshrt_piminus_energy: " << mc_kshrt_piminus_energy << std::endl;
        std::cout << "  mc_kshrt_piminus_px: " << mc_kshrt_piminus_px << std::endl;
        std::cout << "  mc_kshrt_piminus_py: " << mc_kshrt_piminus_py << std::endl;
        std::cout << "  mc_kshrt_piminus_pz: " << mc_kshrt_piminus_pz << std::endl;

        std::cout << "  mc_kshrt_total_energy: " << mc_kshrt_total_energy << std::endl;
        std::cout << "  mc_num_daughters: " << mc_num_daughters << std::endl;

        std::cout << "  mc_neutrino_vertex_x: " << mc_neutrino_vertex_x << std::endl;
        std::cout << "  mc_neutrino_vertex_y: " << mc_neutrino_vertex_y << std::endl;
        std::cout << "  mc_neutrino_vertex_z: " << mc_neutrino_vertex_z << std::endl;

        std::cout << "  mc_kaon_decay_x: " << mc_kaon_decay_x << std::endl;
        std::cout << "  mc_kaon_decay_y: " << mc_kaon_decay_y << std::endl;
        std::cout << "  mc_kaon_decay_z: " << mc_kaon_decay_z << std::endl;
        std::cout << "  mc_kaon_decay_distance: " << mc_kaon_decay_distance << std::endl;

        std::cout << "  mc_piplus_impact_param: " << mc_piplus_impact_param << std::endl;
        std::cout << "  mc_piminus_impact_param: " << mc_piminus_impact_param << std::endl;

        std::cout << "  mc_piplus_endprocess: " << *mc_piplus_endprocess << std::endl;
        std::cout << "  mc_piminus_endprocess: " << *mc_piminus_endprocess << std::endl;
    }

    // Clean up
    file->Close();
    delete file;
}
