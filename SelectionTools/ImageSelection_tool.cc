#ifndef SELECTION_IMAGE_CXX
#define SELECTION_IMAGE_CXX

#include "SelectionToolBase.h"
#include <vector>
#include "art/Framework/Principal/Event.h"
#include "ImageProcessor.h"

namespace selection 
{
    class ImageSelection : public SelectionToolBase {
    public:
        ImageSelection(const fhicl::ParameterSet& pset);
        ~ImageSelection() {}

        void configure(fhicl::ParameterSet const& pset);
        bool selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& neutrino_slice, const std::vector<image::Image>& calo_images, const std::vector<image::Image>& reco_images, const std::vector<image::Image>& label_images);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);

    private:
        // Calorimetry images
        int                 calo_U_nRows, calo_U_nCols;
        int                 calo_V_nRows, calo_V_nCols;
        int                 calo_W_nRows, calo_W_nCols;
        std::vector<float>  calo_U_pixels;
        std::vector<float>  calo_V_pixels;
        std::vector<float>  calo_W_pixels;

        // Reconstructed images
        int                 reco_U_nRows, reco_U_nCols;
        int                 reco_V_nRows, reco_V_nCols;
        int                 reco_W_nRows, reco_W_nCols;
        std::vector<float>  reco_U_pixels;
        std::vector<float>  reco_V_pixels;
        std::vector<float>  reco_W_pixels;

        // Label images
        int                 label_U_nRows, label_U_nCols;
        int                 label_V_nRows, label_V_nCols;
        int                 label_W_nRows, label_W_nCols;
        std::vector<float>  label_U_pixels;
        std::vector<float>  label_V_pixels;
        std::vector<float>  label_W_pixels;
    };

    ImageSelection::ImageSelection(const fhicl::ParameterSet& pset) {}

    void ImageSelection::configure(fhicl::ParameterSet const& pset) {}

    void ImageSelection::setBranches(TTree* _tree) {
        _tree->Branch("calo_U_nRows", &calo_U_nRows, "calo_U_nRows/I");
        _tree->Branch("calo_U_nCols", &calo_U_nCols, "calo_U_nCols/I");
        _tree->Branch("calo_U_pixels", &calo_U_pixels);
        _tree->Branch("calo_V_nRows", &calo_V_nRows, "calo_V_nRows/I");
        _tree->Branch("calo_V_nCols", &calo_V_nCols, "calo_V_nCols/I");
        _tree->Branch("calo_V_pixels", &calo_V_pixels);
        _tree->Branch("calo_W_nRows", &calo_W_nRows, "calo_W_nRows/I");
        _tree->Branch("calo_W_nCols", &calo_W_nCols, "calo_W_nCols/I");
        _tree->Branch("calo_W_pixels", &calo_W_pixels);

        _tree->Branch("reco_U_nRows", &reco_U_nRows, "reco_U_nRows/I");
        _tree->Branch("reco_U_nCols", &reco_U_nCols, "reco_U_nCols/I");
        _tree->Branch("reco_U_pixels", &reco_U_pixels);
        _tree->Branch("reco_V_nRows", &reco_V_nRows, "reco_V_nRows/I");
        _tree->Branch("reco_V_nCols", &reco_V_nCols, "reco_V_nCols/I");
        _tree->Branch("reco_V_pixels", &reco_V_pixels);
        _tree->Branch("reco_W_nRows", &reco_W_nRows, "reco_W_nRows/I");
        _tree->Branch("reco_W_nCols", &reco_W_nCols, "reco_W_nCols/I");
        _tree->Branch("reco_W_pixels", &reco_W_pixels);

        _tree->Branch("label_U_nRows", &label_U_nRows, "label_U_nRows/I");
        _tree->Branch("label_U_nCols", &label_U_nCols, "label_U_nCols/I");
        _tree->Branch("label_U_pixels", &label_U_pixels);
        _tree->Branch("label_V_nRows", &label_V_nRows, "label_V_nRows/I");
        _tree->Branch("label_V_nCols", &label_V_nCols, "label_V_nCols/I");
        _tree->Branch("label_V_pixels", &label_V_pixels);
        _tree->Branch("label_W_nRows", &label_W_nRows, "label_W_nRows/I");
        _tree->Branch("label_W_nCols", &label_W_nCols, "label_W_nCols/I");
        _tree->Branch("label_W_pixels", &label_W_pixels);
    }

    void ImageSelection::resetTTree(TTree* _tree) {
        calo_U_pixels.clear();
        calo_V_pixels.clear();
        calo_W_pixels.clear();
        reco_U_pixels.clear();
        reco_V_pixels.clear();
        reco_W_pixels.clear();
        label_U_pixels.clear();
        label_V_pixels.clear();
        label_W_pixels.clear();
    }

    bool ImageSelection::selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& neutrino_slice, const std::vector<image::Image>& calo_images, const std::vector<image::Image>& reco_images, const std::vector<image::Image>& label_images) {
        if (calo_images.size() != 3 || reco_images.size() != 3 || label_images.size() != 3) {
            return false;
        }

        for (const auto& img : calo_images) {
            if (img.view() == geo::kU) {
                calo_U_nRows = img.height();
                calo_U_nCols = img.width();
                calo_U_pixels = img.data();
            } else if (img.view() == geo::kV) {
                calo_V_nRows = img.height();
                calo_V_nCols = img.width();
                calo_V_pixels = img.data();
            } else if (img.view() == geo::kW) {
                calo_W_nRows = img.height();
                calo_W_nCols = img.width();
                calo_W_pixels = img.data();
            } else {
                return false;
            }
        }

        for (const auto& img : reco_images) {
            if (img.view() == geo::kU) {
                reco_U_nRows = img.height();
                reco_U_nCols = img.width();
                reco_U_pixels = img.data();
            } else if (img.view() == geo::kV) {
                reco_V_nRows = img.height();
                reco_V_nCols = img.width();
                reco_V_pixels = img.data();
            } else if (img.view() == geo::kW) {
                reco_W_nRows = img.height();
                reco_W_nCols = img.width();
                reco_W_pixels = img.data();
            } else {
                return false;
            }
        }

        for (const auto& img : label_images) {
            if (img.view() == geo::kU) {
                label_U_nRows = img.height();
                label_U_nCols = img.width();
                label_U_pixels = img.data();
            } else if (img.view() == geo::kV) {
                label_V_nRows = img.height();
                label_V_nCols = img.width();
                label_V_pixels = img.data();
            } else if (img.view() == geo::kW) {
                label_W_nRows = img.height();
                label_W_nCols = img.width();
                label_W_pixels = img.data();
            } else {
                return false;
            }
        }

        return true;
    }

    DEFINE_ART_CLASS_TOOL(ImageSelection)
}

#endif