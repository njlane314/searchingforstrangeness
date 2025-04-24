#ifndef SELECTION_IMAGE_CXX
#define SELECTION_IMAGE_CXX

#include "SelectionToolBase.h"
#include <vector>
#include "art/Framework/Principal/Event.h"
#include "../CommonDefs/Image.h"

namespace selection 
{
    class ImageSelection : public SelectionToolBase {
    public:
        ImageSelection(const fhicl::ParameterSet& pset);
        ~ImageSelection() {}

        void configure(fhicl::ParameterSet const& pset);
        bool selectEvent(art::Event const& e, 
                        const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                        const std::vector<image::Image<float>>& calo_images, 
                        const std::vector<image::Image<int>>& reco_images, 
                        const std::vector<image::Image<int>>& label_images);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);

    private:
        // Calorimetry images
        std::vector<float>  calo_pixels_u;
        std::vector<float>  calo_pixels_v;
        std::vector<float>  calo_pixels_w;

        // Reconstructed images
        std::vector<int> reco_pixels_u;
        std::vector<int> reco_pixels_v;
        std::vector<int> reco_pixels_w;

        // Label images
        std::vector<int> label_pixels_u;
        std::vector<int> label_pixels_v;
        std::vector<int> label_pixels_w;
    };

    ImageSelection::ImageSelection(const fhicl::ParameterSet& pset) {}

    void ImageSelection::configure(fhicl::ParameterSet const& pset) {}

    void ImageSelection::setBranches(TTree* _tree) {
        _tree->Branch("calo_pixels_u", &calo_pixels_u);
        _tree->Branch("calo_pixels_v", &calo_pixels_v);
        _tree->Branch("calo_pixels_w", &calo_pixels_w);

        _tree->Branch("reco_pixels_u", &reco_pixels_u);
        _tree->Branch("reco_pixels_v", &reco_pixels_v);
        _tree->Branch("reco_pixels_w", &reco_pixels_w);

        _tree->Branch("label_pixels_u", &label_pixels_u);
        _tree->Branch("label_pixels_v", &label_pixels_v);
        _tree->Branch("label_pixels_w", &label_pixels_w);
    }

    void ImageSelection::resetTTree(TTree* _tree) {
        calo_pixels_u.clear();
        calo_pixels_v.clear();
        calo_pixels_w.clear();
        reco_pixels_u.clear();
        reco_pixels_v.clear();
        reco_pixels_w.clear();
        label_pixels_u.clear();
        label_pixels_v.clear();
        label_pixels_w.clear();
    }

    bool ImageSelection::selectEvent(art::Event const& e, 
                                    const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                    const std::vector<image::Image<float>>& calo_images, 
                                    const std::vector<image::Image<int>>& reco_images, 
                                    const std::vector<image::Image<int>>& label_images) {
        if (calo_images.size() != 3 || reco_images.size() != 3 || label_images.size() != 3) {
            return false;
        }

        for (const auto& img : calo_images) {
            if (img.view() == geo::kU) {
                calo_pixels_u = img.data();
            } else if (img.view() == geo::kV) {
                calo_pixels_v = img.data();
            } else if (img.view() == geo::kW) {
                calo_pixels_w = img.data();
            } else {
                return false;
            }
        }

        for (const auto& img : reco_images) {
            if (img.view() == geo::kU) {
                reco_pixels_u = img.data();
            } else if (img.view() == geo::kV) {
                reco_pixels_v = img.data();
            } else if (img.view() == geo::kW) {
                reco_pixels_w = img.data();
            } else {
                return false;
            }
        }

        for (const auto& img : label_images) {
            if (img.view() == geo::kU) {
                label_pixels_u = img.data();
            } else if (img.view() == geo::kV) {
                label_pixels_v = img.data();
            } else if (img.view() == geo::kW) {
                label_pixels_w = img.data();
            } else {
                return false;
            }
        }

        return true;
    }

    DEFINE_ART_CLASS_TOOL(ImageSelection)
}

#endif