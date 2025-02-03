#ifndef IMAGEINFERENCER_H
#define IMAGEINFERENCER_H

#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Wire.h"
#include "Geometry/GeometryCore.h"
#include "ImageProcessor.h" 
#include <torch/torch.h>
#include <torch/script.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace image {

class ImageInferencer {
public:
  explicit ImageInferencer(const fhicl::ParameterSet &pset, const geo::GeometryCore* geo)
    : _geo(geo)
  {
    _image_width    = pset.get<int>("ImageWidth", 256);
    _image_height   = pset.get<int>("ImageHeight", 256);
    _veto_bad_channels = pset.get<bool>("VetoBadChannels", true);
    _bad_channel_file  = pset.get<std::string>("BadChannelFile", "badchannels.txt");

    auto _det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    _drift_step = (_det_props->SamplingRate() / 1000.) * 
                  _det_props->DriftVelocity(_det_props->Efield(), _det_props->Temperature());

    _wire_pitch_u = _geo->WirePitch(geo::kU);
    _wire_pitch_v = _geo->WirePitch(geo::kV);
    _wire_pitch_w = _geo->WirePitch(geo::kW);

    try {
      _model_u = std::make_shared<torch::jit::script::Module>(
                     torch::jit::load(pset.get<std::string>("ModelFileU")));
      _model_v = std::make_shared<torch::jit::script::Module>(
                     torch::jit::load(pset.get<std::string>("ModelFileV")));
      _model_w = std::make_shared<torch::jit::script::Module>(
                     torch::jit::load(pset.get<std::string>("ModelFileW")));
    } catch (const c10::Error & e) {
      throw cet::exception("ImageInferencer") 
            << "Error loading torch models: " << e.what();
    }

    if (_veto_bad_channels)
      initialiseBadChannelMask();
  }

  ImageInferencer(const ImageInferencer&) = delete;
  ImageInferencer(ImageInferencer&&) = delete;
  ImageInferencer& operator=(const ImageInferencer&) = delete;
  ImageInferencer& operator=(ImageInferencer&&) = delete;

  std::vector<torch::Tensor> inferWireImages(const std::vector<art::Ptr<recob::Wire>>& wires,
                                             const std::vector<image::ImageMeta>& metas)
  {
    std::vector<art::Ptr<recob::Wire>> filteredWires = wires;
    if (_veto_bad_channels)
      filterBadChannels(filteredWires);

    std::vector<image::Image> images = image::WiresToImages(metas, filteredWires, *_geo);

    std::vector<torch::Tensor> outputs;

    for (const auto &img : images) {
      torch::Tensor tensor = img.tensor();
      geo::View_t view = img.getMeta().view();
      torch::Tensor output;
      if (view == geo::kU) {
        output = _model_u->forward({tensor}).tensor();
      } else if (view == geo::kV) {
        output = _model_v->forward({tensor}).tensor();
      } else if (view == geo::kW) {
        output = _model_w->forward({tensor}).tensor();
      } else {
        output = tensor;
      }
      outputs.push_back(output);
    }
    return outputs;
  }

  std::shared_ptr<torch::jit::script::Module> modelU() const { return _model_u; }
  std::shared_ptr<torch::jit::script::Module> modelV() const { return _model_v; }
  std::shared_ptr<torch::jit::script::Module> modelW() const { return _model_w; }

private:
  int _image_width;
  int _image_height;
  float _drift_step;
  float _wire_pitch_u;
  float _wire_pitch_v;
  float _wire_pitch_w;
  bool _veto_bad_channels;
  std::string _bad_channel_file;
  std::vector<bool> _bad_channel_mask;
  const geo::GeometryCore* _geo;

  std::shared_ptr<torch::jit::script::Module> _model_u;
  std::shared_ptr<torch::jit::script::Module> _model_v;
  std::shared_ptr<torch::jit::script::Module> _model_w;

  void initialiseBadChannelMask() {
    size_t n_channels = _geo->Nchannels();
    _bad_channel_mask.resize(n_channels, false);
    if (!_bad_channel_file.empty()) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string fullname;
      sp.find_file(_bad_channel_file, fullname);
      if (fullname.empty()) {
        throw cet::exception("ImageInferencer") 
              << "**bad channel file not found: " << _bad_channel_file;
      }
      std::ifstream inFile(fullname, std::ios::in);
      std::string line;
      while (std::getline(inFile, line)) {
        if (line.find("#") != std::string::npos) continue;
        std::istringstream ss(line);
        int ch1, ch2;
        ss >> ch1;
        if (!(ss >> ch2)) ch2 = ch1;
        for (int i = ch1; i <= ch2; ++i) {
          if (i < static_cast<int>(_bad_channel_mask.size()))
            _bad_channel_mask[i] = true;
        }
      }
    }
  }

  void filterBadChannels(std::vector<art::Ptr<recob::Wire>>& wires) {
    wires.erase(std::remove_if(wires.begin(), wires.end(),
      [this](const art::Ptr<recob::Wire>& wire) {
        return _bad_channel_mask[wire->Channel()];
      }),
      wires.end());
  }
};

} 

#endif // IMAGEINFERENCER_H
