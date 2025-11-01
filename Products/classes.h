#include <vector>

#include "canvas/Persistency/Common/Wrapper.h"

#include "Products/ImageProducts.h"
#include "Products/InferencePerf.h"
#include "Products/InferencePred.h"

// Image products
template class std::vector<image::ImageProduct>;
template class art::Wrapper<image::ImageProduct>;
template class art::Wrapper<std::vector<image::ImageProduct>>;

// Inference performance products
template class std::vector<image::ModelPerf>;
template class art::Wrapper<std::vector<image::ModelPerf>>;
template class art::Wrapper<image::InferencePerfProduct>;

// Inference prediction products
template class std::vector<image::ModelPred>;
template class art::Wrapper<std::vector<image::ModelPred>>;
template class art::Wrapper<image::InferencePredProduct>;

