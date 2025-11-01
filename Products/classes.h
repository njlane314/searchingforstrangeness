#include "canvas/Persistency/Common/Wrapper.h"

#include "Products/ImageProducts.h"
#include "Products/InferencePerf.h"
#include "Products/InferencePred.h"

#include <vector>

template class std::vector<image::PlaneImage>;
template class art::Wrapper<image::PlaneImage>;
template class art::Wrapper<std::vector<image::PlaneImage>>;

template class art::Wrapper<image::InferencePerfProduct>;

template class art::Wrapper<image::InferencePredProduct>;

template class std::vector<image::ModelPerf>;
template class art::Wrapper<std::vector<image::ModelPerf>>;

template class std::vector<image::ModelPred>;
template class art::Wrapper<std::vector<image::ModelPred>>;

