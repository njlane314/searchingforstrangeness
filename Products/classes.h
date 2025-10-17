#include "canvas/Persistency/Common/Wrapper.h"

#include "Products/ImageProducts.h"
#include "Products/InferencePerf.h"
#include "Products/RandomFeatures.h"
#include "Products/SegmentationProducts.h"

#include <vector>

template class std::vector<image::PlaneImage>;
template class art::Wrapper<std::vector<image::PlaneImage>>;

template class std::vector<image::PlaneSegmentation>;
template class art::Wrapper<std::vector<image::PlaneSegmentation>>;

template class art::Wrapper<image::InferencePerfProduct>;

template class std::vector<image::ModelPerf>;
template class art::Wrapper<std::vector<image::ModelPerf>>;

template class art::Wrapper<image::RandomFeatures>;

template class std::vector<image::RandomFeatures>;
template class art::Wrapper<std::vector<image::RandomFeatures>>;
