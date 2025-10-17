#include "canvas/Persistency/Common/Wrapper.h"
#include "Products/ImageProducts.h"
#include "Products/SegmentationProducts.h"
#include "Products/InferencePerf.h"
#include "Products/RandomFeatures.h"
#include <vector>

template class std::vector<image::PlaneImage>;
template class std::vector<image::PlaneSegmentation>;
template class std::vector<image::ModelPerf>;
template class std::vector<image::RandomFeatures>;
