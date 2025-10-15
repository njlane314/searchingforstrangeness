#include "canvas/Persistency/Common/Wrapper.h"
#include "Products/ImageProducts.h"
#include "Products/SegmentationProducts.h"
#include "Products/InferencePerf.h"
#include <vector>

template class std::vector<analysis::PlaneImage>;
template class std::vector<analysis::PlaneSegmentation>;
template class std::vector<analysis::ModelPerf>;
