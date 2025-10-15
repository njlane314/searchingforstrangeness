#include "canvas/Persistency/Common/Wrapper.h"
#include "Analysis/ImageProducts.h"
#include "Analysis/SegmentationProducts.h"
#include <vector>

template class std::vector<analysis::PlaneImageProduct>;
template class std::vector<analysis::PlaneSegmentationProduct>;
