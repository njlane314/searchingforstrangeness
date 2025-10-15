#include "canvas/Persistency/Common/Wrapper.h"
#include "Products/ImageProducts.h"
#include "Products/SegmentationProducts.h"
#include <vector>

template class std::vector<analysis::PlaneImageProduct>;
template class std::vector<analysis::PlaneSegmentationProduct>;
