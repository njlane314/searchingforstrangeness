#include "canvas/Persistency/Common/Wrapper.h"
#include "DataProducts/ImageProducts.h"
#include "DataProducts/SegmentationProducts.h"
#include <vector>

template class std::vector<analysis::PlaneImageProduct>;
template class std::vector<analysis::PlaneSegmentationProduct>;
