#include <vector>

#include "canvas/Persistency/Common/Wrapper.h"

#include "Products/ImageFeatures.h"
#include "Products/InferenceMetrics.h"
#include "Products/InferencePredictions.h"

// Image feature products
template class std::vector<image::ImageFeatures>;
template class art::Wrapper<image::ImageFeatures>;
template class art::Wrapper<std::vector<image::ImageFeatures>>;

// Inference metrics products
template class std::vector<image::ModelMetrics>;
template class art::Wrapper<std::vector<image::ModelMetrics>>;
template class art::Wrapper<image::InferenceMetrics>;

// Inference prediction products
template class std::vector<image::ModelPrediction>;
template class art::Wrapper<std::vector<image::ModelPrediction>>;
template class art::Wrapper<image::InferencePredictions>;
