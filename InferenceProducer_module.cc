#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ImagePipeline/InferenceConfiguration.h"
#include "ImagePipeline/InferencePlaneTriplet.h"
#include "ImagePipeline/InferenceProduction.h"
#include "Products/ImageFeatures.h"
#include "Products/InferenceMetrics.h"
#include "Products/InferencePredictions.h"

#include <exception>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class InferenceProducer : public art::EDProducer {
  public:
    explicit InferenceProducer(const fhicl::ParameterSet &parameters);
    void produce(art::Event &event) override;

  private:
    art::InputTag planes_tag_;
    image::InferenceConfiguration configuration_;
};

InferenceProducer::InferenceProducer(
    const fhicl::ParameterSet &parameters)
    : planes_tag_{parameters.get<std::string>("PlanesTag")},
      configuration_{parameters} {
    produces<image::InferenceMetrics>();
    produces<image::InferencePredictions>();
}

void InferenceProducer::produce(art::Event &event) {
    auto metrics = std::make_unique<image::InferenceMetrics>();
    auto predictions = std::make_unique<image::InferencePredictions>();
    const auto publish = [&]() {
        event.put(std::move(metrics));
        event.put(std::move(predictions));
    };

    const auto planes =
        event.getValidHandle<std::vector<image::ImageFeatures>>(planes_tag_);
    if (planes->size() < 3U) {
        mf::LogDebug("InferenceProduction")
            << "Need three planes (U,V,W), got " << planes->size()
            << "; skipping inference for event " << event.id();
        publish();
        return;
    }

    image::ValidatedPlaneTriplet detector_images = [&]() {
        try {
            return image::ValidatedPlaneTriplet(*planes);
        } catch (const std::exception &error) {
            throw cet::exception("InferenceProduction")
                << "Invalid image input for event " << event.id() << ": "
                << error.what();
        }
    }();

    if (detector_images.empty()) {
        mf::LogDebug("InferenceProduction")
            << "Detector images are empty; skipping inference for event "
            << event.id();
        publish();
        return;
    }

    std::ostringstream event_context;
    event_context << "event " << event.id();
    for (const auto &model : configuration_.models()) {
        const auto run_config = configuration_.makeRunConfig(
            model, event_context.str());

        image::InferenceProduction::Result result;
        try {
            result = image::InferenceProduction::runInference(
                detector_images, run_config);
        } catch (const std::exception &error) {
            throw cet::exception("InferenceProduction") << error.what();
        }

        image::ModelMetrics model_metrics;
        model_metrics.model = model.productName();
        model_metrics.t_exec_total_ms =
            static_cast<float>(result.perf.t_exec_total_ms);
        model_metrics.t_child_infer_ms =
            static_cast<float>(result.perf.t_child_infer_ms);
        model_metrics.child_max_rss_mb =
            static_cast<float>(result.perf.child_max_rss_mb);
        metrics->per_model.push_back(std::move(model_metrics));

        image::ModelPrediction prediction;
        prediction.model = model.productName();
        prediction.scores = std::move(result.cls);
        predictions->per_model.push_back(std::move(prediction));
    }

    publish();
}

DEFINE_ART_MODULE(InferenceProducer)
