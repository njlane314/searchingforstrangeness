#include <iostream>
#include <string>
#include <vector>

#include "gallery/Event.h"
#include "evwgh/MCEventWeight.h"

void list_weights(std::string filename,
                  std::string eventweight_producer_label="eventweight")
{
  std::vector<std::string> filenames { filename };
  gallery::Event ev(filenames);

  ev.toBegin();
  if (ev.atEnd()) { std::cout << "Empty file?\n"; return; }

  auto weights_handle =
    ev.getValidHandle<std::vector<evwgh::MCEventWeight>>(eventweight_producer_label);

  if (weights_handle->empty()) {
    std::cout << "No MCEventWeight objects in this event.\n";
    return;
  }

  // Print keys from the first neutrino vertex in the first event
  const auto & mw = weights_handle->at(0);

  std::cout << "Found " << mw.fWeight.size() << " weight keys:\n";
  for (const auto & kv : mw.fWeight) {
    std::cout << "  " << kv.first << "  (Nuni=" << kv.second.size() << ")\n";
  }
}
