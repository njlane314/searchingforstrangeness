"""Utility to load plugin configuration for analysis and plotting.

This module reads a JSON configuration file that defines analysis and plot
plugins. Each entry is keyed by the plugin name and its configuration is the
value, without referencing concrete plugin classes. An example configuration
is::

{
  "analysis_plugins": [
    {"variables": [...]},
    {"regions": [...]}
  ],
  "plot_plugins": [
    {
      "histograms": [
        {"type": "stacked_histogram", "variable": "v", "region": "r"}
      ]
    }
  ]
}
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, Any, List

@dataclass
class PluginConfig:
    name: str
    config: Any


def load_plugins(path: str) -> Dict[str, Dict[str, PluginConfig]]:
    """Load analysis and plot plugins from a JSON file.

    Parameters
    ----------
    path: str
        Path to the JSON configuration.

    Returns
    -------
    Dict[str, Dict[str, PluginConfig]]
        Dictionary with two keys, ``analysis`` and ``plot`` mapping to dictionaries
        of plugin name to :class:`PluginConfig` instances.
    """
    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)

    def _parse(entries: List[Dict[str, Any]]) -> Dict[str, PluginConfig]:
        parsed: Dict[str, PluginConfig] = {}
        for entry in entries:
            if len(entry) != 1:
                raise ValueError("Each plugin entry must have exactly one key")
            name, config = next(iter(entry.items()))
            parsed[name] = PluginConfig(name, config)
        return parsed

    analysis_cfgs = _parse(raw.get("analysis_plugins", []))
    plot_cfgs = _parse(raw.get("plot_plugins", []))

    return {"analysis": analysis_cfgs, "plot": plot_cfgs}


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Load plugin configuration")
    parser.add_argument("config", help="Path to plugin JSON configuration")
    args = parser.parse_args()

    cfgs = load_plugins(args.config)
    for section, plugins in cfgs.items():
        print(f"{section} plugins:")
        for plugin in plugins.values():
            print(f"  - {plugin.name}: {json.dumps(plugin.config)}")


if __name__ == "__main__":
    main()
