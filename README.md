# **searchingforstrangeness**

_This [searchingforstrangeness](https://github.com/njlane314/searchingforstrangeness) code is a [Pandora](https://github.com/PandoraPFA/larpandora) neutrino selection framework, built with [LArSoft](https://github.com/LArSoft/larsoft), using the design implemented in [searchingfornues](https://github.com/ubneutrinos/searchingfornues)._

_The framework loads a candidate neutrino slice given by Pandora, runs a single selection tool for each event, and a series of analysis tools for each slice. The selection and analysis tools are configurable at run, and the output is to a single root file. This implementation is specifically designed to search for rare neutrino processes that give distinct topologies or identifiable pattern signatures within the MicroBooNE detector; these signatures are identified using a visual deep-learning network that processes each of the wire planes independently._



lar -c eventdump.fcl -s /pnfs/uboone/overlay/uboone/reconstructed/prod_v08_00_00_82/prod_strange_resample_fhc/run2_fhc_reco2/00/01/18/45/PhysicsRun-2017_6_26_2_28_39-0011845-00142_20190227T070705_ext_unbiased_3_20190227091437_merged_20240806T065608_simmxd_detsim_mix_r1a_r1b_r1c_20240807_5e545787-fcd5-46d6-ac4e-5595c3305bed.root | grep POT