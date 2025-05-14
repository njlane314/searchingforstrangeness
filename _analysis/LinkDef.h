// LinkDef.h
#ifdef __CLING__

#pragma link C++ namespace AnalysisFramework;
#pragma link C++ namespace AnalysisFramework::DataLoader; // If you want to use its functions directly
#pragma link C++ namespace AnalysisFramework::Selections; // If you want to use its const strings directly

#pragma link C++ enum AnalysisFramework::NTupleFileType;

#pragma link C++ struct AnalysisFramework::FilePathConfig+;
#pragma link C++ struct AnalysisFramework::DatasetInfo+;
#pragma link C++ struct AnalysisFramework::RunConfig+;
#pragma link C++ struct AnalysisFramework::DetVarDatasetInfo+;
#pragma link C++ struct AnalysisFramework::DetVarRunConfig+;
#pragma link C++ struct AnalysisFramework::CategoryDisplayInfo+; // From EventCategories.h

#pragma link C++ class AnalysisFramework::VariableManager+;
#pragma link C++ class AnalysisFramework::ConfigurationManager+;
#pragma link C++ class AnalysisFramework::EventCategories+;
#pragma link C++ class AnalysisFramework::RDataFrameExporter+; // Your renamed class

// If you have global functions in DataLoader namespace you want to call from ROOT prompt:
// #pragma link C++ function AnalysisFramework::DataLoader::LoadSample; // Example, might be complex due to RDataFrame return

// For std::vector and std::map of your custom types if used as function arguments or return types directly in ROOT
#pragma link C++ class std::vector<AnalysisFramework::RunConfig>+;
#pragma link C++ class std::map<std::string, AnalysisFramework::DatasetInfo>+;
// Add other std::vector/map instantiations if ROOT complains about them for your public interfaces

#endif