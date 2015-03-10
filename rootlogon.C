{
  // include path for RooFit
  TString rfitpath("/software/tier3/osg/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms4/include");
  TString path = gSystem->GetIncludePath();
  path += "-I. -I$ROOTSYS/src -I";
  path += rfitpath;
  gSystem->SetIncludePath(path.Data());

  // for plots
  gROOT->Macro("CPlot.cc++");
  gROOT->Macro("KStyle.cc++");
}
