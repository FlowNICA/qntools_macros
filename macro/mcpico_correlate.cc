#include "QnDataFrame.hpp"

vector <vector<string>> u1Q1=
{
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
  {"tru_pionP_PLAIN", "psi_rp_PLAIN"},
  {"tru_pionM_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u2Q2{
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
  {"tru_pionP_PLAIN", "psi_rp_PLAIN"},
  {"tru_pionM_PLAIN", "psi_rp_PLAIN"},
};

void mcpico_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.}});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  ROOT::RDataFrame d( "tree", inputFiles.c_str() );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wUnity3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return 1; };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
  for (auto &corr:u1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
  }
  for (auto &corr:u2Q2)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wSumWu, wy, qn, qn);
  }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = TFile::Open(outputFile.c_str(), "RECREATE");
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  corrFile->Close();
}
