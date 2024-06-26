#include "QnDataFrame.hpp"

vector <vector<string>> Q1Q1=
{
  {"F1_RESCALED", "F2_RESCALED"},
  {"F1_RESCALED", "F3_RESCALED"},
  {"F2_RESCALED", "F3_RESCALED"},

  {"Tp_RESCALED", "F1_RESCALED"},
  {"Tp_RESCALED", "F2_RESCALED"},
  {"Tp_RESCALED", "F3_RESCALED"},

  {"Tneg_RESCALED", "F1_RESCALED"},
  {"Tneg_RESCALED", "F2_RESCALED"},
  {"Tneg_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> u1Q1=
{
  {"lambda_signal_RECENTERED",      "psi_rp_PLAIN"},
  {"lambda_signal_PLAIN",           "psi_rp_PLAIN"},
  {"lambda_background_RECENTERED",  "psi_rp_PLAIN"},
  {"lambda_background_PLAIN",       "psi_rp_PLAIN"},
  {"lambda_good_RECENTERED",        "psi_rp_PLAIN"},
  {"lambda_good_PLAIN",             "psi_rp_PLAIN"},
  {"lambda_true_PLAIN",             "psi_rp_PLAIN"},

  {"lambda_signal_RECENTERED",      "F1_RESCALED"},
  {"lambda_background_RECENTERED",  "F1_RESCALED"},
  {"lambda_good_RECENTERED",        "F1_RESCALED"},
  {"lambda_true_PLAIN",        "F1_RESCALED"},

  {"lambda_signal_RECENTERED",      "F2_RESCALED"},
  {"lambda_background_RECENTERED",  "F2_RESCALED"},
  {"lambda_good_RECENTERED",        "F2_RESCALED"},
  {"lambda_true_PLAIN",        "F2_RESCALED"},

  {"lambda_signal_RECENTERED",      "F3_RESCALED"},
  {"lambda_background_RECENTERED",  "F3_RESCALED"},
  {"lambda_good_RECENTERED",        "F3_RESCALED"},
  {"lambda_true_PLAIN",        "F3_RESCALED"},
};

void lambda_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 8, 0, 40});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  ROOT::RDataFrame d( "tree", inputFiles.c_str() );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
  for (auto &corr:Q1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);
  }
  for (auto &corr:u1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
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
