#include "QnDataFrame.hpp"

vector <vector<string>> Q1Q1{
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

vector <vector<string>> u1Q1{
  {"proton_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED"},

  {"pi_pos_RESCALED", "F1_RESCALED"},
  {"pi_pos_RESCALED", "F2_RESCALED"},
  {"pi_pos_RESCALED", "F3_RESCALED"},

  {"pi_neg_RESCALED", "F1_RESCALED"},
  {"pi_neg_RESCALED", "F2_RESCALED"},
  {"pi_neg_RESCALED", "F3_RESCALED"},

  {"proton_RESCALED", "psi_rp_PLAIN"},
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},

  {"pi_pos_RESCALED", "psi_rp_PLAIN"},
  {"tru_pi_pos_PLAIN", "psi_rp_PLAIN"},

  {"pi_neg_RESCALED", "psi_rp_PLAIN"},
  {"tru_pi_neg_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u2Q1Q1{
  {"proton_RESCALED", "F1_RESCALED", "F3_RESCALED"},
  {"pi_pos_RESCALED", "F1_RESCALED", "F3_RESCALED"},
  {"pi_neg_RESCALED", "F1_RESCALED", "F3_RESCALED"},
  {"proton_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"pi_pos_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"pi_neg_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"pi_pos_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"pi_neg_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED", "F3_RESCALED"},
  {"pi_pos_RESCALED", "F3_RESCALED", "F3_RESCALED"},
  {"pi_neg_RESCALED", "F3_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> u2Q2{
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_pos_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_neg_PLAIN", "psi_rp_PLAIN"},
  {"proton_RESCALED", "psi_rp_PLAIN"},
  {"pi_pos_RESCALED", "psi_rp_PLAIN"},
  {"pi_neg_RESCALED", "psi_rp_PLAIN"},
};

vector <vector<string>> u3Q1Q1Q1{
  {"proton_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"pi_pos_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"pi_neg_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"pi_pos_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"pi_neg_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},
  {"pi_pos_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},
  {"pi_neg_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},
};

vector< vector<string>> u3Q3{
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_pos_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_neg_PLAIN", "psi_rp_PLAIN"},
  {"proton_RESCALED", "psi_rp_PLAIN"},
  {"pi_pos_RESCALED", "psi_rp_PLAIN"},
  {"pi_neg_RESCALED", "psi_rp_PLAIN"},
};

vector <vector<string>> u4Q1Q1Q1Q1{
  {"proton_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"pi_pos_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"pi_neg_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"pi_pos_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"pi_neg_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},
  {"pi_pos_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},
  {"pi_neg_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},
};

vector< vector<string>> u4Q4{
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_pos_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_neg_PLAIN", "psi_rp_PLAIN"},
  {"proton_RESCALED", "psi_rp_PLAIN"},
  {"pi_pos_RESCALED", "psi_rp_PLAIN"},
  {"pi_neg_RESCALED", "psi_rp_PLAIN"},
};

void mpd_fixed_target_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  //Qn::AxisD centAxis({"b_norm", 20, 0, 2});
  Qn::AxisD centAxis({"cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.}});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  ROOT::RDataFrame d( "tree", inputFiles.c_str() );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  namespace P4 = Qn::Correlation::FourParticle;
  namespace P5 = Qn::Correlation::FiveParticle;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wUnity3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return 1; };
  auto wSumWu4part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) { return a.sumweights(); };
  auto wUnity4part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) { return 1; };
  auto wSumWu5part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d, const Qn::QVector &k) { return a.sumweights(); };
  auto wUnity5part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d, const Qn::QVector &k) { return 1; };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
  for (auto &corr:Q1Q1){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);
  }
  
  for (auto &corr:u1Q1){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
  }

  for (auto &corr:u2Q2){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wSumWu, wy, qn, qn);
  }

  for (auto &corr : u2Q1Q1 ){
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1y1", P3::yxy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1x1", P3::yyx(2, 1, 1), wSumWu3part, wy, qn, qn);
  }

  for (auto &corr:u3Q3){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x3", P2::xx(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y3", P2::yy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y3", P2::xy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x3", P2::yx(3, 3), wSumWu, wy, qn, qn);
  }

  for (auto &corr : u3Q1Q1Q1 ){
    std::array<std::string, 4> qn{corr.at(0), corr.at(1), corr.at(2), corr.at(3)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2)+"."+corr.at(3);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1x1x1", P4::xxxx(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1y1y1", P4::xxyy(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x1x1y1", P4::yxxy(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x1y1x1", P4::yxyx(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1x1y1", P4::xyxy(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1y1x1", P4::xyyx(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1x1x1", P4::yyxx(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1y1y1", P4::yyyy(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
  }

  for (auto &corr:u4Q4){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4x4", P2::xx(4, 4), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4y4", P2::yy(4, 4), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4y4", P2::xy(4, 4), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4x4", P2::yx(4, 4), wSumWu, wy, qn, qn);
  }

  for (auto &corr : u4Q1Q1Q1Q1 ){
    std::array<std::string, 5> qn{corr.at(0), corr.at(1), corr.at(2), corr.at(3), corr.at(4)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2)+"."+corr.at(3)+"."+corr.at(4);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4x1x1x1x1", P5::xxxxx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4x1x1y1y1", P5::xxxyy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4x1y1x1y1", P5::xxyxy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4x1y1y1x1", P5::xxyyx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4y1x1x1y1", P5::xyxxy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4y1x1y1x1", P5::xyxyx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4y1y1x1x1", P5::xyyxx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4y1y1y1y1", P5::xyyyy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4x1x1x1y1", P5::yxxxy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4x1x1y1x1", P5::yxxyx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4x1y1x1x1", P5::yxyxx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4x1y1y1y1", P5::yxyyy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4y1x1x1x1", P5::yyxxx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4y1x1y1y1", P5::yyxyy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4y1y1x1y1", P5::yyyxy(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4y1y1y1x1", P5::yyyyx(4, 1, 1, 1, 1), wSumWu5part, wy, qn, qn);
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
