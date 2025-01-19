//
// Created by Misha on 3/7/2023.
//

using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PxPyPzE4D<double>>;

void unigen_correct(std::string list, std::string cm_energy="2.5", bool isCms=true)
{
  const double sNN = std::stod( cm_energy ); // in GeV
  const double M = 0.938; // in GeV/c^2
  const double T = sNN*sNN/(2.*M) - 2.*M;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );

  std::cout << "sqrtSnn = " << sNN << " GeV; T = " << T << "A GeV; Y_BEAM = " << Y_BEAM << std::endl;

  // Centrality for 2.5 GeV
  std::vector<float> cent_bins;
  std::vector<float> cent_b;

  // Xe+W @ 2.87 GeV
  // cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 95.};
  // cent_b = { 1.43172, 2.87406, 4.05374, 5.0333, 5.86304, 6.58249, 7.22197, 7.80409, 8.34523, 8.8571, 9.34824, 9.8255, 10.2956, 10.7666, 11.2494, 11.7595, 12.3179, 12.9533, 13.7034};
  // Xe+Xe @ 2.87 GeV
  cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 95.};
  cent_b = { 1.39543, 2.70341, 3.76519, 4.64771, 5.40276, 6.06901, 6.67396, 7.23605, 7.76659, 8.27183, 8.75493, 9.21803, 9.66421, 10.0995, 10.5351, 10.9889, 11.488, 12.0706, 12.7879};

  TStopwatch timer;
  timer.Start();
  
  std::string treename = "events";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
    .Alias( "b", "fB" )
    .Filter("b<20.")
    .Define("psi_rp", [](double _psi){ return (float)_psi; }, {"fPhi"})
    .Define("cent", [cent_bins,cent_b](float b){
      float cent = -1.;
      if (cent_b.size() == 0) return (float)-1.;
      if (cent_bins.size() == 0) return (float)-1.;
      if (b <= cent_b.at(0))
        cent = cent_bins.at(0);
      for (int i=1; i<cent_b.size(); ++i){
        if (b < cent_b.at(i) && b >= cent_b.at(i-1))
          cent = cent_bins.at(i);
      }
      return cent;
    }, {"b"})
    .Define("particles", [](RVec<double> &px, RVec<double> &py, RVec<double> &pz, RVec<double> &e){
      vector<fourVector> pv;
      int Np = px.size();
      for (int i=0; i<Np; i++){
        pv.push_back( {px.at(i), py.at(i), pz.at(i), e.at(i)} );
      }
      return pv;
    }, {"event.fParticles.fPx", "event.fParticles.fPy", "event.fParticles.fPz", "event.fParticles.fE"})
    .Define("pdg", [](RVec<int> &_pdg){ vector<int> pdg; for (auto &val:_pdg){ pdg.push_back((int)val); } return pdg; }, {"event.fParticles.fPdg"})
    .Define("pT", [](RVec<fourVector> &p){ vector<float> pt; for (auto &pi:p){ pt.push_back((float)pi.Pt()); } return pt; }, {"particles"})
    .Define("y", [Y_BEAM,isCms](RVec<fourVector> &p){
      vector<float> y;
      for (auto &pi:p){
        if (isCms) y.push_back((float)pi.Rapidity());
        else y.push_back((float)(pi.Rapidity() - Y_BEAM));
      }
      return y;
    }, {"particles"})
    .Define("phi", [](RVec<fourVector> &p){ vector<float> phi; for (auto &pi:p){ phi.push_back((float)pi.Phi()); } return phi; }, {"particles"})
  ;

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("b|psi_rp|cent"));
  correction_task.SetTrackVariables({std::regex("pT|y|phi|pdg")});

  correction_task.InitVariables();
  correction_task.AddEventAxis( { "cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.} } );

  VectorConfig psi_rp( "psi_rp", "psi_rp", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  std::vector<Qn::AxisD> sim_proton_axes{
          { "sim_y", 15, -1.5, 1.5 },
          { "sim_pT", 20, 0.0, 2.0 },
  };

  std::vector<Qn::AxisD> sim_pion_axes{
          { "sim_y", 15, -1.5, 1.5 },
          { "sim_pT", 15, 0.0, 1.5 },
  };

  VectorConfig tru_proton( "tru_proton", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_proton.SetHarmonicArray( {1, 2} );
  tru_proton.SetCorrections( {CORRECTION::PLAIN } );
  tru_proton.SetCorrectionAxes( sim_proton_axes );
  tru_proton.AddCut( "pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
  }, "tru_proton cut" );
  tru_proton.AddHisto2D({{"y", 300, -1.5, 1.5}, {"pT", 100, 0.0, 2.0}}, "is_proton");
  correction_task.AddVector(tru_proton);

  VectorConfig tru_pionP( "tru_pionP", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_pionP.SetHarmonicArray( {1, 2} );
  tru_pionP.SetCorrections( {CORRECTION::PLAIN } );
  tru_pionP.SetCorrectionAxes( sim_pion_axes );
  tru_pionP.AddCut( "pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 211;
  }, "tru_pionP cut" );
  tru_pionP.AddHisto2D({{"y", 300, -1.5, 1.5}, {"pT", 100, 0.0, 2.0}}, "is_pionP");
  correction_task.AddVector(tru_pionP);

  VectorConfig tru_pionM( "tru_pionM", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_pionM.SetHarmonicArray( {1, 2} );
  tru_pionM.SetCorrections( {CORRECTION::PLAIN } );
  tru_pionM.SetCorrectionAxes( sim_pion_axes );
  tru_pionM.AddCut( "pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -211;
  }, "tru_pionM cut" );
  tru_pionM.AddHisto2D({{"y", 300, -1.5, 1.5}, {"pT", 100, 0.0, 2.0}}, "is_pionM");
  correction_task.AddVector(tru_pionM);

  correction_task.Run();
}
