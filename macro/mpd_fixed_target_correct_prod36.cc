//
// Created by Misha on 3/7/2023.
//

using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

void mpd_fixed_target_correct_prod36(std::string list, std::string cm_energy="2.5", std::string eff_file="", std::string pid_file=""){

  const double sNN = std::stod( cm_energy ); // in GeV
  const double M = 0.938; // in GeV/c^2
  const double T = sNN*sNN/(2.*M) - 2.*M;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );

  std::cout << "sqrtSnn = " << sNN << " GeV; T = " << T << "A GeV; Y_BEAM = " << Y_BEAM << std::endl;

  if (eff_file != "") std::cout << "Efficiency file: " << eff_file.c_str() << std::endl;
  std::unique_ptr<TFile> fiEff{TFile::Open( eff_file.c_str(), "read")};

  TH2D* h_eff_proton{nullptr}; 
  TH2D* h_eff_pionP {nullptr}; 
  TH2D* h_eff_pionM {nullptr}; 
  TH2D* h_eff_pions {nullptr}; 
  TH2D* h_eff_kaonP {nullptr}; 
  TH2D* h_eff_kaonM {nullptr}; 
  TH2D* h_eff_kaons {nullptr}; 
  
  fiEff->GetObject("h2_effYPt_proton", h_eff_proton);
  fiEff->GetObject("h2_effYPt_pionP" , h_eff_pionP );
  fiEff->GetObject("h2_effYPt_pionM" , h_eff_pionM );
  fiEff->GetObject("h2_effYPt_pions" , h_eff_pions );
  fiEff->GetObject("h2_effYPt_kaonP" , h_eff_kaonP );
  fiEff->GetObject("h2_effYPt_kaonM" , h_eff_kaonM );
  fiEff->GetObject("h2_effYPt_kaons" , h_eff_kaons );
    
  if (!h_eff_proton) std::cerr << "Warning: No efficiency histogram h2_effYPt_proton was found in file " << eff_file.c_str() << "\n";
  if (!h_eff_pionP)  std::cerr << "Warning: No efficiency histogram h2_effYPt_pionP was found in file "  << eff_file.c_str() << "\n";
  if (!h_eff_pionM)  std::cerr << "Warning: No efficiency histogram h2_effYPt_pionM was found in file "  << eff_file.c_str() << "\n";
  if (!h_eff_pions)  std::cerr << "Warning: No efficiency histogram h2_effYPt_pions was found in file "  << eff_file.c_str() << "\n";
  if (!h_eff_kaonP)  std::cerr << "Warning: No efficiency histogram h2_effYPt_kaonP was found in file "  << eff_file.c_str() << "\n";
  if (!h_eff_kaonM)  std::cerr << "Warning: No efficiency histogram h2_effYPt_kaonM was found in file "  << eff_file.c_str() << "\n";
  if (!h_eff_kaons)  std::cerr << "Warning: No efficiency histogram h2_effYPt_kaons was found in file "  << eff_file.c_str() << "\n";

  if (pid_file != "") std::cout << "PID file: " << pid_file.c_str() << std::endl;
  std::unique_ptr<TFile> fiPid{TFile::Open( pid_file.c_str(), "read")};

  TF1 *f1_dedx_pip      {nullptr};
  TF1 *f1_dedx_p        {nullptr};
  TF1 *f1_dedx_sigm_pip {nullptr};
  TF1 *f1_dedx_sigm_p   {nullptr};
  TF1 *f1_m2_sigm_pip   {nullptr};
  TF1 *f1_m2_sigm_p     {nullptr};

  fiPid->GetObject("f1_dedx_pip",      f1_dedx_pip);
  fiPid->GetObject("f1_dedx_proton",   f1_dedx_p);
  fiPid->GetObject("f1_dedx_sigm_pip", f1_dedx_sigm_pip);
  fiPid->GetObject("f1_dedx_sigm_p",   f1_dedx_sigm_p);
  fiPid->GetObject("f1_m2_sigm_pip",   f1_m2_sigm_pip);
  fiPid->GetObject("f1_m2_sigm_p",     f1_m2_sigm_p);

  if (!f1_dedx_pip) std::cerr << "Warning: No PID function f1_dedx_pip was found in file " << pid_file.c_str() << "\n";
  if (!f1_dedx_p) std::cerr << "Warning: No PID function f1_dedx_proton was found in file " << pid_file.c_str() << "\n";
  if (!f1_dedx_sigm_pip) std::cerr << "Warning: No PID function f1_dedx_sigm_pip was found in file " << pid_file.c_str() << "\n";
  if (!f1_dedx_sigm_p) std::cerr << "Warning: No PID function f1_dedx_sigm_p was found in file " << pid_file.c_str() << "\n";
  if (!f1_m2_sigm_pip) std::cerr << "Warning: No PID function f1_m2_sigm_pip was found in file " << pid_file.c_str() << "\n";
  if (!f1_m2_sigm_p) std::cerr << "Warning: No PID function f1_m2_sigm_p was found in file " << pid_file.c_str() << "\n";
  
  std::vector<int> f1_modules = { 15, 16, 17, 22, 24, 29, 30, 31 };
  std::vector<int> f2_modules = { 7, 8, 9, 10, 11, 14, 18, 21, 25, 28, 32, 35, 36, 37, 38, 39 };
  std::vector<int> f3_modules = { 1, 2, 3, 4, 5, 6, 12, 13, 19, 20, 26, 27, 33, 34, 40, 41, 42, 43, 44, 45 };

  // Centrality for 2.5 GeV
  std::vector<float> cent_bins;
  std:;vector<int> cent_mult;
  std::vector<float> cent_b;
  // if (sNN == 2.5){
  //   cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 45., 55., 65., 75., 90.};
  //   cent_mult = {147,129,112,98,85,73,63,54,39,28,19,13,1};
  //   cent_b = {2.89,4.11,5.10,5.92,6.63,7.26,7.85,8.39,9.39,10.25,11.03,11.94,16.34};
  // }
  // if (sNN == 3.0){
  //   cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 45., 55., 65., 75., 90.};
  //   cent_mult = {182,159,139,121,106,92,80,69,50,36,24,16,1};;
  //   cent_b = {2.88,4.10,5.08,5.89,6.59,7.21,7.79,8.33,9.33,10.21,11.00,11.87,15.83};
  // }
  // if (sNN == 3.5){
  //   cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 45., 55., 65., 75., 90.};
  //   cent_mult = {219,191,166,145,127,110,95,82,60,42,29,18,1};
  //   cent_b = {2.84,4.04,4.99,5.79,6.47,7.08,7.65,8.19,9.17,10.03,10.82,11.69,15.76};
  // }
  // Xe+W - prod 35
  // cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 95.};
  // cent_mult = {202, 132, 112, 95, 81, 68, 57, 47, 39, 32, 26, 21, 16, 12, 9, 6, 4, 2, 1};
  // cent_b = { 1.43172, 2.87406, 4.05374, 5.0333, 5.86304, 6.58249, 7.22197, 7.80409, 8.34523, 8.8571, 9.34824, 9.8255, 10.2956, 10.7666, 11.2494, 11.7595, 12.3179, 12.9533, 13.7034};
  // Xe+Xe - prod 36
  cent_bins = {2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 95.};
  cent_mult = {162, 103, 87, 74, 63, 53, 44, 36, 30, 24, 19, 15, 12, 9, 7, 5, 3, 2, 1};
  cent_b = { 1.39543, 2.70341, 3.76519, 4.64771, 5.40276, 6.06901, 6.67396, 7.23605, 7.76659, 8.27183, 8.75493, 9.21803, 9.66421, 10.0995, 10.5351, 10.9889, 11.488, 12.0706, 12.7879};
  
  TStopwatch timer;
  timer.Start();
  
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
    .Filter("mcB < 20.") // at least one filter is mandatory!!!
    .Define("b", [](double b){ return (float)b;}, {"mcB"})
    //.Define( "b_norm", [NUCLEUS_RADIUS](double b){ return (float)b / (float)NUCLEUS_RADIUS; },{"mcB"})
    .Define( "refMult", [](vector<fourVector> _p, RVec<int> _nhits){
      int Mult = 0;
      for (int i=0; i<_p.size(); ++i) {
        auto mom = _p.at(i);
        auto nhits = _nhits.at(i);
        if (nhits<=16) continue;
        if(mom.Eta()>0. && mom.Eta()<2.)
          Mult++;
      }
      return Mult;
    },{"recoGlobalMom", "recoGlobalNhits"})
    .Define( "cent", [cent_mult,cent_bins](int refmult){
      float cent = -1.;
      if (cent_mult.size() == 0) return (float)-1.;
      if (cent_bins.size() == 0) return (float)-1.;
      if (refmult >= cent_mult.at(0))
        cent = cent_bins.at(0);
      for (int i=1; i<cent_mult.size(); ++i){
        if (refmult >= cent_mult.at(i) && refmult < cent_mult.at(i-1))
          cent = cent_bins.at(i);
      }
      return cent;
    },{"refMult"})
    .Define( "psi_rp", "mcRP")
    .Define( "sim_y", [Y_BEAM]( const RVec<int> vec_pdg, vector<fourVector> vec_momentum ){
      std::vector<float> vec_y;
      vec_y.reserve( vec_pdg.size() );
      for( int i=0; i<vec_pdg.size(); ++i ){
        auto pdg = vec_pdg.at(i);
        if( pdg == 0 ){
          vec_y.push_back( -999. );
          continue;
        }
        auto pz = vec_momentum.at(i).Pz();
        auto p = vec_momentum.at(i).P();
        TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
        if (!particle) {
          vec_y.push_back( -999. );
          continue;
        }
        auto m = particle->Mass();
        auto E = sqrt( p*p + m*m );
        auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
        vec_y.push_back(y);
      }
      return vec_y;
    },{"simPdg", "simMom"})
    .Define( "sim_pT", "std::vector<float> vec_pT; for( auto mom : simMom ){ vec_pT.push_back( mom.Pt() ); } return vec_pT;" )
    .Define( "sim_phi", "std::vector<float> vec_phi; for( auto mom : simMom ){ vec_phi.push_back( mom.Phi() ); } return vec_phi;" )
    .Define( "sim_pdg", "std::vector<int> vec_pdg; for( auto pdg : simPdg ){ vec_pdg.push_back( pdg ); } return vec_pdg;" )
    .Define( "sim_is_proton", "simPdg==2212" )
    .Define( "sim_is_pionP", "simPdg==211" )
    .Define( "sim_is_pionM", "simPdg==-211" )
    .Define( "sim_is_pions", "abs(simPdg)==211" )
    .Define( "sim_is_kaonP", "simPdg==321" )
    .Define( "sim_is_kaonM", "simPdg==-321" )
    .Define( "sim_is_kaons", "abs(simPdg)==321" )
    .Define( "sim_mother_id", "std::vector<int> vec_mother_id; for( auto id : simMotherId ){ vec_mother_id.push_back( id ); } return vec_mother_id;" )
    .Define( "tr_pT", "std::vector<float> vec_pT; for( auto mom : recoGlobalMom ){ vec_pT.push_back( mom.Pt() ) ;} return vec_pT;" )
    .Define( "tr_phi", "std::vector<float> vec_phi; for( auto mom : recoGlobalMom ){ vec_phi.push_back( mom.Phi() ); } return vec_phi;" )
    .Define( "tr_eta", "std::vector<float> vec_eta; for( auto mom : recoGlobalMom ){ vec_eta.push_back( mom.Eta() ); } return vec_eta;" )
    .Define( "tr_charge", "std::vector<short> vec_charge; for( auto q : recoGlobalCharge ){ vec_charge.push_back( q ); } return vec_charge;" )
    .Define( "tr_pq", [](vector<fourVector> _p,RVec<short> _q){
      RVec <float> pq_;
      for (int i=0; i<_p.size(); ++i){
        auto mom = _p.at(i);
        auto charge = _q.at(i);
        if (charge != 0)
          pq_.push_back(mom.P()/(float)(charge));
        else
          pq_.push_back(-9999.);
      }
      return pq_;
    }, {"recoGlobalMom", "recoGlobalCharge"})
    .Define( "tr_pdg", "recoGlobalSimPdg" )
    .Define( "tr_TpcDedx", "recoGlobalTpcDedx")
    .Define( "tr_TofMass2", "recoGlobalTofMass2")
    //.Define( "tr_primary", "recoGlobalSimMotherId==-1" )
    .Define( "tr_y", [Y_BEAM]( const RVec<int> vec_pdg, vector<fourVector> vec_momentum ){
      std::vector<float> vec_y;
      vec_y.reserve( vec_pdg.size() );
      for( int i=0; i<vec_pdg.size(); ++i ){
        auto pdg = vec_pdg.at(i);
        if( pdg == 0 ){
          vec_y.push_back( -999. );
          continue;
        }
        auto pz = vec_momentum.at(i).Pz();
        auto p = vec_momentum.at(i).P();
        TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
        if (!particle) {
          vec_y.push_back( -999. );
          continue;
        }
        auto m = particle->Mass();
        auto E = sqrt( p*p + m*m );
        auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
        vec_y.push_back(y);
      }
      return vec_y;
    },{"recoGlobalSimPdg", "recoGlobalMom"} )
    // .Define( "tr_is_proton", "recoGlobalSimPdg==2212" )
    // .Define( "tr_is_pionP", "recoGlobalSimPdg==211" )
    // .Define( "tr_is_pionM", "recoGlobalSimPdg==-211" )
    // .Define( "tr_is_pions", "abs(recoGlobalSimPdg)==211" )
    .Define( "tr_is_kaonP", "recoGlobalSimPdg==321" )
    .Define( "tr_is_kaonM", "recoGlobalSimPdg==-321" )
    .Define( "tr_is_kaons", "abs(recoGlobalSimPdg)==321" )
    .Define( "tr_nhits", "std::vector<int> vec_nhits; for( auto n : recoGlobalNhits ){ vec_nhits.push_back( n ); } return vec_nhits;" )
    .Define( "tr_dca", "vector <float> vec_dca; for ( auto dca : recoGlobalDca ){ vec_dca.push_back( sqrt( dca.Mag2() ) ); } return vec_dca;" )
    .Define( "tr_primary", "RVec <int> vec_prim; for ( auto dca : tr_dca ){ int prim = ( dca < 1. ) ? 1 : 0; vec_prim.push_back( prim ); } return vec_prim;")
    .Define("isTOF", "recoGlobalTofFlag")
    .Define("tr_nsigDedxProton", [f1_dedx_p, f1_dedx_sigm_p]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = vec_Pq.at(i);
        auto dedx = vec_dedx.at(i);
        if (pq<0.3)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_p->Eval(pq))/(f1_dedx_p->Eval(pq)*f1_dedx_sigm_p->Eval(pq)) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TpcDedx"})
    .Define("tr_nsigM2Proton", [f1_m2_sigm_p]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = vec_Pq.at(i);
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(2212)->Mass(), 2);
        if (pq<0.3 || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_p->Eval(pq) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TofMass2", "isTOF"})
    .Define("tr_nsigDedxPionP", [f1_dedx_pip, f1_dedx_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = vec_Pq.at(i);
        auto dedx = vec_dedx.at(i);
        if (abs(pq)<0.1 || pq<0.)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_pip->Eval(pq))/(f1_dedx_pip->Eval(pq)*f1_dedx_sigm_pip->Eval(pq)) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TpcDedx"})
    .Define("tr_nsigM2PionP", [f1_m2_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = vec_Pq.at(i);
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2);
        if (abs(pq)<0.1 || pq<0. || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_pip->Eval(pq) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TofMass2", "isTOF"})
    .Define("tr_nsigDedxPionM", [f1_dedx_pip, f1_dedx_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = (float)-1.*vec_Pq.at(i);
        auto dedx = vec_dedx.at(i);
        if (abs(pq)<0.1 || pq<0.)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_pip->Eval(pq))/(f1_dedx_pip->Eval(pq)*f1_dedx_sigm_pip->Eval(pq)) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TpcDedx"})
    .Define("tr_nsigM2PionM", [f1_m2_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = (float)-1.*vec_Pq.at(i);
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2);
        if (abs(pq)<0.1 || pq<0. || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_pip->Eval(pq) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TofMass2", "isTOF"})
    .Define("tr_nsigDedxPions", [f1_dedx_pip, f1_dedx_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_dedx ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_dedx.size() );
      for( int i=0; i<vec_dedx.size(); ++i){
        auto pq = abs(vec_Pq.at(i));
        auto dedx = vec_dedx.at(i);
        if (pq<0.1)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (dedx - f1_dedx_pip->Eval(pq))/(f1_dedx_pip->Eval(pq)*f1_dedx_sigm_pip->Eval(pq)) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TpcDedx"})
    .Define("tr_nsigM2Pions", [f1_m2_sigm_pip]( RVec<float> vec_Pq, RVec<float> vec_m2, RVec<bool> vec_tof ){
      RVec<float> vec_nsig;
      vec_nsig.reserve( vec_m2.size() );
      for( int i=0; i<vec_m2.size(); ++i){
        auto tof = vec_tof.at(i);
        auto pq = abs(vec_Pq.at(i));
        auto m2 = vec_m2.at(i);
        auto m2_p = pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2);
        if (pq<0.1 || !tof)
          vec_nsig.push_back( -999. );
        else
          vec_nsig.push_back( (m2 - m2_p)/f1_m2_sigm_pip->Eval(pq) );
      }
      return vec_nsig;
    },{"tr_pq", "tr_TofMass2", "isTOF"})
    .Define("tr_is_proton", []( RVec<float> vec_nsigDedx_p, RVec<float> vec_nsigM2_p, RVec<float> vec_nsigDedx_pi, RVec<float> vec_nsigM2_pi ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx_p.size() );
      for( int i=0; i<vec_nsigDedx_p.size(); ++i){
        auto dedx_p = vec_nsigDedx_p.at(i);
        auto m2_p = vec_nsigM2_p.at(i);
        auto dedx_pi = vec_nsigDedx_pi.at(i);
        auto m2_pi = vec_nsigM2_pi.at(i);
        if ( sqrt(pow(dedx_p,2)+pow(m2_p,2)) < 2. &&
             sqrt(pow(dedx_pi,2)+pow(m2_pi,2)) > 3.)
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"tr_nsigDedxProton", "tr_nsigM2Proton", "tr_nsigDedxPionP", "tr_nsigM2PionP"})
    .Define("tr_is_pionP", []( RVec<float> vec_nsigDedx_p, RVec<float> vec_nsigM2_p, RVec<float> vec_nsigDedx_pi, RVec<float> vec_nsigM2_pi ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx_p.size() );
      for( int i=0; i<vec_nsigDedx_p.size(); ++i){
        auto dedx_p = vec_nsigDedx_p.at(i);
        auto m2_p = vec_nsigM2_p.at(i);
        auto dedx_pi = vec_nsigDedx_pi.at(i);
        auto m2_pi = vec_nsigM2_pi.at(i);
        if ( sqrt(pow(dedx_p,2)+pow(m2_p,2)) > 3. &&
             sqrt(pow(dedx_pi,2)+pow(m2_pi,2)) < 2.)
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"tr_nsigDedxProton", "tr_nsigM2Proton", "tr_nsigDedxPionP", "tr_nsigM2PionP"})
    .Define("tr_is_pionM", []( RVec<float> vec_nsigDedx, RVec<float> vec_nsigM2 ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx.size() );
      for( int i=0; i<vec_nsigDedx.size(); ++i){
        auto dedx = vec_nsigDedx.at(i);
        auto m2 = vec_nsigM2.at(i);
        if ( sqrt(pow(dedx,2)+pow(m2,2)) < 2. )
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"tr_nsigDedxPionM", "tr_nsigM2PionM"})
    .Define("tr_is_pions", []( RVec<float> vec_nsigDedx_p, RVec<float> vec_nsigM2_p, RVec<float> vec_nsigDedx_pi, RVec<float> vec_nsigM2_pi ){
      RVec<int> vec_pid;
      vec_pid.reserve( vec_nsigDedx_p.size() );
      for( int i=0; i<vec_nsigDedx_p.size(); ++i){
        auto dedx_p = vec_nsigDedx_p.at(i);
        auto m2_p = vec_nsigM2_p.at(i);
        auto dedx_pi = vec_nsigDedx_pi.at(i);
        auto m2_pi = vec_nsigM2_pi.at(i);
        if ( sqrt(pow(dedx_p,2)+pow(m2_p,2)) > 3. &&
             sqrt(pow(dedx_pi,2)+pow(m2_pi,2)) < 2.)
          vec_pid.push_back(1);
        else
          vec_pid.push_back(0);
      }
      return vec_pid;
    },{"tr_nsigDedxProton", "tr_nsigM2Proton", "tr_nsigDedxPions", "tr_nsigM2Pions"})
    .Define( "fhcal_module_id", "fhcalModId" )
    .Define( "fhcal_module_pos", "fhcalModPos" )
    .Define( "fhcal_module_x", "std::vector<float> vec_x{}; for( auto pos : fhcal_module_pos ){ vec_x.push_back( pos.X() ); }; return vec_x;" )
    .Define( "fhcal_module_y", "std::vector<float> vec_y{}; for( auto pos : fhcal_module_pos ){ vec_y.push_back( pos.Y() ); }; return vec_y;" )
    .Define( "fhcal_module_z", "std::vector<float> vec_z{}; for( auto pos : fhcal_module_pos ){ vec_z.push_back( pos.Z() ); }; return vec_z;" )
    .Define( "fhcal_module_phi", "std::vector<float> vec_phi{}; for( auto pos : fhcal_module_pos ){ vec_phi.push_back( pos.Phi() ); }; return vec_phi;" )
    .Define( "fhcal_module_energy", "std::vector<float> vec_en; for( auto en : fhcalModE){ vec_en.push_back( en ); } return vec_en;" )
    .Define( "tr_WeightProton", [h_eff_proton](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = _pdg.at(i);
        auto prim = std::round(_primary.at(i));
        auto eff = (h_eff_proton) ? (float)h_eff_proton->GetBinContent(h_eff_proton->GetXaxis()->FindBin(y), h_eff_proton->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_proton", "tr_primary"})
    .Define( "tr_WeightPionP", [h_eff_pionP](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = std::round(_pdg.at(i));
        auto prim = std::round(_primary.at(i));
        auto eff =(h_eff_pionP) ? (float)h_eff_pionP->GetBinContent(h_eff_pionP->GetXaxis()->FindBin(y), h_eff_pionP->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_pionP", "tr_primary"})
    .Define( "tr_WeightPionM", [h_eff_pionM](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = std::round(_pdg.at(i));
        auto prim = std::round(_primary.at(i));
        auto eff = (h_eff_pionM) ? (float)h_eff_pionM->GetBinContent(h_eff_pionM->GetXaxis()->FindBin(y), h_eff_pionM->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_pionM", "tr_primary"})
    .Define( "tr_WeightPions", [h_eff_pions](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = std::round(_pdg.at(i));
        auto prim = std::round(_primary.at(i));
        auto eff = (h_eff_pions) ? (float)h_eff_pions->GetBinContent(h_eff_pions->GetXaxis()->FindBin(y), h_eff_pions->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_pions", "tr_primary"})
    .Define( "tr_WeightKaonP", [h_eff_kaonP](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = std::round(_pdg.at(i));
        auto prim = std::round(_primary.at(i));
        auto eff = (h_eff_kaonP) ? (float)h_eff_kaonP->GetBinContent(h_eff_kaonP->GetXaxis()->FindBin(y), h_eff_kaonP->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_kaonP", "tr_primary"})
    .Define( "tr_WeightKaonM", [h_eff_kaonM](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = std::round(_pdg.at(i));
        auto prim = std::round(_primary.at(i));
        auto eff = (h_eff_kaonM) ? (float)h_eff_kaonM->GetBinContent(h_eff_kaonM->GetXaxis()->FindBin(y), h_eff_kaonM->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_kaonM", "tr_primary"})
    .Define( "tr_WeightKaons", [h_eff_kaons](std::vector<float> _pt, std::vector<float> _y, RVec<int> _pdg, RVec<int> _primary){
      std::vector<float> vec_w;
      for (int i=0; i<_pdg.size(); ++i){
        auto pt = _pt.at(i);
        auto y = _y.at(i);
        auto pdg = std::round(_pdg.at(i));
        auto prim = std::round(_primary.at(i));
        auto eff = (h_eff_kaons) ? (float)h_eff_kaons->GetBinContent(h_eff_kaons->GetXaxis()->FindBin(y), h_eff_kaons->GetYaxis()->FindBin(pt)) : (float)1;
        if (eff > 1e-3 && pdg && prim)
          vec_w.push_back(1./eff);
        else
          vec_w.push_back(0.);
      }
      return vec_w;
    },{"tr_pT", "tr_y", "tr_is_kaons", "tr_primary"})
    ;
  
  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("b_norm|psi_rp|cent"));
  correction_task.SetChannelVariables({std::regex("fhcal_module_(id|phi|energy|x|y|z)")});
  correction_task.SetTrackVariables({
      std::regex("tr_(pT|pq|y|eta|phi|charge|pdg|primary|TpcDedx|TofMass2|is_proton|is_pionP|is_pionM|is_pions|nhits|dca|WeightProton|WeightPionP|WeightPionM|WeightPions|WeightKaonP|WeightKaonM|WeightKaons)"),
      std::regex("sim_(pT|y|phi|pdg|is_proton|is_pionP|is_pionM|is_pions|mother_id)")
    });

  correction_task.InitVariables();
  //correction_task.AddEventAxis( { "b_norm", 10, 0, 2.0 } );
  correction_task.AddEventAxis( { "cent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 50., 60., 70., 80.} } );

  VectorConfig f1( "F1", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f1.AddCut( "fhcal_module_id", [f1_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcal_module_x", 100, -100, 100}, {"fhcal_module_y", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f2.AddCut( "fhcal_module_id", [f2_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcal_module_x", 100, -100, 100}, {"fhcal_module_y", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f3.AddCut( "fhcal_module_id", [f3_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcal_module_x", 100, -100, 100}, {"fhcal_module_y", 100, -100, 100}});
  correction_task.AddVector(f3);

  std::vector<Qn::AxisD> proton_axes{
          { "tr_y", 15, -1.5, 1.5 },
          { "tr_pT", 20, 0.0, 2.0 },
  };
  std::vector<Qn::AxisD> pion_axes{
          { "tr_y", 15, -1.5, 1.5 },
          { "tr_pT", 15, 0.0, 1.5 },
  };

  VectorConfig proton( "proton", "tr_phi", "tr_WeightProton", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2} );
  proton.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  proton.SetCorrectionAxes( proton_axes );
  // proton.AddCut( "tr_pdg", [](double pid){
  //   auto pdg_code = std::round(pid);
  //   return pdg_code == 2212;
  // }, "proton cut" );
  proton.AddCut( "tr_is_proton", [](double _pid){
    auto pid = std::round(_pid);
    return pid == 1;
  }, "proton cut" );
  proton.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  proton.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "q > 0" );
  proton.AddCut( "tr_nhits", [](double nhits){ return nhits > 27; }, "Nhits > 27" );
  proton.AddCut( "tr_dca", [](double dca){ return dca < 1.; }, "dca < 1 cm" );
  proton.AddCut( "tr_eta", [](double eta){ return eta < 2.; }, "no tracks in FHCal");
  proton.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightProton");
  proton.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TofMass2", 220, -0.2, 2.0}}, "tr_is_proton");
  proton.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TpcDedx", 500, 0., 5.e3}}, "tr_is_proton");
  correction_task.AddVector(proton);
  
  VectorConfig pi_pos( "pi_pos", "tr_phi", "tr_WeightPionP", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  pi_pos.SetHarmonicArray( {1, 2} );
  pi_pos.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  pi_pos.SetCorrectionAxes( pion_axes );
  // pi_pos.AddCut( "tr_pdg", [](double pid){
  //   auto pdg_code = std::round(pid);
  //   return pdg_code == 211;
  // }, "pi_pos cut" );
  pi_pos.AddCut( "tr_is_pionP", [](double _pid){
    auto pid = std::round(_pid);
    return pid == 1;
  }, "pi_pos cut" );
  pi_pos.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  pi_pos.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "q > 0" );
  pi_pos.AddCut( "tr_nhits", [](double nhits){ return nhits > 22; }, "Nhits > 22" );
  pi_pos.AddCut( "tr_dca", [](double dca){ return dca < 3.5; }, "dca < 3.5 cm" );
  pi_pos.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightPionP");
  pi_pos.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TofMass2", 220, -0.2, 2.0}}, "tr_is_pionP");
  pi_pos.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TpcDedx", 500, 0., 5.e3}}, "tr_is_pionP");
  correction_task.AddVector(pi_pos);
  
  VectorConfig pi_neg( "pi_neg", "tr_phi", "tr_WeightPionM", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  pi_neg.SetHarmonicArray( {1, 2} );
  pi_neg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  pi_neg.SetCorrectionAxes( pion_axes );
  // pi_neg.AddCut( "tr_pdg", [](double pid){
  //   auto pdg_code = std::round(pid);
  //   return pdg_code == -211;
  // }, "pi_neg cut" );
  // pi_neg.AddCut( "tr_is_pionM", [](double _pid){
  //   auto pid = std::round(_pid);
  //   return pid == 1;
  // }, "pi_neg cut" );
  pi_neg.AddCut( "tr_charge", [](double _charge){
      auto charge = std::round(_charge);
      return charge < 0;
    }, "pi_neg cut");
  pi_neg.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  pi_neg.AddCut( "tr_charge", [](double charge){ return charge < 0; }, "q < 0" );
  pi_neg.AddCut( "tr_nhits", [](double nhits){ return nhits > 22; }, "Nhits > 22" );
  pi_neg.AddCut( "tr_dca", [](double dca){ return dca < 3.5; }, "dca < 3.5 cm" );
  pi_neg.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightPionM");
  pi_neg.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TofMass2", 220, -0.2, 2.0}}, "tr_is_pionM");
  pi_neg.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TpcDedx", 500, 0., 5.e3}}, "tr_is_pionM");
  correction_task.AddVector(pi_neg);
  
  VectorConfig Tp( "Tp", "tr_phi", "tr_WeightProton", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1, 2} );
  Tp.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  // Tp.AddCut( "tr_pdg", [](double pid){
  //   auto pdg_code = std::round(pid);
  //   return pdg_code == 2212;
  // }, "proton cut" );
  Tp.AddCut( "tr_is_proton", [](double _pid){
    auto pid = std::round(_pid);
    return pid == 1;
  }, "proton cut" );
  Tp.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  Tp.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "proton cut" );
  Tp.AddCut( "tr_y", [](double ycm){
    return -1. < ycm && ycm < -0.6;
  }, "Tp ycm cut" );
  Tp.AddCut( "tr_nhits", [](double nhits){ return nhits > 27; }, "Nhits > 27" );
  Tp.AddCut( "tr_dca", [](double dca){ return dca < 1.; }, "dca < 1 cm" );
  Tp.AddCut( "tr_eta", [](double eta){ return eta < 2.; }, "no tracks in FHCal");
  Tp.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_is_proton");
  Tp.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TofMass2", 220, -0.2, 2.0}}, "tr_is_proton");
  Tp.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TpcDedx", 500, 0., 5.e3}}, "tr_is_proton");
  correction_task.AddVector(Tp);
  
  VectorConfig Tneg( "Tneg", "tr_phi", "tr_WeightPions", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  // Tneg.AddCut( "tr_pdg", [](double pid){
  //   auto pdg_code = std::round(pid);
  //   return abs(pdg_code) == 211;
  // }, "pion cut" );
  // Tneg.AddCut( "tr_is_pionM", [](double _pid){
  //   auto pid = std::round(_pid);
  //   return pid == 1;
  // }, "pi_neg cut" );
  Tneg.AddCut( "tr_charge", [](double _charge){
      auto charge = std::round(_charge);
      return charge < 0;
  }, "pi_neg cut");
  Tneg.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  Tneg.AddCut( "tr_y", [](double eta){
    return -1.3 < eta && eta < -0.2;
  }, "Tneg y cut" );
  Tneg.AddCut( "tr_nhits", [](double nhits){ return nhits > 22; }, "Nhits > 22" );
  Tneg.AddCut( "tr_dca", [](double dca){ return dca < 3.5; }, "dca < 3.5 cm" );
  Tneg.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightPions");
  Tneg.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TofMass2", 220, -0.2, 2.0}}, "tr_is_pionM");
  Tneg.AddHisto2D({{"tr_pq", 1000, -5., 5.}, {"tr_TpcDedx", 500, 0., 5.e3}}, "tr_is_pionM");
  correction_task.AddVector(Tneg);

  std::vector<Qn::AxisD> sim_proton_axes{
          { "sim_y", 15, -1.5, 1.5 },
          { "sim_pT", 20, 0.0, 2.0 },
  };

  std::vector<Qn::AxisD> sim_pion_axes{
          { "sim_y", 15, -1.5, 1.5 },
          { "sim_pT", 15, 0.0, 1.5 },
  };

  VectorConfig tru_proton( "tru_proton", "sim_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_proton.SetHarmonicArray( {1, 2} );
  tru_proton.SetCorrections( {CORRECTION::PLAIN } );
  tru_proton.SetCorrectionAxes( sim_proton_axes );
  tru_proton.AddCut( "sim_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
  }, "tru_proton cut" );
  tru_proton.AddCut( "sim_mother_id", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -1;
  }, "cut on primary" );
  tru_proton.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "sim_is_proton");
  correction_task.AddVector(tru_proton);

  VectorConfig tru_pi_pos( "tru_pi_pos", "sim_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_pi_pos.SetHarmonicArray( {1, 2} );
  tru_pi_pos.SetCorrections( {CORRECTION::PLAIN } );
  tru_pi_pos.SetCorrectionAxes( sim_pion_axes );
  tru_pi_pos.AddCut( "sim_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 211;
  }, "tru_pi_pos cut" );
  tru_pi_pos.AddCut( "sim_mother_id", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -1;
  }, "cut on primary" );
  tru_pi_pos.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "sim_is_pionP");
  correction_task.AddVector(tru_pi_pos);

  VectorConfig tru_pi_neg( "tru_pi_neg", "sim_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_pi_neg.SetHarmonicArray( {1, 2} );
  tru_pi_neg.SetCorrections( {CORRECTION::PLAIN } );
  tru_pi_neg.SetCorrectionAxes( sim_pion_axes );
  tru_pi_neg.AddCut( "sim_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -211;
  }, "tru_pi_neg cut" );
  tru_pi_neg.AddCut( "sim_mother_id", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -1;
  }, "cut on primary" );
  tru_pi_neg.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "sim_is_pionM");
  correction_task.AddVector(tru_pi_neg);

  VectorConfig psi_rp( "psi_rp", "psi_rp", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  correction_task.Run();
}
