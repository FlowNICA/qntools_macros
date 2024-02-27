//
// Created by Misha on 3/7/2023.
//

using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

void mpd_fixed_target_correct_w_eff(std::string list, std::string cm_energy="2.5", std::string str_nucleus_mass="209", std::string eff_file=""){

  const double sNN = std::stod( cm_energy ); // in GeV
  const double M = 0.938; // in GeV/c^2
  const double T = sNN*sNN/(2.*M) - 2.*M;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );

  std::cout << "sqrtSnn = " << sNN << " GeV; T = " << T << "A GeV; Y_BEAM = " << Y_BEAM << std::endl;
  std::cout << "A = " << nucleus_mass << "; R = " << NUCLEUS_RADIUS << std::endl;

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
  
  std::vector<int> f1_modules = { 15, 16, 17, 22, 24, 29, 30, 31 };
  std::vector<int> f2_modules = { 7, 8, 9, 10, 11, 14, 18, 21, 25, 28, 32, 35, 36, 37, 38, 39 };
  std::vector<int> f3_modules = { 1, 2, 3, 4, 5, 6, 12, 13, 19, 20, 26, 27, 33, 34, 40, 41, 42, 43, 44, 45 };
  
  TStopwatch timer;
  timer.Start();
  
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
    .Filter("mcB < 20.") // at least one filter is mandatory!!!
    .Define( "b_norm", [NUCLEUS_RADIUS](double b){ return (float)b / (float)NUCLEUS_RADIUS; },{"mcB"})
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
    .Define( "tr_pdg", "recoGlobalSimPdg" )
    .Define( "tr_primary", "recoGlobalSimMotherId==-1" )
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
    .Define( "tr_is_proton", "recoGlobalSimPdg==2212" )
    .Define( "tr_is_pionP", "recoGlobalSimPdg==211" )
    .Define( "tr_is_pionM", "recoGlobalSimPdg==-211" )
    .Define( "tr_is_pions", "abs(recoGlobalSimPdg)==211" )
    .Define( "tr_is_kaonP", "recoGlobalSimPdg==321" )
    .Define( "tr_is_kaonM", "recoGlobalSimPdg==-321" )
    .Define( "tr_is_kaons", "abs(recoGlobalSimPdg)==321" )
    .Define( "tr_nhits", "std::vector<int> vec_nhits; for( auto n : recoGlobalNhits ){ vec_nhits.push_back( n ); } return vec_nhits;" )
    .Define( "tr_dca", "vector <float> vec_dca; for ( auto dca : recoGlobalDca ){ vec_dca.push_back( sqrt( dca.Mag2() ) ); } return vec_dca;" )
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
  correction_task.SetEventVariables(std::regex("b_norm|psi_rp"));
  correction_task.SetChannelVariables({std::regex("fhcal_module_(id|phi|energy|x|y|z)")});
  correction_task.SetTrackVariables({
      std::regex("tr_(pT|y|eta|phi|charge|pdg|primary|is_proton|is_pionP|is_pionM|is_pions|nhits|dca|WeightProton|WeightPionP|WeightPionM|WeightPions|WeightKaonP|WeightKaonM|WeightKaons)"),
      std::regex("sim_(pT|y|phi|pdg|is_proton|is_pionP|is_pionM|is_pions|mother_id)")
    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( { "b_norm", 10, 0, 2.0 } );

  VectorConfig f1( "F1", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f1.AddCut( "fhcal_module_id", [f1_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcal_module_x", 100, -100, 100}, {"fhcal_module_y", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f2.AddCut( "fhcal_module_id", [f2_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcal_module_x", 100, -100, 100}, {"fhcal_module_y", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
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
  proton.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
  }, "proton cut" );
  proton.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  proton.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "q > 0" );
  proton.AddCut( "tr_nhits", [](double nhits){ return nhits > 27; }, "Nhits > 27" );
  proton.AddCut( "tr_dca", [](double dca){ return dca < 1.; }, "dca < 1 cm" );
  proton.AddCut( "tr_eta", [](double eta){ return eta < 2.; }, "no tracks in FHCal");
  proton.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightProton");
  correction_task.AddVector(proton);
  
  VectorConfig pi_pos( "pi_pos", "tr_phi", "tr_WeightPionP", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  pi_pos.SetHarmonicArray( {1, 2} );
  pi_pos.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  pi_pos.SetCorrectionAxes( pion_axes );
  pi_pos.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 211;
  }, "pi_pos cut" );
  pi_pos.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  pi_pos.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "q > 0" );
  pi_pos.AddCut( "tr_nhits", [](double nhits){ return nhits > 22; }, "Nhits > 22" );
  pi_pos.AddCut( "tr_dca", [](double dca){ return dca < 3.5; }, "dca < 3.5 cm" );
  pi_pos.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightPionP");
  correction_task.AddVector(pi_pos);
  
  VectorConfig pi_neg( "pi_neg", "tr_phi", "tr_WeightPionM", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  pi_neg.SetHarmonicArray( {1, 2} );
  pi_neg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  pi_neg.SetCorrectionAxes( pion_axes );
  pi_neg.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -211;
  }, "pi_neg cut" );
  pi_neg.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  pi_neg.AddCut( "tr_charge", [](double charge){ return charge < 0; }, "q < 0" );
  pi_neg.AddCut( "tr_nhits", [](double nhits){ return nhits > 22; }, "Nhits > 22" );
  pi_neg.AddCut( "tr_dca", [](double dca){ return dca < 3.5; }, "dca < 3.5 cm" );
  pi_neg.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightPionM");
  correction_task.AddVector(pi_neg);
  
  VectorConfig Tp( "Tp", "tr_phi", "tr_WeightProton", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1, 2} );
  Tp.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tp.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
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
  correction_task.AddVector(Tp);
  
  VectorConfig Tneg( "Tneg", "tr_phi", "tr_WeightPions", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return abs(pdg_code) == 211;
  }, "pion cut" );
  Tneg.AddCut( "tr_primary", [](double primary) { auto prim = std::round(primary); return prim == 1; }, "primary tracks");
  Tneg.AddCut( "tr_y", [](double eta){
    return -1.3 < eta && eta < -0.2;
  }, "Tneg y cut" );
  Tneg.AddCut( "tr_nhits", [](double nhits){ return nhits > 22; }, "Nhits > 22" );
  Tneg.AddCut( "tr_dca", [](double dca){ return dca < 3.5; }, "dca < 3.5 cm" );
  Tneg.AddHisto2D({{"tr_y", 300, -1.5, 1.5}, {"tr_pT", 100, 0.0, 2.0}}, "tr_WeightPions");
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
