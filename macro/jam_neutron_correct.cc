#include <cmath>
#include <vector>

void jam_neutron_correct(std::string list, std::string calib_in_file="qa.root")
{

  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
    .Define("ndetSize","return ndetModPos.size();")
    .Define("ndetTotSig","float sum = 0.; for(auto& sig:ndetModSig) sum += sig; return sum;")
    .Define("ndetModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:ndetModPos) phi.push_back(pos.phi()); return phi;")
    .Define("ndetModX","ROOT::VecOps::RVec<float> x; for(auto& pos:ndetModPos) x.push_back(pos.x()); return x;")
    .Define("ndetModY","ROOT::VecOps::RVec<float> y; for(auto& pos:ndetModPos) y.push_back(pos.y()); return y;")
    .Filter("b < 20. && ndetSize > 0 && ndetTotSig > 0");

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("b|psiRP"));
  correction_task.SetChannelVariables({
    std::regex("ndetMod(Phi|Sig|Id|X|Y)"),
  });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"b", 20, 0., 10.} );

  VectorConfig psi_rp( "psi_rp", "psiRP", "Ones", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  VectorConfig ndet( "ndet", "ndetModPhi", "ndetModSig", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  ndet.SetHarmonicArray( {1, 2} );
  ndet.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  ndet.AddHisto2D({{"ndetModX", 1000, -1000, 1000}, {"ndetModY", 1000, -1000, 1000}});
  correction_task.AddVector(ndet);

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}
