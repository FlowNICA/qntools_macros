//
// Created by Misha on 4/6/2023.
//

#ifndef QNTOOLSINTERFACE_VECTOR_CONFIG_H
#define QNTOOLSINTERFACE_VECTOR_CONFIG_H

#include <Axis.hpp>
#include <CorrectionManager.hpp>

enum class VECTOR_TYPE{
  CHANNEL,
  TRACK
};

enum class NORMALIZATION{
  M,
  MAG,
  NONE,
};

enum class CORRECTION{
  PLAIN,
  RECENTERING,
  RESCALING
};

struct vector_cut;
struct histo1d;
struct histo2d;

struct vector_cut{
  std::string field;
  std::function<bool(double)> function;
  std::string description;
};

struct histo1d{
  Qn::AxisD axis;
  std::string weight{"Ones"};
};

struct histo2d{
  std::vector<Qn::AxisD> axes;
  std::string weight{"Ones"};
};


class VectorConfig {
public:
  VectorConfig(std::string name, std::string phi_field, std::string weight_field,
               VECTOR_TYPE type, NORMALIZATION normalization);
  VectorConfig(const VectorConfig&) = default;
  VectorConfig(VectorConfig&&) = default;
  VectorConfig& operator=(const VectorConfig&) = default;
  VectorConfig& operator=(VectorConfig&&) = default;
  ~VectorConfig() = default;

  void SetHarmonicArray(const std::vector<int> &harmonic_array) {
    harmonic_array_ = harmonic_array;
  }
  void SetCorrectionAxes(const std::vector<Qn::AxisD> &correction_axes) {
    correction_axes_ = correction_axes;
  }
  void SetCorrections(const std::vector<CORRECTION> &corrections) {
    corrections_ = corrections;
  }
  void AddCut( const std::string& field, const std::function<bool(double)>& function, const std::string& description );
  void AddHisto1D( const Qn::AxisD& axis, const std::string& weight= "Ones" );
  void AddHisto2D( const std::vector<Qn::AxisD>& axis, const std::string& weight= "Ones" );
  void Decorate( const std::shared_ptr<Qn::CorrectionManager>& man ) const;

private:
  std::string name_;
  std::string phi_field_;
  std::string weight_field_;
  VECTOR_TYPE type_;
  NORMALIZATION normalization_;
  std::vector<int> harmonic_array_{1};
  std::vector<CORRECTION> corrections_;
  std::vector<Qn::AxisD> correction_axes_{};
  std::vector<vector_cut> cuts_{};
  std::vector<histo1d> vec_histo1d_{};
  std::vector<histo2d> vec_histo2d_{};
};


#endif //QNTOOLSINTERFACE_VECTOR_CONFIG_H
