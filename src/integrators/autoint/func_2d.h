//* indirect incident light field

#include "siren.h"
#include <mitsuba/core/frame.h>
#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class Func_2d {
public:
  Func_2d() = default;

  void fit(const std::vector<std::pair<Vector, Spectrum>> &datas) { //

    // vector to theta, phi
    std::vector<Float> theta_phis;
    for (int i = 0; i < datas.size(); ++i) {
      Vector v          = datas[i].first;
      auto [theta, phi] = Frame::toThetaPhi(v);
      theta_phis.emplace_back(theta);
      theta_phis.emplace_back(phi);
    }
    // std vector to torch tensor
    auto          options = torch::TensorOptions().dtype(at::kFloat);
    torch::Tensor xs =
        torch::from_blob(theta_phis.data(), theta_phis.size(), options);
    xs.reshape({(long)datas.size(), 2});

    //! red
    {

      torch::optim::SGD optimizer(r->parameters(), 0.01f);

      std::vector<Float> rs;
      for (int i = 0; i < datas.size(); ++i) {
        rs.emplace_back(datas[i].second[0]);
      }
      torch::Tensor ys = torch::from_blob(rs.data(), rs.size(), options);
      ys.reshape({(long)rs.size(), 1});

      for (size_t epoch = 1; epoch <= m_epochs; ++epoch) {
        torch::Tensor prediction = r->forward(xs);

        torch::Tensor loss = torch::l1_loss(prediction, ys);
        optimizer.zero_grad();
        loss.backward();
        optimizer.step();
      }
    }

    //* green
    {

      torch::optim::SGD optimizer(g->parameters(), 0.01f);

      std::vector<Float> gs;
      for (int i = 0; i < datas.size(); ++i) {
        gs.emplace_back(datas[i].second[1]);
      }
      torch::Tensor ys = torch::from_blob(gs.data(), gs.size(), options);
      ys.reshape({(long)gs.size(), 1});

      for (size_t epoch = 1; epoch <= m_epochs; ++epoch) {
        torch::Tensor prediction = r->forward(xs);

        torch::Tensor loss = torch::l1_loss(prediction, ys);
        optimizer.zero_grad();
        loss.backward();
        optimizer.step();
      }
    }

    //? blue
    {

      torch::optim::SGD optimizer(b->parameters(), 0.01f);

      std::vector<Float> bs;
      for (int i = 0; i < datas.size(); ++i) {
        bs.emplace_back(datas[i].second[2]);
      }
      torch::Tensor ys = torch::from_blob(bs.data(), bs.size(), options);
      ys.reshape({(long)bs.size(), 1});

      for (size_t epoch = 1; epoch <= m_epochs; ++epoch) {
        torch::Tensor prediction = r->forward(xs);

        torch::Tensor loss = torch::l1_loss(prediction, ys);
        optimizer.zero_grad();
        loss.backward();
        optimizer.step();
      }
    }
  }

  //* Spectrum_1 : integration of siren
  //* Spectrum_2 : derivative of siren at a specific point
  std::pair<Spectrum, Spectrum> control_variate(Vector wo) const {
    //

    return {Spectrum(.0f), Spectrum(.0f)};
  }

private:
  SIREN r, g, b;

  static constexpr int m_epochs = 1000;
};

MTS_NAMESPACE_END