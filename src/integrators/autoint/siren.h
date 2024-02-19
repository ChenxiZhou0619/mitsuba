// clang-format off
#include <torch/torch.h>
#include <mitsuba/mitsuba.h>
// clang-format on

MTS_NAMESPACE_BEGIN

class SIRENImpl : public torch::nn::Module {
public:
  SIRENImpl() : fc1(2, 32), fc2(32, 64), fc3(64, 32), fc4(32, 1) {
    register_module("fc1", fc1);
    register_module("fc2", fc2);
    register_module("fc3", fc3);
    register_module("fc4", fc4);
  }

  torch::Tensor forward(torch::Tensor x) {
    x = torch::sin(fc1->forward(x));
    x = torch::sin(fc2->forward(x));
    x = torch::sin(fc3->forward(x));
    x = fc4->forward(x);
    return x;
  }

  torch::Tensor d2f_dxdy(torch::Tensor xy) {
    //

    auto e1 = torch::tensor({1.0, 0.0}).reshape({1, 2});
    auto e2 = torch::tensor({0.0, 1.0}).reshape({1, 2});

    auto d2v1_dxdy = e1.mm(torch::transpose(fc1->weight, 0, 1)) *
                     e2.mm(torch::transpose(fc1->weight, 0, 1)) *
                     (-torch::sin(fc1->forward(xy)));

    auto v1 = torch::sin(fc1->forward(xy));

    auto dv1_dx = e1.mm(torch::transpose(fc1->weight, 0, 1)) *
                  torch::cos(fc1->forward(xy));
    auto dv1_dy = e2.mm(torch::transpose(fc1->weight, 0, 1)) *
                  torch::cos(fc1->forward(xy));

    //
    auto d2v2_dxdy = d2v1_dxdy.mm(torch::transpose(fc2->weight, 0, 1)) *
                         torch::cos(fc2->forward(v1)) -
                     dv1_dx.mm(torch::transpose(fc2->weight, 0, 1)) *
                         dv1_dy.mm(torch::transpose(fc2->weight, 0, 1)) *
                         torch::sin(fc2->forward(v1));

    auto v2 = torch::sin(fc2->forward(v1));

    auto dv2_dx = dv1_dx.mm(torch::transpose(fc2->weight, 0, 1)) *
                  torch::cos(fc2->forward(v1));
    auto dv2_dy = dv1_dy.mm(torch::transpose(fc2->weight, 0, 1)) *
                  torch::cos(fc2->forward(v1));

    //
    auto d2v3_dxdy = d2v2_dxdy.mm(torch::transpose(fc3->weight, 0, 1)) *
                         torch::cos(fc3->forward(v2)) -
                     dv2_dx.mm(torch::transpose(fc3->weight, 0, 1)) *
                         dv2_dy.mm(torch::transpose(fc3->weight, 0, 1)) *
                         torch::sin(fc3->forward(v2));

    auto d2v4_dxdy = d2v3_dxdy.mm(torch::transpose(fc4->weight, 0, 1));

    return d2v4_dxdy;
  }

private:
  torch::nn::Linear fc1;
  torch::nn::Linear fc2;
  torch::nn::Linear fc3;
  torch::nn::Linear fc4;
};
TORCH_MODULE(SIREN);

MTS_NAMESPACE_END