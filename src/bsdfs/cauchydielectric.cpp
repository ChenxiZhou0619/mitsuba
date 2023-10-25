#include "ior.h"
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/render/bsdf.h>

MTS_NAMESPACE_BEGIN

// Dielectric transmission material with wavelength-dependent refractive index

// TODO configure specular reflectance and transmittance
class SmoothCauchy : public BSDF {
public:
  SmoothCauchy(const Properties &props) : BSDF(props) {
    // Dense flint glass SF10
    int_A = props.getFloat("int_A", 1.7280);
    int_B = props.getFloat("int_B", 13.42); // !int_B should be nanometer

    // Default vaccum
    ext_A = props.getFloat("ext_A", 1);
    ext_B = props.getFloat("ext_B", .0);
  }

  SmoothCauchy(Stream *stream, InstanceManager *manager) : BSDF(stream, manager) {
    // no
  }

  void serialize(Stream *stream, InstanceManager *manager) const {
    // no
  }

  void configure() {
    // no
  }

  void addChild(const std::string &name, ConfigurableObject *child) {
    // no
  }

  inline Vector reflect(const Vector &wi) const { return Vector{-wi.x, -wi.y, wi.z}; }

  inline Vector refract(const Vector &wi, Float cosThetaT, Float eta) const {
    Float scale = -(cosThetaT < 0 ? (1.f / eta) : eta);
    return Vector(scale * wi.x, scale * wi.y, cosThetaT);
  }

  Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const { return Spectrum(.0f); }

  Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const { return .0f; }

  // TODO pdf is wavelength-dependent
  //! In transmission case, only one channel should be valid
  Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
    bool sampleReflection =
        (bRec.typeMask & EDeltaReflection) && (bRec.component == -1 || bRec.component == 0);
    bool sampleTransmission =
        (bRec.typeMask & EDeltaTransmission) && (bRec.component == -1 || bRec.component == 1);

    Float cosThetaT;

    //* eta = int / ext
    Float eta = CauchyRefractiveIndex(bRec.lambda, int_A, int_B) /
                CauchyRefractiveIndex(bRec.lambda, ext_A, ext_B),
          invEta = 1.f / eta;

    Float F = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, eta);

    if (sampleTransmission && sampleReflection) {
      if (sample.x <= F) {
        bRec.sampledComponent = 0;
        bRec.sampledType      = EDeltaReflection;
        bRec.wo               = reflect(bRec.wi);
        bRec.eta              = 1.0f;
        pdf                   = F;

        return Spectrum(1.f);
      } else {
        bRec.sampledComponent = 1;
        bRec.sampledType      = EDeltaTransmission;
        bRec.wo               = refract(bRec.wi, cosThetaT, eta);
        bRec.eta              = cosThetaT < 0 ? eta : invEta;
        pdf                   = 1 - F;

        /* Radiance must be scaled to account for the solid angle compression
           that occurs when crossing the interface. */
        Float factor = (bRec.mode == ERadiance) ? (cosThetaT < 0 ? invEta : eta) : 1.0f;

        Spectrum value(.0f);
        value[bRec.channel] = factor * factor;
        return value;
      }
    } else if (sampleReflection) {
      bRec.sampledComponent = 0;
      bRec.sampledType      = EDeltaReflection;
      bRec.wo               = reflect(bRec.wi);
      bRec.eta              = 1.0f;
      pdf                   = 1.0f;

      return Spectrum(F);
    } else if (sampleTransmission) {
      bRec.sampledComponent = 1;
      bRec.sampledType      = EDeltaTransmission;
      bRec.wo               = refract(bRec.wi, cosThetaT, eta);
      bRec.eta              = cosThetaT < 0 ? eta : invEta;
      pdf                   = 1.0f;

      /* Radiance must be scaled to account for the solid angle compression
         that occurs when crossing the interface. */
      Float    factor = (bRec.mode == ERadiance) ? (cosThetaT < 0 ? invEta : eta) : 1.0f;
      Spectrum value(.0f);
      value[bRec.channel] = factor * factor * (1 - F);
      return value;
    }
    return Spectrum(0.0f);
  }

  Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
    bool sampleReflection =
        (bRec.typeMask & EDeltaReflection) && (bRec.component == -1 || bRec.component == 0);
    bool sampleTransmission =
        (bRec.typeMask & EDeltaTransmission) && (bRec.component == -1 || bRec.component == 1);

    Float cosThetaT;

    //* eta = int / ext
    Float eta = CauchyRefractiveIndex(bRec.lambda, int_A, int_B) /
                CauchyRefractiveIndex(bRec.lambda, ext_A, ext_B),
          invEta = 1.f / eta;

    Float F = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, eta);

    if (sampleTransmission && sampleReflection) {
      if (sample.x <= F) {
        bRec.sampledComponent = 0;
        bRec.sampledType      = EDeltaReflection;
        bRec.wo               = reflect(bRec.wi);
        bRec.eta              = 1.0f;

        return Spectrum(1.f);
      } else {
        bRec.sampledComponent = 1;
        bRec.sampledType      = EDeltaTransmission;
        bRec.wo               = refract(bRec.wi, cosThetaT, eta);
        bRec.eta              = cosThetaT < 0 ? eta : invEta;

        /* Radiance must be scaled to account for the solid angle compression
           that occurs when crossing the interface. */
        Float factor = (bRec.mode == ERadiance) ? (cosThetaT < 0 ? invEta : eta) : 1.0f;

        Spectrum value(.0f);
        value[bRec.channel] = factor * factor;
        return value;
      }
    } else if (sampleReflection) {
      bRec.sampledComponent = 0;
      bRec.sampledType      = EDeltaReflection;
      bRec.wo               = reflect(bRec.wi);
      bRec.eta              = 1.0f;

      return Spectrum(F);
    } else if (sampleTransmission) {
      bRec.sampledComponent = 1;
      bRec.sampledType      = EDeltaTransmission;
      bRec.wo               = refract(bRec.wi, cosThetaT, eta);
      bRec.eta              = cosThetaT < 0 ? eta : invEta;

      /* Radiance must be scaled to account for the solid angle compression
         that occurs when crossing the interface. */
      Float    factor = (bRec.mode == ERadiance) ? (cosThetaT < 0 ? invEta : eta) : 1.0f;
      Spectrum value(.0f);
      value[bRec.channel] = factor * factor * (1 - F);
      return value;
    }
    return Spectrum(0.0f);
  }

  Float getRoughness(const Intersection &its, int component) const { return .0f; }

  std::string toString() const {
    std::ostringstream oss;
    oss << "SmoothCauchy";
    return oss.str();
  }

  Shader *createShader(Renderer *renderer) const;

  MTS_DECLARE_CLASS()

protected:
  //* A + B / lambda^2
  Float CauchyRefractiveIndex(Float lambda, Float A, Float B) const {
    return A + B / (lambda * lambda);
  }

private:
  // Interior and exterior Cauchy coefficients
  Float int_A, int_B;
  Float ext_A, ext_B;
};
class SmoothCauchyShader : public Shader {
public:
  SmoothCauchyShader(Renderer *renderer) : Shader(renderer, EBSDFShader) { m_flags = ETransparent; }

  Float getAlpha() const { return 0.3f; }

  void generateCode(std::ostringstream &oss, const std::string &evalName,
                    const std::vector<std::string> &depNames) const {
    oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
        << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
        << "        return vec3(0.0);" << endl
        << "    return vec3(inv_pi * cosTheta(wo));" << endl
        << "}" << endl
        << endl
        << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
        << "    return " << evalName << "(uv, wi, wo);" << endl
        << "}" << endl;
  }

  MTS_DECLARE_CLASS()
};

Shader *SmoothCauchy::createShader(Renderer *renderer) const {
  return new SmoothCauchyShader(renderer);
}

MTS_IMPLEMENT_CLASS(SmoothCauchyShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothCauchy, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothCauchy, "Smooth cauchy dielectric BSDF");

MTS_NAMESPACE_END