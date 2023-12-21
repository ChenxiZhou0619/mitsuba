#include "shiftmap.h"
#include "../gradient_pt.h"

MTS_NAMESPACE_BEGIN

HalfVectorShiftResult halfVectorShift(Vector3 tangentSpaceMainWi,
                                      Vector3 tangentSpaceMainWo,
                                      Vector3 tangentSpaceShiftedWi,
                                      Float mainEta, Float shiftedEta) {
  HalfVectorShiftResult result;

  if (Frame::cosTheta(tangentSpaceMainWi) *
          Frame::cosTheta(tangentSpaceMainWo) <
      (Float)0) {
    // Refraction.

    // Refuse to shift if one of the Etas is exactly 1. This causes degenerate
    // half-vectors.
    if (mainEta == (Float)1 || shiftedEta == (Float)1) {
      // This could be trivially handled as a special case if ever needed.
      result.success = false;
      return result;
    }

    // Get the non-normalized half vector.
    Vector3 tangentSpaceHalfVectorNonNormalizedMain;
    if (Frame::cosTheta(tangentSpaceMainWi) < (Float)0) {
      tangentSpaceHalfVectorNonNormalizedMain =
          -(tangentSpaceMainWi * mainEta + tangentSpaceMainWo);
    } else {
      tangentSpaceHalfVectorNonNormalizedMain =
          -(tangentSpaceMainWi + tangentSpaceMainWo * mainEta);
    }

    // Get the normalized half vector.
    Vector3 tangentSpaceHalfVector =
        normalize(tangentSpaceHalfVectorNonNormalizedMain);

    // Refract to get the outgoing direction.
    Vector3 tangentSpaceShiftedWo =
        refract(tangentSpaceShiftedWi, tangentSpaceHalfVector, shiftedEta);

    // Refuse to shift between transmission and full internal reflection.
    // This shift would not be invertible: reflections always shift to other
    // reflections.
    if (tangentSpaceShiftedWo.isZero()) {
      result.success = false;
      return result;
    }

    // Calculate the Jacobian.
    Vector3 tangentSpaceHalfVectorNonNormalizedShifted;
    if (Frame::cosTheta(tangentSpaceShiftedWi) < (Float)0) {
      tangentSpaceHalfVectorNonNormalizedShifted =
          -(tangentSpaceShiftedWi * shiftedEta + tangentSpaceShiftedWo);
    } else {
      tangentSpaceHalfVectorNonNormalizedShifted =
          -(tangentSpaceShiftedWi + tangentSpaceShiftedWo * shiftedEta);
    }

    Float hLengthSquared =
        tangentSpaceHalfVectorNonNormalizedShifted.lengthSquared() /
        (D_EPSILON + tangentSpaceHalfVectorNonNormalizedMain.lengthSquared());
    Float WoDotH =
        abs(dot(tangentSpaceMainWo, tangentSpaceHalfVector)) /
        (D_EPSILON + abs(dot(tangentSpaceShiftedWo, tangentSpaceHalfVector)));

    // Output results.
    result.success  = true;
    result.wo       = tangentSpaceShiftedWo;
    result.jacobian = hLengthSquared * WoDotH;
  } else {
    // Reflection.
    Vector3 tangentSpaceHalfVector =
        normalize(tangentSpaceMainWi + tangentSpaceMainWo);
    Vector3 tangentSpaceShiftedWo =
        reflect(tangentSpaceShiftedWi, tangentSpaceHalfVector);

    Float WoDotH = dot(tangentSpaceShiftedWo, tangentSpaceHalfVector) /
                   dot(tangentSpaceMainWo, tangentSpaceHalfVector);
    Float jacobian = abs(WoDotH);

    result.success  = true;
    result.wo       = tangentSpaceShiftedWo;
    result.jacobian = jacobian;
  }

  return result;
}

ReconnectionShiftResult
reconnectShift(const Scene *scene, Point3 mainSourceVertex, Point3 targetVertex,
               Point3 shiftSourceVertex, Vector3 targetNormal, Float time) {
  ReconnectionShiftResult result;

  // Check visibility of the connection.
  if (!testVisibility(scene, shiftSourceVertex, targetVertex, time)) {
    // Since this is not a light sample, we cannot allow shifts through
    // occlusion.
    result.success = false;
    return result;
  }

  // Calculate the Jacobian.
  Vector3 mainEdge    = mainSourceVertex - targetVertex;
  Vector3 shiftedEdge = shiftSourceVertex - targetVertex;

  Float mainEdgeLengthSquared    = mainEdge.lengthSquared();
  Float shiftedEdgeLengthSquared = shiftedEdge.lengthSquared();

  Vector3 shiftedWo = -shiftedEdge / sqrt(shiftedEdgeLengthSquared);

  Float mainOpposingCosine =
      dot(mainEdge, targetNormal) / sqrt(mainEdgeLengthSquared);
  Float shiftedOpposingCosine = dot(shiftedWo, targetNormal);

  Float jacobian =
      std::abs(shiftedOpposingCosine * mainEdgeLengthSquared) /
      (D_EPSILON + std::abs(mainOpposingCosine * shiftedEdgeLengthSquared));

  // Return the results.
  result.success  = true;
  result.jacobian = jacobian;
  result.wo       = shiftedWo;
  return result;
}

ReconnectionShiftResult environmentShift(const Scene *scene, const Ray &mainRay,
                                         Point3 shiftSourceVertex) {
  const Emitter *env = scene->getEnvironmentEmitter();

  ReconnectionShiftResult result;

  // Check visibility of the environment.
  Ray offsetRay = mainRay;
  offsetRay.setOrigin(shiftSourceVertex);

  if (!testEnvironmentVisibility(scene, offsetRay)) {
    // Environment was occluded.
    result.success = false;
    return result;
  }

  // Return the results.
  result.success  = true;
  result.jacobian = Float(1);
  result.wo       = mainRay.d;

  return result;
}

MTS_NAMESPACE_END
