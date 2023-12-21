#pragma once

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/// Result of a half-vector duplication shift.
struct HalfVectorShiftResult {
  bool    success;  ///< Whether the shift succeeded.
  Float   jacobian; ///< Local Jacobian determinant of the shift.
  Vector3 wo;       ///< Tangent space outgoing vector for the shift.
};

/// Result of a reconnection shift.
struct ReconnectionShiftResult {
  bool    success;  ///< Whether the shift succeeded.
  Float   jacobian; ///< Local Jacobian determinant of the shift.
  Vector3 wo;       ///< World space outgoing vector for the shift.
};

/// Calculates the outgoing direction of a shift by duplicating the local
/// half-vector.
HalfVectorShiftResult halfVectorShift(Vector3 tangentSpaceMainWi,
                                      Vector3 tangentSpaceMainWo,
                                      Vector3 tangentSpaceShiftedWi,
                                      Float mainEta, Float shiftedEta);

/// Tries to connect the offset path to a specific vertex of the main path.
ReconnectionShiftResult
reconnectShift(const Scene *scene, Point3 mainSourceVertex, Point3 targetVertex,
               Point3 shiftSourceVertex, Vector3 targetNormal, Float time);

/// Tries to connect the offset path to a the environment emitter.

ReconnectionShiftResult environmentShift(const Scene *scene, const Ray &mainRay,
                                         Point3 shiftSourceVertex);

MTS_NAMESPACE_END