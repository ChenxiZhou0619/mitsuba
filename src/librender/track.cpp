#include <mitsuba/render/track.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

AnimatedTransform::AnimatedTransform(Stream *stream) {
	size_t nTracks = stream->readSize();
	if (nTracks == 0) {
		m_transform = Transform(stream);
	} else {
		for (size_t i=0; i<nTracks; ++i) {
			AbstractAnimationTrack::EType type = 
				(AbstractAnimationTrack::EType) stream->readUInt();
			AbstractAnimationTrack *track = NULL;
			switch (type) {
				case AbstractAnimationTrack::ETranslationX:
				case AbstractAnimationTrack::ETranslationY:
				case AbstractAnimationTrack::ETranslationZ:
				case AbstractAnimationTrack::EScaleX:
				case AbstractAnimationTrack::EScaleY:
				case AbstractAnimationTrack::EScaleZ:
				case AbstractAnimationTrack::ERotationX:
				case AbstractAnimationTrack::ERotationY:
				case AbstractAnimationTrack::ERotationZ:
					track = new FloatTrack(type, stream);
					break;
				case AbstractAnimationTrack::ETranslationXYZ:
				case AbstractAnimationTrack::EScaleXYZ:
					track = new VectorTrack(type, stream);
					break;
				case AbstractAnimationTrack::ERotationQuat:
					track = new QuatTrack(type, stream);
					break;
				default:
					Log(EError, "Encountered an unknown animation track type (%i)!", type);
			}

			track->incRef();
			m_tracks.push_back(track);
		}
	}
}

void AnimatedTransform::addTrack(AbstractAnimationTrack *track) {
	track->incRef();
	m_tracks.push_back(track);
}

AABB1 AnimatedTransform::getTimeBounds() const {
	if (m_tracks.size() == 0) 
#if !defined(__clang__)
		return AABB1(0.0f, 0.0f);
#else
	// HACK Workaround for clang
	{
		AABB1 b;
		b.min = b.max = 0.0f;
		return b;
	}
#endif

	Float min = std::numeric_limits<Float>::infinity();
	Float max = -std::numeric_limits<Float>::infinity();

	for (size_t i=0; i<m_tracks.size(); ++i) {
		const AbstractAnimationTrack *track = m_tracks[i];
		size_t size = track->getSize();
		SAssert(size > 0);
		min = std::min(min, track->getTime(0));
		max = std::max(max, track->getTime(size-1));
	}

#if !defined(__clang__)
	return AABB1(min, max);
#else
	// HACK Workaround for clang
	AABB1 b;
	b.min = min;
	b.max = max;
	return b;
#endif
}

AABB AnimatedTransform::getTranslationBounds() const {
	if (m_tracks.size() == 0) {
		Point p = m_transform(Point(0.0f));
		return AABB(p, p);
	}

	AABB aabb;

	for (size_t i=0; i<m_tracks.size(); ++i) {
		const AbstractAnimationTrack *absTrack = m_tracks[i];
		switch (absTrack->getType()) {
			case AbstractAnimationTrack::ETranslationX:
			case AbstractAnimationTrack::ETranslationY:
			case AbstractAnimationTrack::ETranslationZ: {
					int idx  = absTrack->getType() - AbstractAnimationTrack::ETranslationX;
					const FloatTrack *track = 
						static_cast<const FloatTrack *>(absTrack);
					for (size_t j=0; j<track->getSize(); ++j) {
						Float value = track->getValue(j);
						aabb.max[idx] = std::max(aabb.max[idx], value);
						aabb.min[idx] = std::min(aabb.min[idx], value);
					}
				}
				break;

			case AbstractAnimationTrack::ETranslationXYZ: {
					const VectorTrack *track = 
						static_cast<const VectorTrack *>(absTrack);
					for (size_t j=0; j<track->getSize(); ++j) 
						aabb.expandBy(Point(track->getValue(j)));
				}
				break;
			default:
				break;
		}
	}
	for (int i=0; i<3; ++i) {
		if (aabb.min[i] > aabb.max[i])
			aabb.min[i] = aabb.max[i] = 0.0f;
	}

	return aabb;
}

AABB AnimatedTransform::getSpatialBounds(const AABB &aabb) const {
	AABB result;

	if (m_tracks.size() == 0) {
		for (int j=0; j<8; ++j)
			result.expandBy(m_transform(aabb.getCorner(j)));
	} else {
		/* Compute approximate bounds */
		int nSteps = 100;
		AABB1 timeBounds = getTimeBounds();
		Float step = timeBounds.getExtents().x / (nSteps-1);

		for (int i=0; i<nSteps; ++i) {
			const Transform &trafo = eval(timeBounds.min.x + step * i);
			for (int j=0; j<8; ++j)
				result.expandBy(trafo(aabb.getCorner(j)));
		}
	}

	return result;
}

AnimatedTransform::~AnimatedTransform() {
	for (size_t i=0; i<m_tracks.size(); ++i)
		m_tracks[i]->decRef();
}

void AnimatedTransform::serialize(Stream *stream) const {
	stream->writeSize(m_tracks.size());
	if (m_tracks.size() == 0) {
		m_transform.serialize(stream);
	} else {
		for (size_t i=0; i<m_tracks.size(); ++i)
			m_tracks[i]->serialize(stream);
	}
}

void AnimatedTransform::TransformFunctor::operator()(const Float &t, Transform &trafo) const {
	Vector translation(0.0f);
	Vector scale(1.0f);
	Quaternion rotation;

	for (size_t i=0; i<m_tracks.size(); ++i) {
		AbstractAnimationTrack *track = m_tracks[i];
		switch (track->getType()) {
			case AbstractAnimationTrack::ETranslationX:
				translation.x = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ETranslationY: 
				translation.y = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ETranslationZ: 
				translation.z = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ETranslationXYZ: 
				translation = static_cast<VectorTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleX:
				scale.x = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleY: 
				scale.y = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleZ: 
				scale.z = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleXYZ: 
				scale = static_cast<VectorTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ERotationQuat: 
				rotation = static_cast<QuatTrack *>(track)->eval(t);
				break;
			default:
				Log(EError, "Encountered an unsupported "
					"animation track type: %i!", track->getType());
		}
	}

	trafo = Transform::translate(translation) * 
		rotation.toTransform() *
		Transform::scale(scale);
}

std::string AnimatedTransform::toString() const {
	if (m_tracks.size() == 0) {
		return m_transform.toString();
	} else {
		std::ostringstream oss;
		oss << "AnimatedTransform[tracks=" << m_tracks.size() << "]";
		return oss.str();
	}
}

MTS_IMPLEMENT_CLASS(AbstractAnimationTrack, true, Object)
MTS_IMPLEMENT_CLASS(AnimatedTransform, false, Object)
MTS_NAMESPACE_END
