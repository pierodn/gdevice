#pragma once

#include "utils/glsl.h"

// TODO height occlusion culling (?)
//		http://www.lighthouse3d.com/opengl/viewfrustum/index.php?intro
//		http://www.lighthouse3d.com/opengl/viewfrustum/index.php?gasource
//		http://www.gamedev.net/community/forums/topic.asp?topic_id=487066
//
// NOTE Since this is the camera-centered case, it is possible to save near and far planes tests :)


struct Frustum 
{
	vec4 planes[6];

	void setup(const mat4 &m)
	{
/*
		const vec4 r0(m.elem[0][0], m.elem[0][1], m.elem[0][2], m.elem[0][3]);
		const vec4 r1(m.elem[1][0], m.elem[1][1], m.elem[1][2], m.elem[1][3]);
		const vec4 r2(m.elem[2][0], m.elem[2][1], m.elem[2][2], m.elem[2][3]);
		const vec4 r3(m.elem[3][0], m.elem[3][1], m.elem[3][2], m.elem[3][3]);
*/
		planes[0] = m.r3 + m.r0;
		planes[1] = m.r3 - m.r0;
		planes[2] = m.r3 + m.r1;
		planes[3] = m.r3 - m.r1;
		planes[4] = m.r3 + m.r2;
		planes[5] = m.r3 - m.r2;
	}

	int aabb_vs_frustum(const AABB &aabb, const Frustum &f)
	{
		int result = 1;

		for (int i = 0; i < 6; ++i) 
		{
			const vec4 &plane = f.planes[i];

			const vec3 pv(
				plane.x > 0 ? aabb.max.x : aabb.min.x,
				plane.y > 0 ? aabb.max.y : aabb.min.y,
				plane.z > 0 ? aabb.max.z : aabb.min.z
			);

			const vec3 nv(
				plane.x < 0 ? aabb.max.x : aabb.min.x,
				plane.y < 0 ? aabb.max.y : aabb.min.y,
				plane.z < 0 ? aabb.max.z : aabb.min.z
			);

			const float n = dot( vec4(pv, 1.0f), plane );
			if (n < 0) return -1;

			const float m = dot( vec4(nv, 1.0f), plane );
			if (m < 0) result = 0;
		}

		return result;
	}

};