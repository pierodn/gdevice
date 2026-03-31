compute:
#version 430 core

layout(local_size_x = 16, local_size_y = 16) in;
layout(binding=0) uniform sampler2D detailsTU;
layout(binding=1, rgba32f) uniform writeonly image2D detailsDxIU;
layout(binding=2, rgba32f) uniform writeonly image2D detailsDyIU;

vec4 getDetail(vec2 uv, float lod)
{
    vec4 det4 = texture(detailsTU, uv, lod).xyzw;
    //det4 = pow(det4,vec4(0.45));
    return det4;
}

void getGradient(vec2 uv, float size, out vec4 dx, out vec4 dy)
{
    dx = dy = vec4(0.0);
    float lods = log2(size);
    for(float lod = 0.0; lod<lods; lod++) {
        float delta = pow(2.0, lod)/size;
        vec2 offset = vec2(delta,0);
        dx += (getDetail(uv - offset.xy, lod) - getDetail(uv + offset.xy, lod))/2.0;
        dy += (getDetail(uv - offset.yx, lod) - getDetail(uv + offset.yx, lod))/2.0;
    }
}

void main() 
{
	ivec2 ij = ivec2(gl_GlobalInvocationID.xy);
	
	ivec2 size = textureSize(detailsTU, 0);
    if(ij.x >= size.x || ij.y >= size.y) {
        return;
    }
     
    vec4 dx, dy;
    getGradient(vec2(ij)/size.x, size.x, dx, dy);
    
    imageStore( detailsDxIU, ij, (dx + 1.0)/2.0 );
    imageStore( detailsDyIU, ij, (dy + 1.0)/2.0 );
}