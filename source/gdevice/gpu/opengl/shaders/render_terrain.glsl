VERTEX:
#version 400
// =========================================================
//  _____         _              _____ _         _         
// |  |  |___ ___| |_ ___ _ _   |   __| |_ ___ _| |___ ___ 
// |  |  | -_|  _|  _| -_|_'_|  |__   |   | .'| . | -_|  _|
//  \___/|___|_| |_| |___|_,_|  |_____|_|_|__,|___|___|_|  
//
// Given gl_VertexID, turns a height-quartet into 
// a LOD blended vertex position (in model space).

uniform vec2  tileOffset;
uniform float kernelSize;
uniform sampler2D quartetsTU;
int tileSize = textureSize(quartetsTU, 0).x;

out vec3 position;

void main() 
{	
	// Model-space vertex position in plane XY
	vec2 p = vec2( gl_VertexID % tileSize, gl_VertexID / tileSize );		
	
	// Get the quartet
	vec4 qu	= textureLod( quartetsTU, (p.xy + 0.5)/tileSize, 0 );

	// Compute blend factor (with LOD+1, in [kernelSize-2, kernelSize-1])
	float blend = 0.0;		
	vec2 offset = tileOffset + p/(tileSize-1);
	float chebyshev = max(abs(offset.x), abs(offset.y));
	blend = smoothstep(kernelSize-2, kernelSize-1, chebyshev);	
	
	// Blend position with coarser vertex position.
	p = mix(p, p-mod(p,2), blend);
	
	// Blend height with coarser vertex height.
	float height = mix(mix(qu.x, qu.y, blend), mix(qu.z, qu.w, blend), blend);
	
	position = vec3(p/(tileSize-1), height);
}

CONTROL:
#version 400  
// ===========================================================================================================                                                                                                        
//   _____                 _ _     _   _            _____         _           _    _____ _         _         
//  |_   _|___ ___ ___ ___| | |___| |_|_|___ ___   |     |___ ___| |_ ___ ___| |  |   __| |_ ___ _| |___ ___ 
//    | | | -_|_ -|_ -| -_| | | .'|  _| | . |   |  |   --| . |   |  _|  _| . | |  |__   |   | .'| . | -_|  _|
//    |_| |___|___|___|___|_|_|__,|_| |_|___|_|_|  |_____|___|_|_|_| |_| |___|_|  |_____|_|_|__,|___|___|_|  
//                                                                                                           
// * Patch frustum culling 
// * Tessellation control of edges and big patches (including vertical) 
// * Displace mapping

uniform vec2  tileOffset;
uniform float kernelSize;
uniform int   Tessellator = 1;
uniform float tessellation = 0.15;
uniform mat4  ModelViewProjectionMatrix;
uniform float displacement_range = 1.0;
uniform float scale;
uniform float lod_factor = 100.0;   // TEMP
uniform vec2  viewport;             // TEMP

in vec3 position[];

out vec3 tPosition[];
layout(vertices = 4) out;

vec4 projectToNDS(vec3 position) {
   vec4 result = ModelViewProjectionMatrix * vec4(position, 1.0);
   result /= result.w; // Scale to [-1,+1]
   return result;
}

bool isOffScreen(vec4 ndsPosition) {
    //return (ndsPosition.z < -0.5) || any(lessThan(ndsPosition.xy, vec2(-1.0)) || greaterThan(ndsPosition.xy, vec2(1.0)));
    return (ndsPosition.z < -1.0) || any(lessThan(ndsPosition.xy, vec2(-2.0)) || greaterThan(ndsPosition.xy, vec2(2.0)));
}
    
vec2 projectToScreen(vec4 ndsPosition) {
	return (clamp(ndsPosition.xy, -1.3, 1.3) + 1.0) * (viewport.xy*0.5);
}

float computeTessellationLevel(vec2 ss0, vec2 ss1) {
	return clamp(distance(ss0, ss1)/lod_factor, 0, 1);
}

void main()
{
	if( gl_InvocationID == 0 ) {
		vec4 p0 = projectToNDS(position[0]);
		vec4 p1 = projectToNDS(position[1]);
		vec4 p2 = projectToNDS(position[2]);
		vec4 p3 = projectToNDS(position[3]);
	
        if(all(bvec4(isOffScreen(p0), isOffScreen(p1), isOffScreen(p2), isOffScreen(p3)))) {
            // Discard patch (late frustum culling).
            gl_TessLevelOuter[0] = 0;
            gl_TessLevelOuter[1] = 0;
            gl_TessLevelOuter[2] = 0;
            gl_TessLevelOuter[3] = 0;
            gl_TessLevelInner[0] = 0;
            gl_TessLevelInner[1] = 0;
        } else {
            // TODO: Fetch (and lod blend) gradient, submap, mixmap
            // TODO: Tessellation control for plane change, tessellation of edges and big patch (also vertical).
            
        //// TEMP
			vec2 ss0 = projectToScreen(p0);
			vec2 ss1 = projectToScreen(p1);
			vec2 ss2 = projectToScreen(p2);
			vec2 ss3 = projectToScreen(p3);

			float e0 = computeTessellationLevel(ss1, ss2);
			float e1 = computeTessellationLevel(ss0, ss1);
			float e2 = computeTessellationLevel(ss3, ss0);
			float e3 = computeTessellationLevel(ss2, ss3);
        ////
			vec4 ox = tileOffset.x + vec4( position[0].x, position[1].x, position[2].x, position[3].x );
			vec4 oy = tileOffset.y + vec4( position[0].y, position[1].y, position[2].y, position[3].y );
			vec4 d = sqrt( ox*ox + oy*oy ); 
					//max(abs(ox), abs(oy));				
			d = 1.0 - smoothstep(0.0, displacement_range*(kernelSize-1), d);
				// clamp( (1-d/(kernelSize-1.0))/(1+d), 0, 1); 
			vec4 t = 1.0 + Tessellator * 63.0 * tessellation * float(scale == 1.0) * d ;//* (e0+e1+e2+e3)/4; // TEMP
			t = mix(t.yxwz, t.zyxw, 0.5); 
		
            gl_TessLevelOuter[0] = t.x;
            gl_TessLevelOuter[1] = t.y;
            gl_TessLevelOuter[2] = t.z;
            gl_TessLevelOuter[3] = t.w;
            gl_TessLevelInner[0] = mix(t.y, t.z, 0.5);
            gl_TessLevelInner[1] = mix(t.x, t.w, 0.5);
        }
    }

	tPosition[gl_InvocationID] = position[gl_InvocationID];
}

EVALUATION:
#version 400        
// ====================================================================================================================                                                                                                      
//   _____                 _ _     _   _            _____         _         _   _            _____ _         _         
//  |_   _|___ ___ ___ ___| | |___| |_|_|___ ___   |   __|_ _ ___| |_ _ ___| |_|_|___ ___   |   __| |_ ___ _| |___ ___ 
//    | | | -_|_ -|_ -| -_| | | .'|  _| | . |   |  |   __| | | .'| | | | .'|  _| | . |   |  |__   |   | .'| . | -_|  _|
//    |_| |___|___|___|___|_|_|__,|_| |_|___|_|_|  |_____|\_/|__,|_|___|__,|_| |_|___|_|_|  |_____|_|_|__,|___|___|_|  
//                                                                                                                     
// Computes position, gradient, color, mixmap (of the new vertex)
// interpolating within the generated patch (gl_TessCoord.xy).
// (position is still in model space).
// NOTE: The height is displaced but it takes lots of fetches.

// UNIFORMS:
uniform sampler2D gradientsTU;
uniform sampler2D colorsTU;
uniform sampler2D mixmapsTU;
uniform float scale;
int tileSize = textureSize(mixmapsTU, 0).x;

// TEMP: These are necessary for Fog and Displacement:
    uniform mat4 ModelViewMatrix; //
    
// TEMP: These are necessary for Fog:
    uniform float fog_density = 0.0005; //

// TEMP: These are necessary for displacement:
    uniform sampler2D detailsTU;
    uniform sampler2D detailsDxTU;
    uniform sampler2D detailsDyTU;
    uniform float kernelSize; //
    uniform float displacementRange = 1.0; //
    
uniform vec4 defaultColorR = vec4(0.36, 0.30, 0.26, 0.0); // Light stone
uniform vec4 defaultColorG = vec4(0.28, 0.24, 0.20, 0.0); // Pebbles
uniform vec4 defaultColorB = vec4(0.34, 0.26, 0.22, 0.0); // Brownish stone
uniform vec4 defaultColorA = vec4(0.50, 0.38, 0.30, 0.0); // Dirty sand

// INPUT:
layout(quads, fractional_odd_spacing, ccw) in;	
in vec3 tPosition[];

// OUTPUT:
out VertexData {
	vec4 position;
	vec4 gradient;
	vec4 color;
	vec4 mixmap;
} tVertex;
// gl_Position




////////////////////////////////////////////////////////////
// Common functions
//
vec3 getTriplanarWeightVector(vec3 N) {
    vec3 w = pow(abs(N), vec3(64.0)); 
    // w = max(w, 0.000000000000001); 
    return w / dot(w, vec3(1.0));
}

vec4 textureTriplanar(sampler2D detailsTU, vec3 position, vec3 weight, float lodBias)
{
    return weight.x*textureLod(detailsTU, position.zy, lodBias) 
         + weight.y*textureLod(detailsTU, position.zx, lodBias) 
         + weight.z*textureLod(detailsTU, position.xy, lodBias);
}

// https://www.shadertoy.com/view/Ms3yzS
vec4 blendMixmap(vec4 mixmap) 
{
    const float MixmapBlending = 0.5; // 0.0 sharp, 0.5 smooth, 1.0 smoother
    float maxima = max(max(mixmap.x, mixmap.y), max(mixmap.z, mixmap.w));
    mixmap = max(mixmap - maxima + MixmapBlending, vec4(0.0));
    mixmap /= dot(mixmap, vec4(1.0));
    return mixmap;
}

vec4 blendColor(vec4 color, vec4 mixmap, float x1, float x2)
{    
    color = mix(color, defaultColorR, smoothstep(x1, x2, mixmap.r))
          + mix(color, defaultColorG, smoothstep(x1, x2, mixmap.g))
          + mix(color, defaultColorB, smoothstep(x1, x2, mixmap.b))
          + mix(color, defaultColorA, smoothstep(x1, x2, mixmap.a));
    return color/4.0;
}

// http://weber.itn.liu.se/~stegu/TNM084-2017/bumpmapping/bumpmapping.pdf
// Also check doBumpMap() in Desert Canyon (https://www.shadertoy.com/view/Xs33Df)
void displaceVertexAndRecomputeNormal(inout vec3 p, inout vec3 N, float H, vec3 dH)
{
    p = p + H*N;
    N = normalize( (1.0 + dot(dH,N))*N - dH );
}
////////////////////////////////////////////////////////////


void main()
{
	float u = gl_TessCoord.x;
    float v = gl_TessCoord.y;
    
    vec3 position = mix(mix(tPosition[1],tPosition[0],u), mix(tPosition[2],tPosition[3],u), v);
	vec2 p = (position.xy*(tileSize-1) + 0.5)/tileSize;
	vec4 gradient = textureLod(gradientsTU, p, 0); 
	vec4 color = textureLod(colorsTU,  p, 0); 
	vec4 mixmap	= textureLod(mixmapsTU, p, 0);
	vec2 coords = position.xy * scale;	// TODO needed?

    // Used for Displacement and also for Fog later.
    // TODO Compute at vertex shader (it's E vector).
    float povDistance = length((ModelViewMatrix * vec4(position,1)).xyz);
		
#if 1   // Displace position along the normal in LOD zero.
	    float displacement = 1.0 - smoothstep(0.0, displacementRange*(kernelSize-1), povDistance);
	    if( displacement > 0.0 // Only within the 1st LOD (otherwise you get cracks).
           && color.a == 0.0 ) // Only for non water surfaces
	    { 
	    
	    	const float LodBias = 0.0;
	    	vec3 normal = normalize(vec3(gradient.xy/scale, 1));
            vec3 w = getTriplanarWeightVector(normal);
            vec4 luma4 = textureTriplanar(detailsTU, position.xyz, w, LodBias);    
            vec4 mixmap = blendMixmap(mixmap + luma4);
            float luma = 0.2;
            vec4 luma4Dx = textureTriplanar(detailsDxTU, position.xyz, w, LodBias);
		    vec4 luma4Dy = textureTriplanar(detailsDyTU, position.xyz, w, LodBias);

            float H = 0.01*dot(mixmap, luma4);
			vec3 dH = vec3(
			    dot(mixmap, luma4Dx),
			    dot(mixmap, luma4Dy),
			    0.0 // Dz doesn't matter because it's texture mapping, no relief here.
			);

			vec3 position1 = position.xyz;
			vec3 normal1   = normal.xyz;
            displaceVertexAndRecomputeNormal(position1, normal1, H, dH);
	    
	        color += 0.00000000001*blendColor(color, mixmap, 0.3, 06);
	    
/*	        vec4 details = pow(texture(detailsTU, coords.xy, +4.5), vec4(1.0));
		    
	        vec4 displacementWeights = vec4(0.120, 0.040, 0.030, 0.050);
	        float extrusion = dot(displacementWeights, mixmap);
	        float luma = extrusion * dot(details, mixmap);
	        luma = (luma - extrusion/2) * displacement;
            vec3 N = normalize(vec3(gradient.xy, 1));
		    position += luma * N;
            #if 0  // Re-adjusting gradient after displacement
                // FIX Discontinuity
		        gradient.xy += (1.0 + dot(luma*gradient.xy,N.xy))*N.xy - luma*gradient.xy;
            #endif 
            */
	    }	    
#endif  


	// Output vertex
	tVertex.position = vec4(coords, position.z,	povDistance ); // TODO: povDistance is not currently used down the pipeline
	tVertex.gradient = gradient;
    tVertex.color	 = color;
	tVertex.mixmap	 = mixmap;
    gl_Position = vec4(position, 1.0);
}


GEOMETRY:
#version 400     
// ===================================================================                                                          
//   _____                   _              _____ _         _         
//  |   __|___ ___ _____ ___| |_ ___ _ _   |   __| |_ ___ _| |___ ___ 
//  |  |  | -_| . |     | -_|  _|  _| | |  |__   |   | .'| . | -_|  _|
//  |_____|___|___|_|_|_|___|_| |_| |_  |  |_____|_|_|__,|___|___|_|  
//                                  |___|   
//                          
// Transforms from model space to world space.
// TODO: Insert detail.
// 
uniform mat3 NormalMatrix;
uniform mat4 ModelViewMatrix;
uniform mat4 ModelViewProjectionMatrix;
uniform vec4 Light0_position; //
uniform float scale; //

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;
 
in VertexData {
	vec4 position;
	vec4 gradient;
	vec4 color;
	vec4 mixmap;
} tVertex[];

out VertexData {
	vec4 position;  // Position in model space
	vec4 gradient;  // Gradients (T0 and T1)
	vec4 color;
	vec4 mixmap;
	vec3 N0;        // TEMP?
	vec3 E;		    // Position in world space
	vec3 L;         // Light in world space
	vec3 barycentric;
} gVertex;
	
void main() 
{
	for(int i=0; i<gl_in.length(); i++) 
	{
		gVertex.position    = tVertex[i].position;
		gVertex.gradient	= tVertex[i].gradient; 
		gVertex.color		= tVertex[i].color;
		gVertex.mixmap		= tVertex[i].mixmap;
		
		gVertex.N0	= NormalMatrix * normalize(vec3(tVertex[i].gradient.zw, 1));
		gVertex.E	= (ModelViewMatrix * gl_in[i].gl_Position).xyz; // length of it, is povDistance
		gVertex.L	= (ModelViewMatrix * vec4(Light0_position.xy/scale, Light0_position.zw)).xyz;	// TODO precompute?

		gVertex.barycentric =	i == 0 ? vec3(1,0,0) : 
								i == 1 ? vec3(0,1,0) : 
										 vec3(0,0,1) ;
										 
		gl_Position	= ModelViewProjectionMatrix * gl_in[i].gl_Position;
		EmitVertex();
	}
	EndPrimitive();
}


FRAGMENT:
#version 400  
// ===================================================================                                                                
//   _____                           _      _____ _         _         
//  |   __|___ ___ ___ _____ ___ ___| |_   |   __| |_ ___ _| |___ ___ 
//  |   __|  _| .'| . |     | -_|   |  _|  |__   |   | .'| . | -_|  _|
//  |__|  |_| |__,|_  |_|_|_|___|_|_|_|    |_____|_|_|__,|___|___|_|  
//                |___|     
// 
// Parallax mapping, bump mapping, texture mapping, AO, 
// diffuse, specular, indirect light, atmospheric scattering, 
// post effects.
	 
// UNIFORMS:
uniform int Wireframe   = 0;
uniform int DebugMode   = 0; const int ColorView = 1, HeightBlendView = 2, NormalView = 3, LightView = 4, FresnelView = 5; 
uniform int Bumps       = 1; //const int BumpsNormal = 1, BumpsMicro = 2;
uniform int Tessellator	= 1; // TEMP 
uniform int Diffuse		= 1;
uniform int Specular	= 1;
uniform int Indirect	= 1;
uniform int Sky			= 1;	
uniform int Fresnel		= 1;
uniform int Shadows		= 1;
uniform int Desaturate	= 1;
uniform int Scattering	= 1;
uniform int Gamma		= 1;
uniform int Contrast	= 1;
uniform int Unsaturate	= 1;
uniform int Tint		= 1;
uniform int Vignetting	= 1;

uniform sampler2D detailsTU;
uniform sampler2D detailsDxTU;
uniform sampler2D detailsDyTU;

uniform float bump_intensity = 0.400;

uniform vec4 defaultColorR = vec4(0.36, 0.30, 0.26, 0.0); // Light stone
uniform vec4 defaultColorG = vec4(0.28, 0.24, 0.20, 0.0); // Pebbles
uniform vec4 defaultColorB = vec4(0.34, 0.26, 0.22, 0.0); // Brownish stone
uniform vec4 defaultColorA = vec4(0.50, 0.38, 0.30, 0.0); // Dirty sand

//////////////////////////// ROCK GRIT BONE SAND
//uniform vec4  fresmap = vec4(1.0, 1.0, 0.2, 0.5);
//uniform vec4  frespow = vec4(5.0, 5.0, 2.0, 1.0);
uniform vec4  specmap = vec4(1.0, 1.0, 1.0, 1.0);
uniform vec4  specpow = vec4(20,  20,  20,  20 );

uniform vec2 viewport;
uniform mat4 InverseRotationProjection;
uniform vec4 Light0_position;
uniform float visibileDistance = 100;
//uniform float AbsoluteTime;

uniform mat3 NormalMatrix; //// TEMP ////
uniform float scale; // TEMP
	
// INPUT:
in VertexData {
	vec4 position; // Position in model space
    vec4 gradient; // Gradients (T0 and T1)
	vec4 color;
	vec4 mixmap;
	vec3 N0;       // TEMP
	vec3 E;        // Position in world space
	vec3 L;        // Light in world space
	vec3 barycentric;
} gVertex;
 
// OUTPUT:
out vec4 fragColor;


////////////////////////////////////////////////////////////
// Common functions
//
vec3 getTriplanarWeightVector(vec3 N) {
    vec3 w = pow(abs(N), vec3(64.0)); 
    // w = max(w, 0.000000000000001); 
    return w / dot(w, vec3(1.0));
}

vec4 textureTriplanar(sampler2D detailsTU, vec3 position, vec3 weight, float lodBias)
{
    return weight.x*texture(detailsTU, position.zy, lodBias) 
         + weight.y*texture(detailsTU, position.zx, lodBias) 
         + weight.z*texture(detailsTU, position.xy, lodBias);
}

// https://www.shadertoy.com/view/Ms3yzS
vec4 blendMixmap(vec4 mixmap) 
{
    const float MixmapBlending = 0.5; // 0.0 sharp, 0.5 smooth, 1.0 smoother
    float maxima = max(max(mixmap.x, mixmap.y), max(mixmap.z, mixmap.w));
    mixmap = max(mixmap - maxima + MixmapBlending, vec4(0.0));
    mixmap /= dot(mixmap, vec4(1.0));
    return mixmap;
}

vec4 blendColor(vec4 color, vec4 mixmap, float x1, float x2)
{
    color = mix(color, defaultColorR, smoothstep(x1, x2, mixmap.r))
          + mix(color, defaultColorG, smoothstep(x1, x2, mixmap.g))
          + mix(color, defaultColorB, smoothstep(x1, x2, mixmap.b))
          + mix(color, defaultColorA, smoothstep(x1, x2, mixmap.a));
    return color/4.0;
}

// http://weber.itn.liu.se/~stegu/TNM084-2017/bumpmapping/bumpmapping.pdf
// Also check doBumpMap() in Desert Canyon (https://www.shadertoy.com/view/Xs33Df)
void displaceVertexAndRecomputeNormal(inout vec3 p, inout vec3 N, float H, vec3 dH)
{
    p = p + H*N;
    N = normalize( (1.0 + dot(dH,N))*N - dH );
}
////////////////////////////////////////////////////////////
	
// TODO The whole environment sphere texture should be precomputed when time (and location?) changes significantly.
vec3 getSkyDomeColor(vec3 E, vec3 L, vec3 sunColor, vec3 zenithColor, vec3 horizonColor, vec3 groundColor)
{
    vec3 color  = mix(zenithColor, horizonColor, pow(1.0-E.z, 4.0));
		
	// sun and halo
	float sunHaloWidth = mix(2, 30, smoothstep(0.0, 0.5, L.z));
	float sunDiscWidth = 5000; 
	float EdotL = max(dot(E,L),0.0);
	color = mix(color, sunColor, pow(EdotL,sunHaloWidth));		
	color += pow(EdotL, sunDiscWidth); //											
	
	// clouds
	//color += 0.03*value1((E.xy+0.1)/(E.z+0.1)*2.0);
	
	// horizon
    color = mix(color, groundColor, smoothstep(-0.1, 0.0, -E.z));

	// NOTE no postprocessing
	return color;
}

vec3 tone(vec3 color, float t) 
{
    return smoothstep(0.0, 1.0, color/(color + t)*(1.0 + t));
    //return color/(color + t)*(1.0 + t);
}
	

void main()
{	
	//
	// Texture generation
	//
	// NOTE: Here height, luma, occlusion are assimilated to the same physical entity.
	float vertexLuma = 0.2;
	float luma      = vertexLuma;
	vec3 normal     = normalize(vec3(gVertex.gradient.xy, 1));
	vec4 color      = gVertex.color; 
	vec4 mixmap     = gVertex.mixmap;

    float dist      = length(gVertex.E);
	float textureFadingFactor = clamp(dist/visibileDistance, 0.0, 1.0);
	
	if( textureFadingFactor <= 1.0 )
	{
        const float LodBias = -0.666; // 0.0 smooth, -0.5 sharper
        vec3 w = getTriplanarWeightVector( normalize(vec3(gVertex.gradient.xy/scale, 1)) );
        vec4 luma4 = textureTriplanar(detailsTU, gVertex.position.xyz, w, LodBias); 
        // TODO: de-offst for mixmap>1, meaning: offset [1..2]=>[1..0] and amplitude 1=>0.
            
        mixmap = blendMixmap(mixmap + luma4); // moveable to tessellator ?
        color  = blendColor(gVertex.color, mixmap, 0.3, 0.6); // moveable to tessellator ?
		
		// Graciously fade texture mapping out with distance.
		luma = mix( dot(mixmap, luma4), vertexLuma, textureFadingFactor);
		color = mix(color, gVertex.color, textureFadingFactor);

		if( Bumps > 0.0 ) 
		{    
		   vec4 luma4Dx = textureTriplanar(detailsDxTU, gVertex.position.xyz, w, LodBias);
		   vec4 luma4Dy = textureTriplanar(detailsDyTU, gVertex.position.xyz, w, LodBias);

            // Compute H and dH
            // NOTE: Height displacement is disabled here
            float H = 0.0; // Height alteration doesn't matter in bump mapping.
			vec3 dH = vec3(
			    dot(mixmap, luma4Dx),
			    dot(mixmap, luma4Dy),
			    0.0 // Dz doesn't matter because it's texture mapping, no relief here.
			);
			
			// Energy adujustment across LODs (empiric).
			float scale00 = 1 + 0.70*log(scale); 
			
			// Graciously fade bumps out with distance.
		    float _bump_intensity = mix(bump_intensity, 0.0, textureFadingFactor);
			
			//H  *= _bump_intensity * scale00; 
			dH *= _bump_intensity * scale00;
			
			vec3 P = gVertex.position.xyz;
            displaceVertexAndRecomputeNormal(P, normal, H, dH);
		}
	}
	
	//
	// Global Illumination
	//
	
	// Model space
    vec2 ndcoords = gl_FragCoord.xy/viewport*2.0 - 1.0;
	vec3 EE = normalize((InverseRotationProjection * vec4(ndcoords,0,1)).xyz);
	vec3 LL = normalize(Light0_position.xyz);
	vec3 LI = normalize(vec3(-LL.x, -LL.y, 0.0));
	vec3 NN = normalize(vec3(normal.xy/scale, normal.z));
	vec3 RR = reflect(LL,NN);
	
	// atmospheric scattering
	vec3 sunColor	  = mix(vec3(0.80, 0.24, 0.20), vec3(1.00, 0.90, 0.74), smoothstep( 0.0, 0.4, LL.z));
	vec3 zenithColor  = mix(vec3(0.01, 0.02, 0.04), vec3(0.35, 0.48, 0.60), smoothstep(-0.8, 0.0, LL.z));
	vec3 horizonColor = mix(vec3(0.02, 0.03, 0.04), sunColor,               smoothstep(-0.4, 0.5, LL.z));
	vec3 groundColor  = vec3( dot( mix(0.03*zenithColor, 1.4*zenithColor, smoothstep(0.0, 0.4, LL.z)), vec3(0.33)) );
	vec3 skyColorByEE = getSkyDomeColor(EE, LL, sunColor, zenithColor, horizonColor, groundColor);
	
	//
	// Lighting
	//
	vec3 col;
	vec3 light;
    
    vec3 E = normalize(gVertex.E); 
	vec3 L = normalize(gVertex.L);
    vec3 N = normalize(NormalMatrix * normal);
	vec3 R = reflect(L,N);
	
	float shadowing = clamp(20.0*dot(normalize(gVertex.N0),L), 1.0 - Shadows, 1.0);
	float indirectL  = smoothstep(-0.5, 0.5, LL.z);
    float occlusion = pow(luma, 1.0);
	float glossmap  = smoothstep(0.2, 0.7, luma); // TODO add sub glossmap for sand
	float fresnmap  = smoothstep(0.2, 0.7, luma);
	
	float diffuse   = shadowing * occlusion * max(0.0, dot(N,L));
	float specular  = shadowing * glossmap  * dot(mixmap, specmap) * pow(max(0.0, dot(E,R)), dot(mixmap, specpow));
	float indirect  = indirectL * occlusion * max(0.0, dot(NN,LI)); 
	float zenith    = 1.0       * occlusion * clamp(0.5 + 0.5*N.z, 0.0, 1.0);
	float fresnel   = 1.0       * fresnmap  * mix( 1.0, pow(1.0-abs(dot(E,N)), 4.0), 0.9 );

    // Energy conservation (?)
    //float reflectivity = 0.40; // reflectivity
    light  = 0.60 * Diffuse  * diffuse  * sunColor;
    light += 0.40 * Specular * specular * sunColor; //getSkyDomeColor(RR, LL, sunColor, zenithColor, horizonColor, groundColor);
    
    // Ambient
    light += 0.04 * Indirect * indirect * sunColor;
    light += 0.12 * Sky	     * zenith   * zenithColor;
    light = mix(light, getSkyDomeColor(reflect(EE,NN), LL, sunColor, zenithColor, horizonColor, groundColor), 0.10 * Fresnel * fresnel);
    
    // Tone mapping
    vec3 desaturate = 0.70 * Desaturate * smoothstep(0.1, 0.4, LL.z) * tone(vec3(diffuse + specular), 0.02);
	col = light * mix(color.rgb, sunColor, desaturate);
    col = tone(col, 0.2); 

    //
	// Scattering (halo and fog)
	//
	if( Scattering > 0 ) {
	    float EdotL = max(0.0, dot(E,L));
	    float dayTime = smoothstep(-0.05, 0.5, LL.z);
	    float scattering = dayTime*(1.0-exp(-0.0030*dist))*pow(EdotL, 8.0); // TODO: and visibileDistance?
		col += 8.0 * scattering * zenithColor;
		
		float volumetric = 100.0 * Bumps * pow(0.01*gVertex.position.z, 6.0); // TODO
		col = mix(col, groundColor, smoothstep(visibileDistance*0.2, visibileDistance, dist + 0.0*volumetric));
	}
	
	// FIXME: unsaturate a bit at dawn because of too much red light.
	col = mix(col, vec3(dot(col,vec3(0.33))), 0.5 * (1.0 - smoothstep(-0.1, 0.5, LL.z)) );

	//
	// Postprocessing
	//
	
	float gamma = Gamma > 0.0 ? 2*2.2 : 3.7;
	col = pow(col, vec3(1.0/gamma));	
    col = mix(col, col*col*(3.0-2.0*col), 0.3 * Contrast);
	col = mix(col, vec3(dot(col,vec3(0.33))), 0.4 * Unsaturate);
	col *= mix(vec3(1.0), vec3(1.06, 1.05, 1.00), 1.0 * Tint);	
	
	col *= mix(1.0, pow(2.0*(ndcoords.x*ndcoords.x-1.0)*(ndcoords.y*ndcoords.y-1), 0.20), 0.5 * Vignetting);	
	col += 0.00000000001 * Gamma * Contrast * Unsaturate * Tint * Vignetting*Desaturate; // TEMP

	// banding removal
	//col.xyz += hash1b( ndcoords + fract(AbsoluteTime*0.0001) )/255.0;
	//col.xyz += 10*hash12( ndcoords*123 + fract(AbsoluteTime*0.01) )/255.0;
	
	
	//
	// Debug controls
	//
   	if( DebugMode == ColorView ) {
	    col = 4.0 * color.rgb * occlusion;
	} else if( DebugMode == LightView ) {
		col = 3.0 * light; // NOTE: Fresnel is not included.
	} else if( DebugMode == NormalView ) {
		col = normalize(vec3(normal.xy/scale, normal.z));
	} else if( DebugMode == FresnelView ) {
		col.xyz = mix(col.xyz, vec3(0.0,1.0,1.0), fresnel);
	} 
	
	if( Wireframe > 0 ) {
	    // Wireframe
		vec3 d = fwidth(gVertex.barycentric);
		vec3 t = smoothstep(vec3(0.0), d*1.5, gVertex.barycentric);
		float edgeFactor = min(min(t.x, t.y), t.z);
		col = mix(pow(col.xyz,vec3(0.8)), pow(col.xyz,vec3(1.2)), edgeFactor);	    
	    // TileFrame   
	    vec2 p = fract(gVertex.position.xy/scale);
	    float m = min(min(p.x, p.y), min(1.0-p.x, 1.0-p.y));
	    col.r += 0.1*(1.0-smoothstep(0.0, 0.05, m)); 
	}
	
	if( gl_FragCoord.x < 700.0 ) {
	    //col = vec3(1,0,0);
	}

	fragColor = vec4( col.rgb, 1 );
}
