VERTEX:
#version 400
	
	void main() 
	{	
	}

GEOMETRY:
#version 400
	 
	layout(points) in;
	layout(triangle_strip, max_vertices = 4) out;

	void main() 
	{
		gl_Position = vec4( 1.0, 1.0, 0.5, 1.0 );
		EmitVertex();
		gl_Position = vec4(-1.0, 1.0, 0.5, 1.0 );
		EmitVertex();
		gl_Position = vec4( 1.0,-1.0, 0.5, 1.0 );
		EmitVertex();
		gl_Position = vec4(-1.0,-1.0, 0.5, 1.0 );
		EmitVertex();
		EndPrimitive(); 
	}

FRAGMENT:
#version 400

    uniform int Gamma		= 1;
    uniform int Contrast	= 1;
    uniform int Unsaturate	= 1;
    uniform int Tint		= 1;
    uniform int Vignetting	= 1;

	uniform vec2 viewport;
	uniform mat4 InverseRotationProjection;
	uniform vec4 Light0_position;
	uniform float AbsoluteTime;
	
	out vec4 fragColor;
	
	float hash1(vec2 x) {
		return fract(sin(dot(x, vec2(12.9898, 78.233)))*43758.5453);
	}
	
	float value1(in vec2 x) {		
		vec2 i = floor(x);
		float a = hash1(i+vec2(0,0));
		float b = hash1(i+vec2(1,0));
		float c = hash1(i+vec2(0,1));
		float d = hash1(i+vec2(1,1));
		vec2 f = smoothstep(0.0, 1.0, fract(x));
		return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
	}

	void main()
	{
		vec2 ndcoords = gl_FragCoord.xy/viewport*2.0 - 1.0;
		vec3 E = normalize((InverseRotationProjection * vec4(ndcoords,0,1)).xyz);
		vec3 L = normalize(Light0_position.xyz );
		
		// Atmospheric scattering
		vec3 sunColor	  = mix(vec3(0.80, 0.40, 0.20), vec3(1.00, 0.90, 0.75), smoothstep( 0.0, 0.3, L.z));
		vec3 zenithColor  = mix(vec3(0.01, 0.02, 0.04), vec3(0.35, 0.48, 0.60), smoothstep(-0.8, 0.0, L.z));
		vec3 horizonColor = mix(vec3(0.02, 0.03, 0.04), sunColor,               smoothstep(-0.4, 0.2, L.z));
		vec3 groundColor  = vec3( dot( mix(0.03*zenithColor, 1.4*zenithColor, smoothstep(0.0, 0.4, L.z)), vec3(0.22,0.33,0.45)) );
		
		vec3 color  = mix(zenithColor, horizonColor, pow(1.0-E.z, 4.0));
		
		// Sun and halo
		float sunHaloWidth = mix(2, 30, smoothstep(0.0, 0.4, L.z));
		float sunDiscWidth = 20000; 
		float EdotL = max(dot(E,L),0.0);
		color = mix(color, sunColor, pow(EdotL,sunHaloWidth));		
		color += 6.0*pow(EdotL, sunDiscWidth);										
		
		// Clouds
		color += 0.04*value1((E.xy + 0.01)/(E.z + 0.01)*2.0);

        // Horizon
        color = mix(color, groundColor, smoothstep(-0.1, 0.0, -E.z));
        
        // Postprocessing
	    float gamma = Gamma > 0.0 ? 2*2.2 : 3.7;
	    color = pow(color, vec3(1.0/gamma));
        color = mix(color, color*color*(3.0-2.0*color), 0.3 * Contrast);
	    color = mix(color, vec3(dot(color,vec3(0.33))), 0.4 * Unsaturate);
	    color *= mix(vec3(1.0), vec3(1.06, 1.05, 1.00), 1.0 * Tint);		
	    color *= mix(1.0, pow(2.0*(ndcoords.x*ndcoords.x-1.0)*(ndcoords.y*ndcoords.y-1), 0.20), 0.5 * Vignetting);
	
	    // Banding removal
		color.xyz += 0.02*(-0.5 + hash1(ndcoords + fract(AbsoluteTime)));

		fragColor = vec4(color.rgb, 1.0);
		gl_FragDepth = 1.0;
	}
