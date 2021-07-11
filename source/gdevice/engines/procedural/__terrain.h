#pragma once

// TODO move gradients as external params, at world level
// TODO saturate to snow to emulate poles
// TODO saturate to desert to emulate equator line


struct Terrain
{
	inline vec3 compute_vertex( int i, int j, Uniforms& uniforms, VertexBuffer* vbo = NULL )
	{
		//
		// World space point
		//
		float u = (float)i/(uniforms.in_width-1);
		float v = (float)j/(uniforms.in_height-1);
		vec2 point = uniforms.in_offset + uniforms.in_size*vec2(u,v);
		point /= scale.xy;

		//
		// Looping a world
		//
		point = amod(point,circumference);
		vec2 p = 2.0f*point/circumference - 1.0f;
		vec3 clipper = vec3(
			-2*PI/360*p.x*sin(p.x*p.x*PI/2)*cos(p.y*p.y*PI/2),
			-2*PI/360*p.y*sin(p.y*p.y*PI/2)*cos(p.x*p.x*PI/2),
			cos(p.x*p.x*PI/2)*cos(p.y*p.y*PI/2));

		//
		// Shaping continents
		// 
		vec3 plate1 = warp<value>( point, mat2(+0.007, -0.027, +0.031, +0.006), 0.3 );
		vec3 plate2 = warp<value>( point, rotate(14.0f)*mat2(0.0013,-0.019,0.017,0.0011) );
		vec3 lands  = 1.0f * mul(plate1,plate2);

		//
		// Heigh frequency details and erosion features
		//
		//vec3 mounts = 7.0f*fbm2_eroded_s2<value>( point, 16, mat2(0.05), mat2(1.6,-1.2,1.2,1.6), 2.1, 0.38, 0.0001 );
/*
		float lod = log2(uniforms.in_size);
		bool maybeOnHigherLod = (i==0 || j==0 || i==uniforms.in_width-1 || j==uniforms.in_height-1) && (i%2==1 || j%2==1);
		float octaves = maybeOnHigherLod ? 10 + (20-10) * (1 - smoothstep(0.0f,9.0f,lod-1)) : 
									       10 + (20-10) * (1 - smoothstep(0.0f,9.0f,lod));
*/
		vec3 mounts = 10.0f*fbm2_eroded3<value>( point, 20, mat2(0.05), mat2(1.6,-1.2,1.2,1.6), 2.1, 0.38, 0.0001 );
		

		//
		// Composition
		//
		vec3 mixture = 1.9f*plate1 - 0.5f*mul(mounts,2.0f*plate2);
		vec3 n = lerp(crop(mixture), bias(-2.0,lands), lands+mounts) + mounts;

		// river level
		vec3 river = abs( bias(river_level,n) );
		// NOTE given for abs(n+river_level)==1
		vec3 nr = vec3( n.xy, -river_level + 1 );
		nr = lerp( plate2-plate1, nr, vec3(n.xy,1) );
		nr = lerp( crop(lands,0.4f*sea,1.0f), lands, lands + nr );
		nr = mul(nr,clipper);
		// joining river/sea
		//float ss = clamp( (nr.z-sea)/sea );
		//nr = mix( vec3(nr.xy,sea), nr, 1-ss );
		///////////////

		n = lerp( plate2-plate1, n, river );
		n = lerp( crop(lands,0.4f*sea,1.0f), lands, lands + n );

		//
		// Clipping
		//
		n = mul(n,clipper);

		float altitude = n.altitude; // + mix_factor; //(mix_factor+random)*flatness;

		vec3 normal = normalize( vec3(-scale.z * n.xy / scale.xy, 1.0f) );
		normal.xy *= uniforms.in_size;

		// surface properties	
		float slope = 1-normal.z;
		float pattern = rand( (n.y+1)/(n.x+1) * (float)0x10000000 );

		//vec4 rock_or_cobs = vec4( mix(ROCK.xyz,SAND.xyz,clamp(pattern-0.0f)), 
		//			0.5f*(1+smoothstep(0.05f,0.10f,slope)) );

		vec4 rock_or_cobs = mix(ROCK, COBS, 0.5f*(1+smoothstep(0.05f,0.10f,slope)) );

		// TODO forget about vec8
		vec8 terrain = slope<0.20 ? 
			vec8( BROWN, rock_or_cobs) :		// <== COBS!!
			vec8( STONE, ROCK);

		// Grass
		if( slope-pattern*0.02f<0.08 && slope*altitude<sea*1.0140f ) 
		{
			vec8 MLAWN = mix( vec8(GREEN4_TEMP, 
							mix( vec4(ROCK.rgb,0), // 1
								LAWN, 
								clamp(pattern+0.4f)) ), 
					vec8(GREEN4,LAWN), clamp(altitude*pattern-0.0f) );
			terrain = mix( vec8(GREEN3,LAWN), MLAWN, smoothstep(0.00f,0.02f,slope*altitude) );
		}


		//
		// Beach
		//
		if( altitude<=sea*1.03f + slope*0.001f )
		{
			terrain = mix( vec8(BEACH,SAND), terrain, smoothstep(sea*1.00f, sea*1.03f, (1-slope*0.5f)*altitude) ); // +sand
		}

		// 
		// River
		//
		if( river.z<=1 )
		{
			float k = 1 - (river.z - 0.9)*10;
			//terrain.color  = mix( mix(2.6f*RIVER,terrain.color,0.3f), RIVER, clamp(6*k+0.2) ); // foam
			terrain.color  = mix( mix(2.6f*RIVER,terrain.color,0.9f), RIVER, clamp(6*k+0.2) );
			terrain.detail = mix( terrain.detail, SAND, clamp(8*k+0.2) );
			n.altitude = nr.z;
		}

		//
		// Sea
		//
		if( altitude<=sea )
		{
			float t = altitude/sea;
			float k = 1 - t*t*t*t*t*t*t*t*t*t;
			terrain.color  = mix( terrain.color, WATER,	clamp(1*k+0.2) );
			terrain.detail = mix( terrain.detail, SAND, clamp(3*k+0.2) );
			n.altitude = sea;
		} 
		else // Snow
		{
			float melting = 70.0f;
			//terrain = mix( terrain, vec8(WHITE,SAND), smoothstep(snow+0.1f*n.y, (snow+0.4f)-0.1f*n.x, (1-slope*slope*melting)*altitude ));
			//terrain = mix( terrain, vec8(SLATE,COBS), smoothstep( snow-(0.01f-0.05f*n.y), snow-(0.01f+0.03f*n.x), clamp((slope-0.5f)*altitude) ));
			terrain = mix( terrain, vec8(WHITE,SAND), smoothstep( snow+0.0f*n.y, (snow+0.4f)-0.0f*n.x, (1-slope*slope*melting)*altitude ));
		}		
/*
		//if( river.z>=1 ) // Haze
		{
			//terrain.color = clamp( terrain.color * 0.4f*log2(uniforms.in_size)/(river.z*river.z), vec4(0), vec4(1) );
			//terrain.color = clamp( terrain.color * 1.3f/(river.z*river.z), vec4(0), vec4(1) );
			float haze = 0.8f;
			float k = river.z-0.0f;
			k = haze/(k*k);
			//k = clamp(k,0.8f,haze);
			//terrain.color = clamp( terrain.color * k, vec4(0), vec4(1) );
			terrain.color = mix( terrain.color, WHITE, clamp(k) );
			
		}
*/

		// attempt to insert artifacts such as buildings
		if( add_artifacts && length(mounts.xy)<0.15 && altitude>sea && river.z>1)
		{
			if( (i+j+int(mod(mounts.z*1000.0f,2.0f))) ^1 )
			{
				n.altitude += 0.01f + mod(mounts.z*1000.0f,0.02f);
				terrain = vec8(GREEN3/4.0f,SAND);
			}
			else
			{
				//terrain.color = vec4(1,0,0,0);
			}
		}


		// TEMP attempting horizontal displacement
		//u += mod(pattern*n.x+altitude*n.y, 0.02f);
		//v += mod(altitude*n.x-pattern*n.y, 0.02f);

		vec4 position = vec4( u,v, scale.z * n.altitude, 1.0 );

		if( show_grid && (i==0 || j==0) ) terrain.color *= 0.40f;
		
		rgba color_ff;
		color_ff.r = 255*terrain.color.r;
		color_ff.g = 255*terrain.color.g;
		color_ff.b = 255*terrain.color.b;
		color_ff.a = 255*terrain.color.a;

		rgba mixmap_20;
		mixmap_20.r = 255*terrain.detail.r;
		mixmap_20.g = 255*terrain.detail.g;
		mixmap_20.b = 255*terrain.detail.b;
		mixmap_20.a = 255*terrain.detail.a;

		rgba mixmap_ff = mixmap_20;
		if( 1 )
		{
			float a = clamp( ff_4th_detail_decay * uniforms.in_size, 0.0f, 1.0f );
			vec4  w = vec4( 0.1, 0.9, 1.0, 0.0 );

			rgba m0 = mixmap_20;
			m0.r = 128 + m0.r/2;
			m0.g = 128 + m0.g/2;
			m0.b = 128 + (m0.b + m0.a)/2;  // since cobs is mixed in by modulation, it must be accompained by a neutral detail like sand

			rgba m1 = rgba(	128 + m0.r*w.r/2, 
							128 + m0.g*w.g/2, 
							128 + (255 - m0.r*w.r - m0.g*w.g)/2, 
							m0.a*w.a );

			mixmap_ff = mix( m0, m1, a );
		}

	}


};





