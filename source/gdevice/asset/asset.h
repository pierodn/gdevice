#pragma once

#include "os/io/image.h"
//#include "type/interpolation.h"

#define TEXTURE_PATH "../../asset/textures/"
#define MAX_MIPMAPS 10

/////////////////////////////////////////////////////////
// HOW TO MAKE TEXTURE HEIGHTMAPS OUT OF NORMAL TEXTURES
//
// 1. Find a 1024x1024 seamless tileable texture
//    Ideally with no shadows but with fine ambient.
// 2. Open it with IrfanView
// 3. Auto adjust colors (SHIFT+U)
// 4. Sharpen (optional) (SHIFT+S)
// 5. Color corrections (SHIFT+G): set gamma to 0.45
// 6. Covert to grayscale (CTRL+G) [ALWAYS right before saving]
// 7. Save as BMP (CTRL+S)
/////////////////////////////////////////////////////////

Texture<vec4> details;
Texture<vec4> detailsDx[MAX_MIPMAPS]; 
Texture<vec4> detailsDy[MAX_MIPMAPS];

void initializeTerrainDetails()
{
    Texture<byte> rock = Texture<byte>(TEXTURE_PATH "terrain/rock16b.bmp"); // rock1b
    Texture<byte> grit = Texture<byte>(TEXTURE_PATH "terrain/asphalt.bmp"); // cobs1b soil01b cobs12b
	Texture<byte> bone = Texture<byte>(TEXTURE_PATH "terrain/rock1b.bmp");  // rock1b rock11
    Texture<byte> sand = Texture<byte>(TEXTURE_PATH "terrain/sand1c.bmp");  // sand1c sand11

    DEBUG_ASSERT(rock.size.x == rock.size.y);
	int size = rock.size.x; 

    // Output file for external tools.
    Texture<rgba> output;
    output.resample<0>(size, size);
    for(int y = 0; y<size; y++)
	for(int x = 0; x<size; x++) {
		output.at(x,y) = rgba(rock.at(x,y), grit.at(x,y), bone.at(x,y), sand.at(x,y));
	}
    output.store(TEXTURE_PATH "terrain/output.bmp");

    // Encode the 4 details into one RGBA texture.
	details.resample<0>(size, size);
	for(int y = 0; y<size; y++)
	for(int x = 0; x<size; x++) {
		details.at(x,y) = vec4(rock.at(x,y), grit.at(x,y), bone.at(x,y), sand.at(x,y))/255.0f;
	}

    //
    // Grandient mipmaps are computed as FBM of the texture mipmaps,
    //

	int lods = log2(float(size));
	DEBUG_ASSERT(lods <= MAX_MIPMAPS);

    // Compute the texture mipmaps.
    Texture<vec4> mipmap[MAX_MIPMAPS];
	for(int lod = 0; lod<lods; lod++) {
		int size = 1<<(lods-lod);
        mipmap[lod].resample<0>(size, size);
		if(lod == 0) {  // Lod 0 => from files.
			for(int y=0; y<size; y++) 
			for(int x=0; x<size; x++) {
				mipmap[lod].at(x,y) = vec4(rock.at(x,y), grit.at(x,y), bone.at(x,y), sand.at(x,y))/255.0f;
			}
		} else {        // Lod>0 => from previous lod.
            for(int y=0; y<size; y++) 
			for(int x=0; x<size; x++)
				mipmap[lod].at(x,y) = mipmap[lod-1].box2x2(x*2, y*2);
		}
    }

    // Compute the gradient mipmaps as FBM of the texture mipmaps
	for(int lod = 0; lod<lods; lod++) {
		int size = 1<<(lods-lod);
		detailsDx[lod].resample<0>(size, size);
		detailsDy[lod].resample<0>(size, size);
		for(int y=0; y<size; y++)
		for(int x=0; x<size; x++) {
            // Fractional Brownian Motion
            vec4 dx = vec4(0.0f);
            vec4 dy = vec4(0.0f);
            for(int i = lod; i<lods; i++) {
                float s = 1<<(lods-i);
                float xx = x*s/size + 0.5;
                float yy = y*s/size + 0.5;
#if 0
                dx += mipmap[i].centralDx(xx,yy); 
                dy += mipmap[i].centralDy(xx,yy);
#elif 0
                dx += mipmap[i].prewittDx(xx,yy);
                dy += mipmap[i].prewittDy(xx,yy);
#else
                dx += mipmap[i].sobelDx(xx,yy);
                dy += mipmap[i].sobelDy(xx,yy);
#endif 
            }
			detailsDx[lod].at(x,y) = dx;
			detailsDy[lod].at(x,y) = dy;
		}
	}

}


