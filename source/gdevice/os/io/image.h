#pragma once

#include "type/gpu.h"
#include "__temp/texture.h"

namespace IO
{
    const uint16 BM = uint16('M' << 8) | 'B';

	template<typename TEXEL>
	void load( Texture<TEXEL>& texture, char* filename )
	{
		FILE* file = fopen(filename, "rb");		
        DEBUG_ASSERT(file!=NULL);

		BITMAPFILEHEADER bmfh;
		fread( &bmfh, sizeof(BITMAPFILEHEADER), 1, file );
		DEBUG_ASSERT(bmfh.bfType==BM);

		BITMAPINFOHEADER bmih;
		fread( &bmih, sizeof(BITMAPINFOHEADER), 1, file );
		int pixelSize = bmih.biBitCount/8;
		int width = bmih.biWidth;
		int height = bmih.biHeight;
        DEBUG_ASSERT(pixelSize == sizeof(TEXEL));

		// There is a palette to skip with 8bit bitmaps.
		if( pixelSize==1 ) {
			fseek(file, 256 * sizeof(rgba), SEEK_CUR);
		}

		texture.resample<0>(width, height);
		DEBUG_ASSERT(texture.bytes() == fread( texture.array, 1, texture.bytes(), file));

        DEBUG_TRACE(filename);
	}

    template<typename TEXEL>
	void save( Texture<TEXEL>& texture, char* filename )
	{
		FILE* file = fopen(filename, "wb");		
        DEBUG_ASSERT(file!=NULL);

        long biSize = sizeof(BITMAPV4HEADER);
        long bfOffBits = sizeof(BITMAPFILEHEADER) + biSize;
        long biSizeImage = texture.bytes();
        long bfSize = bfOffBits + biSizeImage;

        BITMAPFILEHEADER bmfh;
        bmfh.bfType             = BM;
        bmfh.bfSize             = bfSize;
        bmfh.bfReserved1        = 0;
        bmfh.bfReserved2        = 0;
        bmfh.bfOffBits          = bfOffBits;
        
		BITMAPV4HEADER bmv4h;
        bmv4h.bV4Size           = biSize;
        bmv4h.bV4Width          = texture.size.x;
        bmv4h.bV4Height         = texture.size.y;
        bmv4h.bV4Planes         = 1;
        bmv4h.bV4BitCount       = sizeof(TEXEL)*8;
        bmv4h.bV4V4Compression  = sizeof(TEXEL)==4 ? BI_BITFIELDS : BI_RGB;
        bmv4h.bV4SizeImage      = biSizeImage;
        bmv4h.bV4XPelsPerMeter  = 2835;
        bmv4h.bV4YPelsPerMeter  = 2835;
        bmv4h.bV4ClrUsed        = 0;
        bmv4h.bV4ClrImportant   = 0;
        bmv4h.bV4RedMask        = 0x000000FF;
        bmv4h.bV4GreenMask      = 0x0000FF00;
        bmv4h.bV4BlueMask       = 0x00FF0000;
        bmv4h.bV4AlphaMask      = 0xFF000000;
        bmv4h.bV4CSType         = 0;
        //bmv4h.bV4Endpoints      = {{0,0,0}, {0,0,0}, {0,0,0}};
        //    for(int i = 0; i < 3 * 3; i++) {
        //        ((long*)(bmv4h.bV4Endpoints))[i] = 0;
        //    }
        bmv4h.bV4GammaRed       = 0;
        bmv4h.bV4GammaGreen     = 0;
        bmv4h.bV4GammaBlue      = 0;

        fwrite(&bmfh,  1, sizeof(BITMAPFILEHEADER), file);
        fwrite(&bmv4h, 1, sizeof(BITMAPV4HEADER),   file);
        fwrite(texture.array, 1,  texture.bytes(),  file);
        fclose(file);

        DEBUG_TRACE(filename);
	}

 /*
    template<typename TEXEL>
	void saveRGB( Texture<TEXEL>& texture, char* filename )
	{
		FILE* file = fopen(filename, "wb");		
        DEBUG_ASSERT(file!=NULL);

		BITMAPINFOHEADER bmih;
        bmih.biCompression = BI_RGB;
        bmih.biWidth = texture.size.x;
        bmih.biHeight = texture.size.y;
        bmih.biPlanes = 1;
        bmih.biBitCount = sizeof(TEXEL)*8;
        bmih.biSizeImage = texture.bytes();
        bmih.biSize = sizeof(BITMAPINFOHEADER);

        BITMAPFILEHEADER bmfh;
        bmfh.bfType = BM;
        bmfh.bfOffBits = sizeof(BITMAPFILEHEADER) + bmih.biSize;
        bmfh.bfSize = bmfh.bfOffBits + bmih.biSizeImage;
        bmfh.bfReserved1 = bmfh.bfReserved2 = 0;

        fwrite(&bmfh, 1, sizeof(BITMAPFILEHEADER), file);
        fwrite(&bmih, 1, sizeof(BITMAPINFOHEADER), file);
        fwrite(texture.array, 1,  texture.bytes(), file);
        fclose(file);

        DEBUG_TRACE(filename);
	}
*/

};