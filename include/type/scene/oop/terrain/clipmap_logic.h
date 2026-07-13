#pragma once

#include <type/cpu.h>


struct ClipmapScroll
{
	vec2 position;
	int deltaI;
	int deltaJ;
	int quadrantInCoarserTile;
	bool invalidated;
};

inline ClipmapScroll CalculateClipmapScroll(
	double locationX,
	double locationY,
	double previousLocationX,
	double previousLocationY,
	float tileSize,
	int oddity )
{
	ClipmapScroll scroll;
	float span = 2.0f * tileSize;
	float scrollX = scrollValue(locationX, span);
	float scrollY = scrollValue(locationY, span);
	bool bx = scrollX >= tileSize;
	bool by = scrollY >= tileSize;

	scroll.position.x = scrollX / tileSize - float(oddity);
	scroll.position.y = scrollY / tileSize - float(oddity);
	scroll.deltaI = -int((tileValue(locationX, span) -
		tileValue(previousLocationX, span)) * 2.0f);
	scroll.deltaJ = -int((tileValue(locationY, span) -
		tileValue(previousLocationY, span)) * 2.0f);
	scroll.quadrantInCoarserTile = int(bx) + int(by) * 2;
	scroll.invalidated = scroll.deltaI != 0 || scroll.deltaJ != 0;
	return scroll;
}
