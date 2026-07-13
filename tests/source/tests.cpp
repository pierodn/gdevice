#include <stdio.h>
#include <string.h>

#include <gdevice.h>


int failures = 0;

void Expect(bool condition, const char* expression, int line)
{
	if (!condition)
	{
		printf("Failed at line %d: %s\n", line, expression);
		failures++;
	}
}

#define EXPECT(expression) Expect((expression), #expression, __LINE__)

void TestPublicApi()
{
	EXPECT(strcmp(GdeviceName(), "gdevice") == 0);
}

void TestClipmapScroll()
{
	ClipmapScroll scroll = CalculateClipmapScroll(1.0, 0.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.position.x == 0.0f);
	EXPECT(scroll.position.y == -1.0f);
	EXPECT(scroll.deltaI == 0);
	EXPECT(scroll.deltaJ == 0);
	EXPECT(scroll.quadrantInCoarserTile == 1);
	EXPECT(!scroll.invalidated);

	scroll = CalculateClipmapScroll(2.0, 0.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.position.x == -1.0f);
	EXPECT(scroll.deltaI == -2);
	EXPECT(scroll.invalidated);

	scroll = CalculateClipmapScroll(-0.1, -2.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.position.x == 0.9f);
	EXPECT(scroll.position.y == -1.0f);
	EXPECT(scroll.deltaI == 2);
	EXPECT(scroll.deltaJ == 2);
	EXPECT(scroll.quadrantInCoarserTile == 1);
	EXPECT(scroll.invalidated);
}

void TestClipmapQuadrants()
{
	ClipmapScroll scroll = CalculateClipmapScroll(0.0, 0.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.quadrantInCoarserTile == 0);

	scroll = CalculateClipmapScroll(1.0, 0.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.quadrantInCoarserTile == 1);

	scroll = CalculateClipmapScroll(0.0, 1.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.quadrantInCoarserTile == 2);

	scroll = CalculateClipmapScroll(1.0, 1.0, 0.0, 0.0, 1.0f, 1);
	EXPECT(scroll.quadrantInCoarserTile == 3);
}

void TestCpuMath()
{
	EXPECT(scrollValue(-0.1, 2.0) == 1.9f);
	EXPECT(tileValue(-0.1, 2.0) == -1.0f);
	EXPECT(tileValue(3.9, 2.0) == 1.0f);
	EXPECT(amod(-0.1, 2.0) == 1.9);
	EXPECT(average(rgba(0, 10, 20, 30), rgba(10, 20, 30, 40)) == rgba(5, 15, 25, 35));
}

int main()
{
	TestPublicApi();
	TestClipmapScroll();
	TestClipmapQuadrants();
	TestCpuMath();

	if (failures != 0)
	{
		printf("%d test(s) failed.\n", failures);
		return 1;
	}

	printf("Tests passed.\n");
	return 0;
}
