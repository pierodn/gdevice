#include "gdevice.h"
#include <gpu/renderer.h>
//#include "zen.h"

int main()
{
	Renderer* x = new Renderer();

	return GDevice().run();
}
