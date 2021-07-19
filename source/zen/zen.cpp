#include "gdevice.h"
#include <engines/graphic/renderer.h>
//#include "zen.h"

int main()
{
	Renderer* x = new Renderer();

	return GDevice().run();
}
