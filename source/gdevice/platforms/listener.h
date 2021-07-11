#pragma once

#include "platforms/window.h"

template<class Application>
class Window;

template <class Application>
struct Listener
{
	void onOpen( Window<Application>& window )
	{
		static_cast<Application*>(this)->onOpen( window );
	}

	void onDraw( Window<Application>& window, double elapsed )
	{
		static_cast<Application*>(this)->onDraw( window, elapsed );
	}

	void onSize( Window<Application>& window )
	{
		static_cast<Application*>(this)->onSize( window );
	}
};