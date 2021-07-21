#pragma once

#include "window.h"
#include "keyboard.h"
#include "listener.h"


template<class ConcreteApplication>
class Application : public Listener<ConcreteApplication>
{
 public:   
	Window<ConcreteApplication> window;

    int run()
    {
		while(window.isOpen() && !Key(ESCAPE).isPressed()) {
			doUpdate();
		}
		return 0;
	}

    Application()
	{
        window.setListener(this);
	}
/*
	Application( int width, int height )
	{
		window.width = width;
		window.height = height;
        window.setListener( this );
	}
*/
};
