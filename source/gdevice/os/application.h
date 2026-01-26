#pragma once

#include "window.h"
#include "keyboard.h"
#include "listener.h"


namespace gd {
    template<class ConcreteApplication> class Application;
}

template<class ConcreteApplication>
class gd::Application : public Listener<ConcreteApplication>
{
 public:   
	Window<ConcreteApplication> window;
/*
    int run()
    {
        while(window.isOpen()) 
        {
			doUpdate();
		}
		return 0;
	}
*/
    int run() 
    {
        return window.runPlainMessageLoop();
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
