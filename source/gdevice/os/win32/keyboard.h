#pragma once

//////////////////////////////////////
// TODO: Remove Win32 dependency.
#if defined(_WIN32)
#include <windows.h>
#define LSHIFT		VK_LSHIFT
#define RSHIFT		VK_RSHIFT
#define LCONTROL	VK_LCONTROL
#define RCONTROL	VK_RCONTROL
#define ESCAPE		VK_ESCAPE
#define END			VK_END
#define HOME		VK_HOME
#define LEFT		VK_LEFT
#define UP			VK_UP
#define RIGHT		VK_RIGHT
#define	DOWN		VK_DOWN
#define INS 		VK_INSERT
#define DEL  		VK_DELETE
#define	NUMLOCK		VK_NUMLOCK
#define SCROLL		VK_SCROLL
#define NUMPAD0		VK_NUMPAD0
#define NUMPAD1		VK_NUMPAD1
#define NUMPAD2		VK_NUMPAD2
#define	NUMPAD3		VK_NUMPAD3
#define	NUMPAD4		VK_NUMPAD4
#define NUMPAD5		VK_NUMPAD5
#define NUMPAD6		VK_NUMPAD6
#define NUMPAD7		VK_NUMPAD7
#define NUMPAD8		VK_NUMPAD8
#define NUMPAD9		VK_NUMPAD9
#define F1			VK_F1
#define F2			VK_F2
#define	F3			VK_F3
#define F4			VK_F4
#define F5			VK_F5
#define F6			VK_F6
#define F7			VK_F7
#define F8			VK_F8
#define F9			VK_F9
#define	F10			VK_F10
#define F11			VK_F11
#define F12			VK_F12
#endif 
////////////////////////////////////


class Key
{
private:
    int key;
    static int counters[0xFF];

public:
    Key(int key) { this->key = key; }

    enum { RELEASED, JUST_PRESSED, STILL_PRESSED, JUST_RELEASED };

	static int getCode( int key )
    {
        static int state[0xFF];

	    int previous = state[key];
	    int current = GetAsyncKeyState( key ); // WIN32
	    state[key] = current;
	    int code = 
		    !previous && !current ? RELEASED :
		    !previous &&  current ? JUST_PRESSED :
		     previous &&  current ? STILL_PRESSED :
								    JUST_RELEASED; 
        if( code == JUST_PRESSED ) Key::counters[key]++;
	    return code;
    }

    bool isReleased()       { return getCode(key) == RELEASED; }
    bool isJustPressed()    { return getCode(key) == JUST_PRESSED; }
    bool isStillPressed()   { return getCode(key) == STILL_PRESSED; }
    bool isJustReleased()   { return getCode(key) == JUST_RELEASED; }
    bool isPressed()        { return getCode(key) != RELEASED; }

	static int  getCounter( int key ) { return counters[key]; }
    static int* getCounters() { return counters; }
	static int cursorDeltaX() { return getCode(RIGHT)? +1 : getCode(LEFT)? -1 : 0; }
	static int cursorDeltaY() { return getCode(DOWN) ? +1 : getCode(UP)  ? -1 : 0; }
};

int Key::counters[0xFF] = {};