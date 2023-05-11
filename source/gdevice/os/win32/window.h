#pragma once

#include "os/window.h" // ??
#include "os/keyboard.h"
#include "os/listener.h"
#include "os/win32/timer.h" // ??
#include "type/glsl.h"
#include "gpu/renderer.h"

//#include "os/log.h"


template<class Application>
class Window //: public Window
{
public:
	Listener<Application>* listener;
	Renderer* renderer;

    bool isActive;
    //bool isIconified;
	bool fullscreen;
	vec2 windowedSize;
	vec2 size;
	int bits;

	bool isPointerVisible;
	int mouseX;
	int mouseY;
	int mouseWheel;

	float mouseDeltaX()
	{
		static int previousMouseX = mouseX;
		int delta = mouseX - previousMouseX;
		previousMouseX = mouseX;
		return float(delta)/size.width;
	}

	float mouseDeltaY()
	{
		static int previousMouseY = mouseY;
		int delta = mouseY - previousMouseY;
		previousMouseY = mouseY;
		return float(delta)/size.height;
	}

	int mouseDeltaWheel()
	{
		static int previousWheel = mouseWheel;
		int delta = mouseWheel - previousWheel;
		previousWheel = mouseWheel;
		return delta;
	}

private:
	HWND	hWnd;
	HDC		hDC; // TEMP?
	HGLRC	hRC;

public:

    Window( Listener<Application>* listener = 0 ) : hWnd(NULL), hDC(NULL), hRC(NULL), renderer(NULL) //
	{
        DEBUG_ASSERT(!listener); // TODO Windows without listener
		this->listener = listener;
		renderer = new Renderer();
        renderer->controls.values = Key::getCounters(); // TEMP?
        isActive = false;
		fullscreen = false;
		isPointerVisible = false;
		windowedSize = vec2(WINDOW_WIDTH, WINDOW_HEIGHT); 
		bits = 32;
		size = vec2(0,0);
	}

	~Window()
	{
		destroy();
	}

	bool create()
	{
		PIXELFORMATDESCRIPTOR pfd = { sizeof(pfd), 1, 
			PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER, PFD_TYPE_RGBA,
			bits, 0,0,0,0,0,0,0,0, 0,0,0,0,0, 24, 8, 0, PFD_MAIN_PLANE, 0,0,0,0 };

		char* name = "gdevice";
		WNDCLASS windowClass = { CS_HREDRAW | CS_VREDRAW | CS_OWNDC, StaticWindowProc, 0,0,0,0,0,0,0, name };

        int	pixelFormat;
		if( !(RegisterClass(&windowClass))
         || !(hWnd = CreateWindow(name,0,0,0,0,windowedSize.width,windowedSize.height,0,0,0,0))
		 || !(hDC = GetDC(hWnd))
		 || !(pixelFormat = ChoosePixelFormat(hDC, &pfd))
		 || !(SetPixelFormat(hDC, pixelFormat, &pfd))
		 || !(hRC = wglCreateContext(hDC))
		 || !(wglMakeCurrent(hDC, hRC)))
		{
            DEBUG_WARNING(GetLastErrorAsString());
			destroy();
			return false;
		}	

        SetWindowLongPtr( hWnd, GWLP_USERDATA, (LONG_PTR) this );
		ShowWindow(hWnd, SW_NORMAL);
		ShowCursor( FALSE ); 

        // Silent crash (no dialog)
        SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
        _set_abort_behavior(0, _WRITE_ABORT_MSG);

        DEBUG_ASSERT(hRC);
		return true;
	}

    int runPlainMessageLoop()
    {
        if(create()) {
            MSG msg;
		    while(GetMessage(&msg, hWnd, 0, 0) > 0) {
			    if(msg.message == WM_CLOSE) 
                    return 0;

                if(isActive) {
                    TranslateMessage(&msg);   // Translates virtual key codes into WM_CHAR messages.
			        DispatchMessage(&msg);    // Dispatch the message to the WindowProc thread.
                }
		    }
        }
        return 0;
    }

	void show()
	{
		if( fullscreen ) {
			goFullscreen();
		} else {
			goWindowed();
		}
	}

	void goFullscreen()
	{
		DEVMODE screen;
		for( int mode=0; EnumDisplaySettings(NULL, mode, &screen); mode++ );
		
		ChangeDisplaySettings( &screen, CDS_FULLSCREEN );
		windowedSize = size;
		size = vec2( screen.dmPelsWidth, screen.dmPelsHeight );

		//SetWindowLongPtr( hWnd, GWL_EXSTYLE, WS_EX_APPWINDOW | WS_EX_TOPMOST );
		SetWindowLongPtr( hWnd, GWL_STYLE, WS_POPUP | WS_VISIBLE );
		SetWindowPos(     hWnd, HWND_TOPMOST, 0, 0, size.width, size.height, SWP_SHOWWINDOW );
		//isChangeSuccessful = ChangeDisplaySettings(&fullscreenSettings, CDS_FULLSCREEN) == DISP_CHANGE_SUCCESSFUL;
		//	ChangeDisplaySettings( &screen, CDS_FULLSCREEN );
		ShowWindow( hWnd, SW_MAXIMIZE );
	}

	void goWindowed()
	{
		//ChangeDisplaySettings( NULL, 0 );

		SetWindowLongPtr( hWnd, GWL_STYLE, WS_OVERLAPPEDWINDOW | WS_VISIBLE );
		//ChangeDisplaySettings( NULL, CDS_RESET );

		RECT rc = { 0, 0, windowedSize.width, windowedSize.height };
		AdjustWindowRect( &rc, WS_OVERLAPPEDWINDOW, FALSE );   
		size = vec2( rc.right-rc.left, rc.bottom-rc.top );

		SetWindowPos( hWnd, HWND_NOTOPMOST, 0, 0, size.width, size.height, SWP_SHOWWINDOW );

		ChangeDisplaySettings( NULL, 0 ); //
		ShowWindow( hWnd, SW_RESTORE );
	}

	void togglePointer()
	{
		setPointerVisibility( !isPointerVisible );
	}

	void setPointerVisibility( bool visibility )
	{
		this->isPointerVisible = visibility;
		ShowCursor( visibility );
	}


	void setListener( Listener<Application>* listener )
	{
		this->listener = listener;
	}

	void setTitle( const char *format, ... )
	{
		char buffer[1024];

		if(format == NULL) {
			strcpy(buffer, "") ;
		} else {
			va_list ap;
			va_start(ap, format);
				vsprintf(buffer, format, ap);
			va_end(ap);
		}

		SetWindowText( hWnd, buffer );
	}

	void destroy()
	{
		if( hWnd ) {
			if( hDC ) {
				wglMakeCurrent( hDC, 0 );
				if( hRC ) {
					wglDeleteContext( hRC );
				}
				ReleaseDC( hWnd, hDC );
			}
			DestroyWindow( hWnd );
		}

		if( fullscreen ) {
			ChangeDisplaySettings( NULL, 0 );
			ShowCursor( TRUE );
		}
	}

	static LRESULT CALLBACK StaticWindowProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
	{
		Window* window = (Window*) GetWindowLongPtr( hWnd, GWLP_USERDATA );
		return window->WindowProc( hWnd, uMsg, wParam, lParam );
	}

	LRESULT CALLBACK WindowProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
	{
		switch( uMsg )
		{
			
            case WM_NCCREATE:
			case WM_CREATE:
                break;

			 case WM_ACTIVATE: 
                isActive = LOWORD(wParam) != WA_INACTIVE;
				//isIconified = HIWORD(wParam);
				break;

            case WM_SIZE:
				switch (wParam)
				{
					case SIZE_MINIMIZED:
						return 0;

					case SIZE_MAXIMIZED:
					case SIZE_RESTORED:
					{
                        DEBUG_RUN_ONCE(
                            renderer->initialize();
			                listener->onOpen(*this);
                        );
                        
						size = vec2( LOWORD(lParam), HIWORD(lParam) );
						listener->onSize( *this );
					}
					return 0;
				}
				break;

			case WM_PAINT:
                DEBUG_CHECKPOINT_ONCE(WM_PAINT);
                DEBUG_ASSERT(listener);
				{
                    static Timer timer;
                    listener->onDraw( *this, timer.elapsed() );
                    HDC hdc = GetDC(hWnd);
				    SwapBuffers( hDC /*GetDC(hWnd)*/ );
				}
				return 0;

            case WM_ERASEBKGND:				// remove. useless since HBRUSH is null in WNDCLASS
				return 1;					// Prevent flickering while resizing

			case WM_DISPLAYCHANGE:
				// TODO
				break;

		   case WM_DESTROY:
				// TODO
				break;

			case WM_CLOSE:
				//PostMessage( hWnd, WM_CLOSE, 0, 0 ); // ?
				return 0;

            //////////////////////////
            // Controls
            case WM_KEYDOWN:
			case WM_SYSKEYDOWN:
				break;

			case WM_KEYUP:
			case WM_SYSKEYUP:
				break;

			case WM_CHAR:
				break;

			case WM_LBUTTONDBLCLK:
				fullscreen = !fullscreen;
				show();
				return 0;

			case WM_MOUSEMOVE:
				if( isPointerVisible ) {
					mouseX = LOWORD(lParam);
					mouseY = HIWORD(lParam);
				} else {
					RECT clientRect;
					GetClientRect( hWnd, &clientRect);

					POINT clientCenter;
					clientCenter.x = clientRect.right / 2;
					clientCenter.y = clientRect.bottom / 2;

					int dx = LOWORD(lParam) - clientCenter.x;
					int dy = HIWORD(lParam) - clientCenter.y;
					
					if( dx != 0 || dy != 0 ) {
						mouseX += dx;
						mouseY += dy;
						ClientToScreen( hWnd, &clientCenter );
						SetCursorPos( clientCenter.x, clientCenter.y );
					}
				}
				return 0;

			case WM_MOUSEWHEEL:
				mouseWheel += ((int)wParam) >> 16;
				return 0;
	
			case WM_MOVE:
				// TODO
				return 0;

			case WM_SYSCOMMAND:
                //DEBUG_CHECKPOINT(WM_SYSCOMMAND);
				switch( wParam ) {
					case SC_SCREENSAVE:     // Screensaver trying to start
					case SC_MONITORPOWER:   // Monitor entering powersave
						if( fullscreen ) {
							return 0;	    // Forbids when fullscreen
						}
						break;

					case SC_KEYMENU:        // User accessing menu using ALT
						return 0;
				}
				break;                      // Let to happen otherwise

		}

		return DefWindowProc (hWnd, uMsg, wParam, lParam);
	}

	void close()
	{
		PostMessage( hWnd, WM_CLOSE, 0, 0 );
	}

	void toogleFullscreen()
	{
		PostMessage( hWnd, WM_LBUTTONDBLCLK, 0, 0 );
	}

};
