#pragma once

// TODO Linux, MacOS
#include <windows.h>
#include <stdio.h>		// vsprintf_s, va_list
#pragma comment( lib, "opengl32.lib" ) // wglDeleteContext


class Window;

struct Listener
{
	virtual void OnOpen( Window& window ) = 0;
	virtual void OnDraw( Window& window ) = 0;
	virtual void OnSize( Window& window ) = 0;
};

class Window
{
public:
	Listener* listener;

    bool isActive;
    //bool isIconified;
	bool fullscreen;
    struct Size {
	    long x;
	    long y;
    } windowSize, clientSize;
	int bits;

	bool isPointerVisible;
	int mouseX;
	int mouseY;
	int mouseWheel;

private:
	HWND	hWnd;
	HDC		hDC;
	HGLRC	hRC;

public:

    Window(long width = 640, long height = 480) : hWnd(NULL), hDC(NULL), hRC(NULL)
	{
		this->listener = listener;
        isActive = false;
		fullscreen = false;
		isPointerVisible = false;
        windowSize.x = width;
        windowSize.y = height;
        clientSize.x = 0;
        clientSize.y = 0;
		bits = 32;	
	}

	~Window()
	{
		Destroy();
	}

    void SetListener(Listener* listener)
	{
		this->listener = listener;
	}

	bool Open()
	{
        const bool resizeable = true;
        const char* name = "gdevice";

		PIXELFORMATDESCRIPTOR pfd = { sizeof(pfd), 1, 
			PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER, PFD_TYPE_RGBA,
			bits, 0,0,0,0,0,0,0,0, 0,0,0,0,0, 24, 8, 0, PFD_MAIN_PLANE, 0,0,0,0 };

		WNDCLASS windowClass = { CS_HREDRAW | CS_VREDRAW | CS_OWNDC, StaticWindowProc, 0,0,0,0,0,0,0, name };

        DWORD dwStyle = 0;
        if(resizeable) dwStyle += WS_OVERLAPPEDWINDOW;

        int	pixelFormat;
		if( !(RegisterClass(&windowClass))
         || !(hWnd = CreateWindow(name,0,dwStyle,0,0,windowSize.x,windowSize.y,0,0,0,0))
		 || !(hDC = GetDC(hWnd))
		 || !(pixelFormat = ChoosePixelFormat(hDC, &pfd))
		 || !(SetPixelFormat(hDC, pixelFormat, &pfd))
		 || !(hRC = wglCreateContext(hDC))
		 || !(wglMakeCurrent(hDC, hRC)))
		{
            DWORD errorMessageID = ::GetLastError();
            if(errorMessageID != 0) 
            {
                char* messageBuffer = "";
                size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                    NULL, errorMessageID, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);
                printf("Error: %.*s\n", size, messageBuffer);
                LocalFree(messageBuffer);
            }

			Destroy();
			return false;
		}	

        SetWindowLongPtr( hWnd, GWLP_USERDATA, (LONG_PTR) this );
		ShowWindow(hWnd, SW_NORMAL);
		ShowCursor( FALSE ); 

        // Silent crash (no dialog)
        SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
        _set_abort_behavior(0, _WRITE_ABORT_MSG);
		return true;
	}

    int RunDefaultMessageLoop()
    {
        if(Open()) 
        {
            MSG msg;
		    while((GetMessage(&msg, hWnd, 0, 0) > 0) && (msg.message != WM_CLOSE)) 
            {
                if(isActive) 
                {
                    TranslateMessage(&msg);   // Translates virtual key codes into WM_CHAR messages.
			        DispatchMessage(&msg);    // Dispatch the message to the WindowProc thread.
                }
		    }
        }
        return 0;
    }

	void Show()
	{
		if(fullscreen)
        {
			ShowFullscreen();
		} 
        else 
        {
			ShowWindowed();
		}
	}

	void ShowFullscreen()
	{
		DEVMODE screen;
		for( int mode=0; EnumDisplaySettings(NULL, mode, &screen); mode++ );
		
		ChangeDisplaySettings( &screen, CDS_FULLSCREEN );
		windowSize = clientSize;
		clientSize.x = screen.dmPelsWidth;
        clientSize.y = screen.dmPelsHeight;

		//SetWindowLongPtr( hWnd, GWL_EXSTYLE, WS_EX_APPWINDOW | WS_EX_TOPMOST );
		SetWindowLongPtr( hWnd, GWL_STYLE, WS_POPUP | WS_VISIBLE );
		SetWindowPos(     hWnd, HWND_TOPMOST, 0, 0, clientSize.x, clientSize.y, SWP_SHOWWINDOW );
		//isChangeSuccessful = ChangeDisplaySettings(&fullscreenSettings, CDS_FULLSCREEN) == DISP_CHANGE_SUCCESSFUL;
		//	ChangeDisplaySettings( &screen, CDS_FULLSCREEN );
		ShowWindow( hWnd, SW_MAXIMIZE );
	}

	void ShowWindowed()
	{
		//ChangeDisplaySettings( NULL, 0 );

		SetWindowLongPtr( hWnd, GWL_STYLE, WS_OVERLAPPEDWINDOW | WS_VISIBLE );
		//ChangeDisplaySettings( NULL, CDS_RESET );

		RECT rc = { 0, 0, windowSize.x, windowSize.y };
		AdjustWindowRect( &rc, WS_OVERLAPPEDWINDOW, FALSE );   
		clientSize.x = rc.right-rc.left;
        clientSize.y = rc.bottom-rc.top;

		SetWindowPos( hWnd, HWND_NOTOPMOST, 0, 0, clientSize.x, clientSize.y, SWP_SHOWWINDOW );

		ChangeDisplaySettings( NULL, 0 );
		ShowWindow( hWnd, SW_RESTORE );
	}

	void TogglePointer()
	{
		SetPointerVisibility( !isPointerVisible );
	}

	void SetPointerVisibility( bool visibility )
	{
		this->isPointerVisible = visibility;
		ShowCursor( visibility );
	}

	void SetTitle(const char *format, ...)
	{
		char buffer[1024] = "";

		va_list ap;
		va_start(ap, format);
			vsprintf_s(buffer, format, ap);
		va_end(ap);

		SetWindowText(hWnd, buffer);
	}

	void Destroy()
	{
		if(hWnd) 
        {
			if(hDC) 
            {
				wglMakeCurrent(hDC, 0);
				if(hRC)
                {
					wglDeleteContext(hRC);
				}
				ReleaseDC(hWnd, hDC);
			}
			DestroyWindow(hWnd);
		}

		if(fullscreen)
        {
			ChangeDisplaySettings(NULL, 0);
			ShowCursor(TRUE);
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
                        static bool hasRanOnce = false;
                        if( !hasRanOnce )
                        {
			                listener->OnOpen(*this);
                            hasRanOnce = true;
                        }
                        
						clientSize.x = LOWORD(lParam);
                        clientSize.y = HIWORD(lParam);
						listener->OnSize(*this);
					}
					return 0;
				}
				break;

			case WM_PAINT:
				{
                    listener->OnDraw(*this);
                    HDC hdc = GetDC(hWnd);
				    SwapBuffers(hDC);
				}
				return 0;

			case WM_DISPLAYCHANGE:
				break;

		    case WM_DESTROY:
				break;

			case WM_CLOSE:
				PostMessage(hWnd, WM_CLOSE, 0, 0);
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
				Show();
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

                    // TODO Found bugged delta when screen resolution is 4K.
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
				return 0;

			case WM_SYSCOMMAND:
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

    float GetMouseDeltaX()
	{
		static int previousMouseX = mouseX;
		int delta = mouseX - previousMouseX;
		previousMouseX = mouseX;
		return float(delta)/clientSize.x;
	}

	float GetMouseDeltaY()
	{
		static int previousMouseY = mouseY;
		int delta = mouseY - previousMouseY;
		previousMouseY = mouseY;
		return float(delta)/clientSize.y;
	}

	int GetMouseDeltaWheel()
	{
		static int previousWheel = mouseWheel;
		int delta = mouseWheel - previousWheel;
		previousWheel = mouseWheel;
		return delta;
	}

	void Close()
	{
		PostMessage(hWnd, WM_CLOSE, 0, 0);
	}

	void ToggleFullscreen()
	{
		PostMessage(hWnd, WM_LBUTTONDBLCLK, 0, 0);
	}

};
