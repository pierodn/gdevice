#pragma once

#include <gl/link.h>
#include <gl/cacheable.h>
#include <vector>

#define GLSL_(source) #source
#define GLSL(type, version, source) "\n" #type ":\n#version " #version "\n" #source "\n"

typedef Cacheable Shader;

struct Program : public std::vector<Shader>, Cacheable
{
    Program()
	{
        id = 0;
	}

	~Program() 
	{ 
		for( unsigned int i = 0; i < size(); i++ )
		{
			at(i).deallocate();
            //delete at(i);
		}
	}

private:
    static void deallocateProgram( uint& id )
    {
        if( id ) glDeleteProgram(id);
        id = 0;
    }
    static void deallocateShader( uint& id )
    {
        if( id ) glDeleteShader(id);
        id = 0;
    }

public:
    void Build(char* source)
    {
        ASSERT(source != NULL);
        ASSERT(strlen(source) < MAX_PROGRAM_SOURCE_LENGTH);
        ASSERT(id == 0);

        id = glCreateProgram();
        deallocator = &Program::deallocateProgram;

        const int ShaderTypesCount = 10;
        int shaders[ShaderTypesCount];
        int shaderCount = BuildShaders(source, shaders);
        
        for(int i=0; i<shaderCount; i++)
        {
            ASSERT(shaders[i]);

            Shader shader(shaders[i], &Program::deallocateShader);
            push_back(shader);

            AttachShaderToProgram(id, shader.id ); 
        }

        LinkProgram(id);
    }
    
    template<class T> 
    void SetUniform(char* name, T value)
    {
        ASSERT(id);
        // TODO ASSERT program is bound
        //ASSERT(GetCurrentProgramId() == id); // FIX
        SetUniform(id, name, value);
    }

    template<class T> 
    void SetUniform(const std::string& name, T value)
    {
        SetUniform(id, (char*)name.c_str(), value);
    }

    void Use()
    {
       glUseProgram(id);
    }

    void Run()
    {
        ASSERT(id);
// TODO
/*
        glPatchParameteri(GL_PATCH_VERTICES, 4);		
	    glDrawElements(
		    GL_PATCHES, 
		    ibo.lods[lod].count, // It is 1024 for tiles 17x17
		    GL_UNSIGNED_SHORT, 
		    0
	    ); 
*/
    }

private:
    int BuildShaders(const char* source, int* shaders)
    {
        DEBUG_PATH;

        const char* shaderTypes[] = { "VERTEX", "CONTROL", "EVALUATION", "GEOMETRY", "FRAGMENT", "COMPUTE" };
        const int   shaderTypeCodes[] = { GL_VERTEX_SHADER, GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER, GL_COMPUTE_SHADER };
        int shadersCount = sizeof(shaderTypeCodes)/sizeof(int);

        int shaderIndex = 0;
        for( int i=0; i<shadersCount; i++ )
        {
            char* shaderPosition = strstr2( (char*)source, shaderTypes[i]);
            if(shaderPosition == NULL) continue;
            
            shaderPosition += strlen(shaderTypes[i]);
            DEBUG_ASSERT(*shaderPosition == ':');
            shaderPosition++;

            // Find next shader (if any) to compute the current shader length in characters.
	        int shaderLenght = 0;
	        for( int j=i+1; j<shadersCount; j++ )
	        {
                char* nextShaderPosition = strstr2( (char*)source, shaderTypes[j]);
                if(nextShaderPosition != NULL)
                {
		            shaderLenght = nextShaderPosition - shaderPosition;
		            break;
                }
	        }

	        if( shaderLenght==0 )
	        {
		        shaderLenght = strlen(shaderPosition);
	        }
		
            color(CMD_YELLOW, 0); DEBUG_PRINT("%s ", shaderTypes[i]);

	        int shaderId = CompileShader(shaderPosition, shaderLenght, shaderTypeCodes[i]);
            DEBUG_ASSERT(shaderId);

            shaders[shaderIndex++] = shaderId;
        }
        DEBUG_PRINT("\n");
        ASSERT(shaderIndex > 0);
        shaders[shaderIndex] = 0;
        
        return shaderIndex;
    }

    void PrintSource(char* source, int length, int focusLine = -1)
	{
        ASSERT( length <= MAX_PROGRAM_SOURCE_LENGTH );

        char buffer[MAX_PROGRAM_SOURCE_LENGTH];
        strncpy(buffer, source, length);

        char* p = buffer;
        char lineText[256];
        for(int line = 2; *p; line++) 
		{
	        char* nIndex = strstr2(p, "\n");
	        if(nIndex == NULL) break;

            if((focusLine < 0) || (focusLine >= 0) && (focusLine-6 <= line) && (line <= focusLine+3)) 
            { 
	            int length = nIndex - p;			
	            strncpy(lineText, p, length);
	            lineText[length] = 0;

                if(line == focusLine) {
                    color(CMD_WHITE, CMD_RED); 
                } else {
                    color(CMD_LIGHTGRAY, 0); 
                }
	            printf("%i: %s\n", line, lineText);
            }

	        p = nIndex + 1;
        }
    }

    int CompileShader(char* source, int len, int type)
    {
        for(; *source>0 && *source<=' '; source++) len--;

        ASSERT(*source);
        ASSERT(len >= 0);

        GLint shaderId = glCreateShader( type );
        glShaderSource(shaderId, 1, (const char**)&source, &len );
        glCompileShader(shaderId);

        int successful = 0;
        glGetShaderiv( shaderId, GL_COMPILE_STATUS, &successful );

        if( !successful )
        {
            color(CMD_LIGHTRED,0); DEBUG_PRINT("\n\n");

            char message[512];
            glGetShaderInfoLog( shaderId, sizeof(message), NULL, message );
			printf("%s\n", message);
            //DEBUG_PRINT(message);

			// Parse to get the error line.
			// NOTE: Different drivers give different error messages
			// https://gamedev.stackexchange.com/questions/38685/is-this-a-reliable-method-of-parsing-glgetshaderinfolog
			int line = -1;
			if(line <= 0)
            {
				// NVIDIA error message
				char *p1 = 0, *p2 = 0;
				p1 = strchr(message, '(');
				if( p1>0 ) p2 = strchr(p1+1, ')');
				if( p1>0 && p2>0 ) line = strtol(p1+1, &p2, 10) + 1; // it was +0 
			}
			if(line <= 0)
            {
				// ATI/Intel error message
				char *p0 = 0, *p1 = 0, *p2 = 0;
				p0 = strchr(message, ':');
				if( p0>0 ) p1 = strchr(p0+1, ':');
				if( p1>0 ) p2 = strchr(p1+1, ':');
				if(p0>0 && p1>0 && p2>0) line = strtol(p1+1, &p2, 10) + 1;
			}
			ASSERT(line > 0);

            #if defined(DEBUG_SHOW_GLSL_SOURCE)
	            PrintSource( source, len, line ); 
                printf("\n");
            #endif

            CRITICAL("Shader compile error");
        }  	
        return shaderId;
    }

    void CheckProgram( int programId )
    {
        int successful;
        glGetProgramiv( programId, GL_LINK_STATUS, &successful);
        if( successful ) return;

        char message[512];
        glGetProgramInfoLog( programId, sizeof(message), 0, message );
        DEBUG_CRITICAL( message );
    }

    void AttachShaderToProgram(int programId, int shaderId)
    {
        glAttachShader(programId, shaderId);
    }

    void LinkProgram(int programId)
    {
	    glLinkProgram(programId);
        CheckProgram(programId);
    }

    int GetUniformLocation( int programId, char* name )
    {
        int location = glGetUniformLocation( programId, name );

        if( location<0 )
        {	
	        char buffer[200];
	        sprintf( buffer, "Uniform '%s' not found.\n", name );
	        DEBUG_CRITICAL( buffer );
        }
        return location;
    }
    void SetUniform( int programId, char* name, int value )
    {
        glUniform1i( GetUniformLocation(programId, name), value );
    }
    void SetUniform( int programId, char* name, float value )
    {
        glUniform1f( GetUniformLocation(programId, name), value );
    }
    void SetUniform( int programId, char* name, vec2 vector )
    {
        glUniform2fv( GetUniformLocation(programId, name), 1, vector.array );
    }
    void SetUniform( int programId, char* name, vec3& vector )
    {
        glUniform3fv( GetUniformLocation(programId, name), 1, vector.array );
    }
    void SetUniform( int programId, char* name, vec4& vector )
    {
        glUniform4fv( GetUniformLocation(programId, name), 1, vector.array );
    }
    void SetUniform( int programId, char* name, mat4& matrix )
    {
        glUniformMatrix4fv( GetUniformLocation(programId, name), 1, 0, matrix.array );
    }
    void SetUniform( int programId, char* name, mat3& matrix )
    {
        glUniformMatrix3fv( GetUniformLocation(programId, name), 1, 0, matrix.array );
    }
/*
    int GetCurrentProgramId()
    {
        GLint _currentProgram;
        glGetIntegerv(GL_CURRENT_PROGRAM, &_currentProgram);
        return _currentProgram;
    }
*/

};
