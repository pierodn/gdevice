#pragma once

#include <gl/link.h>
#include <gl/buffer.h>
#include <vector>

#define GLSL_(source) #source
#define GLSL(type, version, source) "\n" #type ":\n#version " #version "\n" #source "\n"

typedef Buffer Shader;

struct Program : public std::vector<Shader>, Buffer
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
    void Build(const std::string& source)
    {
        ASSERT(source.size() < MAX_PROGRAM_SOURCE_LENGTH);
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
    int BuildShaders(const std::string& source, int* shaders)
    {
        DEBUG_PATH;

        const char* ShaderTypeNames[] = { "VERTEX", "CONTROL", "EVALUATION", "GEOMETRY", "FRAGMENT", "COMPUTE" };
        const int   ShaderTypeCodes[] = { GL_VERTEX_SHADER, GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER, GL_COMPUTE_SHADER };
        int ShadersCount = sizeof(ShaderTypeCodes)/sizeof(int);

        int shaderIndex = 0;
        for( int i=0; i<ShadersCount; i++ )
        {
            size_t shaderPosition = source.find(ShaderTypeNames[i]);
            if(shaderPosition == std::string::npos)
            {
                continue;
            }

            shaderPosition += strlen(ShaderTypeNames[i]);
            DEBUG_ASSERT(source[shaderPosition] == ':');
            shaderPosition++;

            // Find next shader (if any) to compute the current shader length in characters.
            size_t shaderLenght = std::string::npos;
	        for( int j=i+1; j<ShadersCount; j++ )
	        {
                size_t nextShaderPosition = source.find(ShaderTypeNames[j]);
                if(nextShaderPosition != std::string::npos)
                {
		            shaderLenght = nextShaderPosition - shaderPosition;
		            break;
                }
	        }

            std::string shaderSource = source.substr(shaderPosition, shaderLenght);
 
            color(CMD_YELLOW, 0); DEBUG_PRINT("%s ", ShaderTypeNames[i]);

	        int shaderId = CompileShader(shaderSource, ShaderTypeCodes[i]);
            DEBUG_ASSERT(shaderId);

            shaders[shaderIndex++] = shaderId;
        }
        DEBUG_PRINT("\n");
        ASSERT(shaderIndex > 0);
        shaders[shaderIndex] = 0;
        
        return shaderIndex;
    }

    void PrintSource(const std::string& source, int focusLineNumber = -1)
	{
        ASSERT(source.size() <= MAX_PROGRAM_SOURCE_LENGTH);

        size_t p = 0;
        for(int lineNumber = 2; p<=source.size(); lineNumber++) 
		{
            size_t nIndex = p;
            while(nIndex < source.size() && source.at(nIndex) != '\n') nIndex++;
            if(nIndex >= source.size()) {
                break;
            }
            
            if((focusLineNumber < 0)
                || (focusLineNumber >= 0) && (focusLineNumber-6 <= lineNumber) && (lineNumber <= focusLineNumber+3)) 
            { 
                std::string lineText = source.substr(p, nIndex-p);

                if(lineNumber == focusLineNumber) {
                    color(CMD_WHITE, CMD_RED); 
                } else {
                    color(CMD_LIGHTGRAY, 0); 
                }
                printf("%i: %s\n", lineNumber, lineText.c_str());
            }

	        p = nIndex + 1;
        }
    }

    int CompileShader(const std::string& source, int type)
    {
        GLint shaderId = glCreateShader( type );
        const char* shader_cstr = source.c_str();
        int shader_length = source.size();
        glShaderSource(shaderId, 1, &shader_cstr, &shader_length );
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
	            PrintSource( source, line ); 
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
