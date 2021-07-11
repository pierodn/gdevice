#include <stdio.h>
#include <string.h>

#ifdef _MSC_VER
	#pragma warning(disable: 4996)
#endif

// PROP discard // or # lines
// PROP split the one in more appended strings separated by "\n"

int main(int argc, char **argv) 
{

	if( argc<3 )
	{
		fprintf( stderr, "Usage: bin2c <inputfile> <outputfile>\n");
		return -1;
	}

	FILE *input = fopen( argv[1], "rb");
	if( input==NULL )
	{
		fprintf( stderr, "Error: can't read the input file: %s\n", argv[1]);
		return -1;
	}

	FILE *output = fopen( argv[2], "wb");
	if( output==NULL )
	{
		fprintf( stderr, "Error: can't create the output file: %s\n", argv[2]);
		return -1;
	}

	char* varname = argv[1];
	for( char* t = varname-1; t = strpbrk(t+1,"/\\"); varname = t+1 );
	for( char* t = varname-1; t = strpbrk(t+1," !#$%&'()+,-.;=@[]^{}~"); *t = '_' );

	fprintf( output, "char %s[] = {", varname );
	for( int c; (c=fgetc(input))!=EOF; fprintf(output,
		c>=' ' && c<='~' && c!='\\' && c!='\'' ? "'%c'," : "%d,", (unsigned)c ) ); 
	fprintf( output, "0 };" );
	fclose(output);

	fclose(input);

	fprintf( stdout, "Symbol '%s' generated\n", varname);

	{fgetc(stdin);} // TEMP
	return 0;
}
