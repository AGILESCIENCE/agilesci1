/// /////////////////////////////////////////////
///
/// AG_genapp automatically generated source file


#include <cstdio>
#include <cstring>
#include <string>
#include <vector>


using namespace std;



struct PilDesc
{
	char   type;
	string identifier;
	string prompt;
};

typedef vector<PilDesc> PilDescArr;

static PilDescArr* ReadParams(FILE* i)
{
PilDescArr* parArr = new PilDescArr;
char line[1024];
char lineCopy[1024];
while (fgets(line, 1024, i)) {
	strcpy(lineCopy, line);
	char* endline = strchr(line, '\n');
	if (endline)
		*endline = 0;
	char* token[7];
	token[0] = line;
	bool ok = true;
	for (int i=1; i<7 && ok; ++i) {
		char* comma = strchr(token[i-1], ',');
		if (!comma)
			ok = false;
		else {
			*comma = 0;
			token[i] = comma+1;
			}
		}
	if (!ok) {
		fprintf(stderr, "WARNING: Skipping illegal line: %s", lineCopy);
		continue;
		}
	char* prompt = token[6];
	int promptLen = strlen(prompt);
	if (promptLen<3 || prompt[0]!='\"' || prompt[promptLen-1]!='\"') {
		fprintf(stderr, "WARNING: Prompt should be a quoted string: %s", lineCopy);
		continue;
		}
	prompt[promptLen-1] = 0;
	prompt++;
	char type = *token[1];
	if (type!='b' && type!='i' && type!='r' && type!='s') {
		fprintf(stderr, "WARNING: Parameter type '%s' not supported\n", token[1]);
		continue;
		}

	PilDesc pd;
	pd.type = type;
	pd.identifier = string(token[0]);
	pd.prompt = string(prompt);
	parArr->push_back(pd);
	}
return parArr;
}

static void WriteSource(const PilDescArr& parArr, FILE* o)
{
fprintf(o, "/// /////////////////////////////////////////////\n///\n");
fprintf(o, "/// AG_genapp automatically generated source file\n\n\n");
fprintf(o, "#include \"PilParams.h\"\n\n\n");
fprintf(o, "const PilDescription c_params[] = {\n");
int count = parArr.size();
for (int i=0; i<count; ++i) {
	const PilDesc& pilDesc = parArr[i];
	const char* typeName = "";
	switch (pilDesc.type) {
		case 'b': typeName = "PilBool"; break;
		case 'i': typeName = "PilInt"; break;
		case 'r': typeName = "PilReal"; break;
		case 's': typeName = "PilString";
		}
	fprintf(o, "\t{ %s, \"%s\", \"%s\" },\n", typeName, pilDesc.identifier.c_str(), pilDesc.prompt.c_str());
	}

fprintf(o, "\t{ PilNone, \"\", \"\" }\n\t};\n\n\n");
fprintf(o, "class AppParams: public PilParams\n");
fprintf(o, "{\n");
fprintf(o, "public:\n");
fprintf(o, "\tAppParams(): PilParams(c_params) {}\n");
fprintf(o, "\tvoid Print() { PilParams::Print(c_params); }\n");
fprintf(o, "\t};\n\n\n");

fprintf(o, "int main(int argC, char *argV[])\n{\n");
fprintf(o, "AppParams appParams;\n");
fprintf(o, "if (!appParams.Load(argC, argV))\n");
fprintf(o, "\treturn -1;\n");
fprintf(o, "appParams.Print();\n\n");

for (int i=0; i<count; ++i) {
	switch (parArr[i].type) {
		case 'b': fprintf(o, "bool "); break;
		case 'i': fprintf(o, "int "); break;
		case 'r': fprintf(o, "double "); break;
		case 's': fprintf(o, "const char* ");
		}
	const string& id = parArr[i].identifier;
	fprintf(o, "%s = appParams[\"%s\"];\n", id.c_str(), id.c_str());
	}

fprintf(o, "\n/// Write your app code here\n\n\n");

fprintf(o, "return 0;\n}\n");
}



static void WriteHtml(const char* name, const PilDescArr& parArr, FILE* o)
{
fprintf(o, "<html>\n<head>\n");
fprintf(o, "<title>%s parameters</title>\n", name);
fprintf(o, "<meta name=\"generator\" content=\"AG_genapp\">\n", name);
fprintf(o, "</head>\n<body>\n");

fprintf(o, "<h2>%s parameters</h2>\n<br>\n", name);
/*
fprintf(o, "<table style=\"text-align: left;\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n<tbody>\n<tr>\n");
fprintf(o, "<td style=\"vertical-align: top;\"><span style=\"font-weight: bold;\">Name</span></td>\n");
fprintf(o, "<td style=\"vertical-align: top;\"><span style=\"font-weight: bold;\">Type</span></td>\n");
fprintf(o, "<td style=\"vertical-align: top;\"><span style=\"font-weight: bold;\">Description</span></td>\n</tr>\n");
*/
fprintf(o, "<table style=\"text-align: left;\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n");
fprintf(o, "<tr>\n<td><b>Name</b></td>\n");
fprintf(o, "<td><b>Type</b></td>\n");
fprintf(o, "<td><b>Description</b></td>\n</tr>\n");

int count = parArr.size();
for (int i=0; i<count; ++i) {
	const PilDesc& pilDesc = parArr[i];
	const char* typeName = "";
	switch (pilDesc.type) {
		case 'b': typeName = "Boolean"; break;
		case 'i': typeName = "Integer"; break;
		case 'r': typeName = "Real"; break;
		case 's': typeName = "String";
		}
	fprintf(o, "<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n", pilDesc.identifier.c_str(), typeName, pilDesc.prompt.c_str());
	}
fprintf(o, "</table>\n</body>\n</html>\n");
}



int main(int argC, char *argV[])
{
if (argC!=2) {
	printf("Syntax:\n   AG_genapp <name.par>\nto generate name.cpp source code and name.html documentation\n");
	printf("Existing files won't be overwritten\n");
	return 0;
	}
const char* parfile = argV[1];
char name[1024];
strcpy(name, parfile);
int len = strlen(name);
if (strcmp(name+len-4, ".par")==0)
	name[len-4] = 0;
char cppName[1024];
char htmlName[1024];
strcpy(cppName, name);
strcat(cppName, ".cpp");
strcpy(htmlName, name);
strcat(htmlName, ".html");

/// Write your app code here
FILE* i = fopen(parfile, "rt");
if (!i) {
	fprintf(stderr, "%s not found\n", parfile);
	return -1;
	}

/// Reading the .par file
PilDescArr* parArr = ReadParams(i);
fclose(i);

FILE* o = fopen(cppName, "rt");
if (o) {
	fclose(o);
	printf("Skipping %s because already existing\n", cppName);
	}
else {
	FILE* o = fopen(cppName, "wt");
	if (o) {
		printf("Generating %s source file\n", cppName);
		WriteSource(*parArr, o);
		fclose(o);
		}
	else
		fprintf(stderr, "Could not write to %s\n", cppName);
	}

o = fopen(htmlName, "rt");
if (o) {
	fclose(o);
	printf("Skipping %s because already existing\n", htmlName);
	}
else {
	FILE* o = fopen(htmlName, "wt");
	if (o) {
		printf("Generating %s document\n", htmlName);
		WriteHtml(name, *parArr, o);
		fclose(o);
		}
	else
		fprintf(stderr, "Could not write to %s\n", htmlName);
	}

delete parArr;

return 0;
}
