#include <string>
#include <vector>
#include "utils.h"

bool CommandOptions::parseCommandLine(int argc, char **argv) {
	file = "obj/cube.obj";
	FFfile = "c:/temp/FF.csv";
	cr.set(0.0, 0.0, 1.0); // blue
	specular.set(1.0,1.0,1.0,1.0);
	e.set(2.0,2.0,2.0); // eye coordinates
	at.set(0.0,0.0,0.0); // look at coordinates
	l.set(3.0,3.0,3.0); // point light coordinates
	normMul = 1.0;
	argc = argc;
	argv = argv;
	bool returnOptions = true;
    for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "-help")) {
			returnOptions = false;
			break;
		} else if ((arg == "-tc") || (arg == "-testcase")) {
			if (i + 1 < argc) {
				testCase = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-nd") || (arg == "-nodepth")) {
			depthTest = false;
		} else if ((arg == "-debug")) {
			debug = true;
		} else if ((arg == "-shadows")) {
			shadows = true;
		} else if ((arg == "-shadowsOnly")) {
			shadowsOnly = true;
		} else if ((arg == "-rndcolor")) {
			rndcolor = true;
		} else if ((arg == "-smoothcolor")) {
			smoothcolor = true;
		} else if ((arg == "-occlusion")) {
			occlusion = true;
		} else if ((arg == "-patchesOnly")) {
			patchesOnly = true;
		} else if ((arg == "-loadFF")) {
			loadFF = true;
		} else if ((arg == "-addPRAmbient")) {
			addPRAmbient = true;
		} else if ((arg == "-sm") || (arg == "-shademodel")) {
			if (i + 1 < argc) {
				shadeModel = OSU_FLAT;
				int tmp = atoi(argv[++i]);
				if (tmp == OSU_SMOOTH) shadeModel = OSU_SMOOTH;
			} else {
				printf("%s option requires one argument (0|1)\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-s") || (arg == "-scale")) {
			if (i + 1 < argc) {
				scale = atof(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-rd") || (arg == "-reflectionDepth")) {
			if (i + 1 < argc) {
				reflectionDepth = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-iters")) {
			if (i + 1 < argc) {
				iters = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-nparts")) {
			if (i + 1 < argc) {
				nparts = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-ncontour")) {
			if (i + 1 < argc) {
				ncontour = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-scene")) {
			if (i + 1 < argc) {
				scene = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-solver")) {
			if (i + 1 < argc) {
				solver = atoi(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-normMul")) {
			if (i + 1 < argc) {
				normMul = atof(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-ambient")) {
			if (i + 1 < argc) {
				ambient = (float)atof(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-intensity")) {
			if (i + 1 < argc) {
				intensity = (float)atof(argv[++i]);
			} else {
				printf("%s option requires one argument.\n", arg);
				returnOptions = false;
				break;
			}
		} else if ((arg == "-f") || arg == "-file") {
			if (i + 1 < argc) {
				file = argv[++i];
			} else {
				printf("-file option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if (arg == "-FFfile") {
			if (i + 1 < argc) {
				FFfile = argv[++i];
			} else {
				printf("-FFfile option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-diffuse")) {
			if (i + 1 < argc) {
				std::string diffusestr = argv[++i];
				if (!cr.parse(diffusestr)) {
					std::cout << "Invalid color value " << diffusestr << std::endl;
				}
			} else {
				printf("-diffuse option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-specular")) {
			if (i + 1 < argc) {
				std::string specularstr = argv[++i];
				if (!specular.parse(specularstr)) {
					std::cout << "Invalid color value " << specularstr << std::endl;
				}
			} else {
				printf("-specular option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-e")) {
			if (i + 1 < argc) {
				std::string estr = argv[++i];
				if (!e.parse(estr)) {
					std::cout << "Invalid eye coords " << estr << std::endl;
				}
			} else {
				printf("-e option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-at")) {
			if (i + 1 < argc) {
				std::string atstr = argv[++i];
				if (!at.parse(atstr)) {
					std::cout << "Invalid look 'at' coords " << atstr << std::endl;
				}
			} else {
				printf("-at option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-l")) {
			if (i + 1 < argc) {
				std::string estr = argv[++i];
				if (!l.parse(estr)) {
					std::cout << "Invalid light coords " << estr << std::endl;
				}
			} else {
				printf("-l option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else {
			printf("WARN: unrecognized argument '%s' will be ignored.\n",arg);
		}
	}
	if (!shadows && shadowsOnly) {
		std::cout << "WARN: Shadows disabled. Use -shadows with -shadowsOnly" << std::endl;
		shadowsOnly = false;
	}
	if (!(solver == 0 || solver == 1)) {
		std::cout << "solver must be 0 (Gauss-Seidel) or 1 (progressive)" << std::endl;
		returnOptions = false;
	}
	return returnOptions;
}

